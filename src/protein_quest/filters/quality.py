import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

from cyclopts.types import StdioPath
from tqdm.rich import tqdm

from protein_quest.clustering import cluster_structures, structure_sort_key
from protein_quest.pdbe.ws import Scores
from protein_quest.structure.files import LocateStructureFilesByIdResult
from protein_quest.structure.formats import read_structure
from protein_quest.structure.metadata import structure_metadata

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class QualityStructure:
    """Adapter that combines StructureMetadata with geometry_quality for ClusterableStructure.

    Implements the
    [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol
    so that structures filtered by PDBe quality can participate in residue-range
    clustering.

    Note: geometry_quality is guaranteed to be non-None as structures
    with None quality are filtered out earlier in the pipeline.
    """

    id: str
    uniprot_start: int
    uniprot_end: int
    resolution_value: float
    sequence_identity: float
    chain_length: int
    geometry_quality: float

    def __hash__(self) -> int:
        return hash((self.id, self.uniprot_start, self.uniprot_end))


@dataclass(frozen=True, slots=True)
class UnclusteredStructure:
    """Structure without UniProt accession that cannot be clustered.

    Contains the minimum information needed to apply quality filtering
    to structures that lack UniProt metadata.

    Note: geometry_quality is guaranteed to be non-None as structures
    with None quality are filtered out earlier in the pipeline.
    """

    input_file: Path
    pdb_id: str
    geometry_quality: float


@dataclass(frozen=True, slots=True)
class FilterQualityResult:
    """Result of filtering a structure file based on PDBe quality scores.

    Attributes:
        pdb_id: The PDB ID of the structure if available.
        input_file: The path to the input structure file if available.
        geometry_quality: The geometry quality score if available.
        passed: A boolean indicating whether the structure passed the filter.
        reason: The reason for discarding or passing the structure.
            None if passed without any special reason.
    """

    pdb_id: str | None = None
    input_file: Path | None = None
    geometry_quality: float | None = None
    passed: bool = False
    reason: str | None = None


class PdbIdGeometryQualityPair(TypedDict):
    pdb_id: str
    geometry_quality: float | None


def _associate_files_with_sorted_scores(
    scores: dict[str, Scores], located_ids: LocateStructureFilesByIdResult, top: int | None = None
) -> tuple[list[PdbIdGeometryQualityPair], dict[str, list[Path]]]:
    flattened_scores: list[PdbIdGeometryQualityPair] = []
    found_ids = {pdb_id for pdb_id, _ in located_ids.found}
    for pdb_id, score in scores.items():
        if pdb_id not in found_ids:
            continue
        flattened_scores.append(
            PdbIdGeometryQualityPair(
                pdb_id=pdb_id,
                geometry_quality=score.geometry_quality,
            )
        )

    def sorter(x: PdbIdGeometryQualityPair) -> tuple[bool, float | None]:
        # None as worst and highest as best with reverse sort
        return (x["geometry_quality"] is not None, x["geometry_quality"])

    sorted_scores = sorted(flattened_scores, key=sorter, reverse=True)
    top_sorted_scores = sorted_scores[:top] if top is not None else sorted_scores

    # Convert (('a', file_a), ('b', file_b), ('a', file_c), ('b', file_c)) to {'a': [file_a, file_c], 'b': [file_b]}
    found_id2files: dict[str, list[Path]] = {}
    found_files2ids: dict[Path, str] = {}
    for pdb_id, file in located_ids.found:
        if file in found_files2ids:
            logger.info(
                f"File {file} already associated with PDB ID {found_files2ids[file]}, "
                f"skipping association with PDB ID {pdb_id}"
            )
            continue
        found_id2files.setdefault(pdb_id, []).append(file)
        found_files2ids[file] = pdb_id
    return top_sorted_scores, found_id2files


def filter_by_pdbe_quality(
    scores: dict[str, Scores],
    located_ids: LocateStructureFilesByIdResult,
    /,
    *,
    minimal_geometry_quality: float = 0.0,
    top: int | None = None,
    pass_given_resolution: bool = False,
) -> list[FilterQualityResult]:
    """Filter structure files based on PDBe quality scores.

    Args:
        scores: A dictionary mapping PDB IDs to their corresponding Scores objects.
        located_ids: A dictionary containing located structure files by PDB ID.
        minimal_geometry_quality: Minimum geometry quality score to pass the filter.
        top: If set, only consider the top N structures based on geometry quality.
        pass_given_resolution: If set, structures with a valid resolution will pass regardless of other criteria.

    Returns:
        A list of QualityResult objects.
        No files are written to disk; the caller is responsible for copying or moving the files as needed.
    """
    top_sorted_scores, found_id2files = _associate_files_with_sorted_scores(scores, located_ids, top)

    results: list[FilterQualityResult] = []
    for score_dict in tqdm(top_sorted_scores, desc="Filtering on PDBe quality", unit="file"):
        pdb_id = score_dict["pdb_id"]
        score = scores[pdb_id]
        input_files = found_id2files[pdb_id]
        for input_file in input_files:
            if pass_given_resolution and (
                # Only read structure when needed as checking score is cheaper than reading structure
                score.geometry_quality is None or score.geometry_quality < minimal_geometry_quality
            ):
                structure = read_structure(input_file)
                if structure.resolution != 0.0:
                    results.append(
                        FilterQualityResult(
                            pdb_id=pdb_id,
                            input_file=input_file,
                            geometry_quality=score.geometry_quality,
                            passed=True,
                            reason=f"Passed due to valid resolution {structure.resolution}",
                        )
                    )
                    logger.debug(f"Passing {pdb_id} due to valid resolution {structure.resolution}")
                    continue
            if score.geometry_quality is None:
                results.append(
                    FilterQualityResult(
                        pdb_id=pdb_id,
                        input_file=input_file,
                        reason="No geometry quality score",
                    )
                )
                continue
            if score.geometry_quality >= minimal_geometry_quality:
                results.append(
                    FilterQualityResult(
                        pdb_id=pdb_id,
                        input_file=input_file,
                        geometry_quality=score.geometry_quality,
                        passed=True,
                    )
                )
            else:
                results.append(
                    FilterQualityResult(
                        pdb_id=pdb_id,
                        input_file=input_file,
                        geometry_quality=score.geometry_quality,
                        reason=f"Geometry quality score {score.geometry_quality} < {minimal_geometry_quality}",
                    )
                )

    results.extend(
        FilterQualityResult(
            pdb_id=pdb_id,
            geometry_quality=scores[pdb_id].geometry_quality,
            reason="File not found",
        )
        for pdb_id in located_ids.not_found
    )
    results.extend(
        FilterQualityResult(
            input_file=extra,
            reason="File not found in quality scores",
        )
        for extra in located_ids.extras
    )

    return results


def _extract_uniprot_metadata(file: Path) -> tuple[str, int, int, int, float] | None:
    """Extract UniProt accession and residue range from a structure file.

    Returns a tuple of (uniprot_accession, uniprot_start, uniprot_end, chain_length, resolution) or
    None if no UniProt accession is found or parsing fails.
    """
    try:
        structure = read_structure(file)
        metadata = structure_metadata(structure, path=file)
    except Exception:  # noqa: BLE001 - skip structures that fail to parse
        return None
    if metadata.uniprot_accession is None:
        return None
    return (
        metadata.uniprot_accession,
        metadata.uniprot_start,
        metadata.uniprot_end,
        metadata.chain_length,
        metadata.resolution,
    )


def _build_selected_results(
    selected_structures: set[QualityStructure],
    minimal_geometry_quality: float,
) -> list[FilterQualityResult]:
    """Build filter results for selected clustered structures.

    Assumes all structures have non-None geometry_quality.
    """
    results: list[FilterQualityResult] = []
    for qs in selected_structures:
        if qs.geometry_quality >= minimal_geometry_quality:
            results.append(
                FilterQualityResult(
                    pdb_id=qs.id,
                    geometry_quality=qs.geometry_quality,
                    passed=True,
                )
            )
        else:
            results.append(
                FilterQualityResult(
                    pdb_id=qs.id,
                    geometry_quality=qs.geometry_quality,
                    passed=False,
                    reason=f"Geometry quality score {qs.geometry_quality} < {minimal_geometry_quality}",
                )
            )
    return results


@dataclass(frozen=True, slots=True)
class PrepareQualityStructuresResult:
    """Result of preparing quality structures from scores and located files.

    Attributes:
        clusterable_structures: List of QualityStructure objects with UniProt metadata
            and non-None geometry_quality
        unclustered_structures: List of UnclusteredStructure objects without UniProt accession
            but with non-None geometry_quality for filtering
        no_quality_results: List of FilterQualityResult objects with None geometry_quality
        resolution_passed_results: List of FilterQualityResult objects that passed due to
            valid resolution when pass_given_resolution is True
    """

    clusterable_structures: list[QualityStructure]
    unclustered_structures: list[UnclusteredStructure]
    no_quality_results: list[FilterQualityResult]
    resolution_passed_results: list[FilterQualityResult]


def prepare_quality_structures(
    scores: dict[str, Scores],
    located_ids: LocateStructureFilesByIdResult,
    /,
    *,
    top: int | None = None,
    pass_given_resolution: bool = False,
) -> PrepareQualityStructuresResult:
    """Prepare QualityStructure objects from scores and located files.

    This function handles all I/O operations by reading structure files to extract
    UniProt metadata. It returns both clustered structures (with UniProt accession)
    and unclustered structures (without UniProt accession) with quality data.
    Structures with None geometry_quality are separated out.

    Args:
        scores: A dictionary mapping PDB IDs to their corresponding Scores objects.
        located_ids: A dictionary containing located structure files by PDB ID.
        top: If set, only consider the top N structures based on geometry quality.
        pass_given_resolution: If set, structures with a valid resolution will pass
            regardless of other criteria.

    Returns:
        A PrepareQualityStructuresResult object containing:
        - clusterable_structures: List of QualityStructure objects with UniProt metadata
            and non-None geometry_quality
        - unclustered_structures: List of UnclusteredStructure objects without UniProt accession
            but with non-None geometry_quality for filtering
        - no_quality_results: List of FilterQualityResult objects with None geometry_quality
        - resolution_passed_results: List of FilterQualityResult objects that passed due to
            valid resolution when pass_given_resolution is True
    """
    top_sorted_scores, found_id2files = _associate_files_with_sorted_scores(scores, located_ids, top)

    clusterable_structures: list[QualityStructure] = []
    unclustered_structures: list[UnclusteredStructure] = []
    no_quality_results: list[FilterQualityResult] = []
    resolution_passed_results: list[FilterQualityResult] = []

    for score_dict in tqdm(top_sorted_scores, desc="Building clusters from PDBe quality", unit="file"):
        pdb_id = score_dict["pdb_id"]
        score = scores[pdb_id]
        input_files = found_id2files[pdb_id]

        # Early filtering: skip structures with None geometry_quality
        if score.geometry_quality is None:
            # Check pass_given_resolution when geometry_quality is None
            if pass_given_resolution:
                for input_file in input_files:
                    structure = read_structure(input_file)
                    if structure.resolution != 0.0:
                        resolution_passed_results.append(
                            FilterQualityResult(
                                pdb_id=pdb_id,
                                input_file=input_file,
                                geometry_quality=score.geometry_quality,
                                passed=True,
                                reason=f"Passed due to valid resolution {structure.resolution}",
                            )
                        )
                        logger.debug(f"Passing {pdb_id} due to valid resolution {structure.resolution}")
                        continue
            no_quality_results.extend(
                FilterQualityResult(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=None,
                    passed=False,
                    reason="No geometry quality score",
                )
                for input_file in input_files
            )
            continue

        for input_file in input_files:
            # Check pass_given_resolution for all structures when enabled
            if pass_given_resolution:
                structure = read_structure(input_file)
                if structure.resolution != 0.0:
                    resolution_passed_results.append(
                        FilterQualityResult(
                            pdb_id=pdb_id,
                            input_file=input_file,
                            geometry_quality=score.geometry_quality,
                            passed=True,
                            reason=f"Passed due to valid resolution {structure.resolution}",
                        )
                    )
                    logger.debug(f"Passing {pdb_id} due to valid resolution {structure.resolution}")
                    continue

            meta = _extract_uniprot_metadata(input_file)
            if meta is None:
                # No UniProt accession found — create UnclusteredStructure with quality data
                unclustered_structures.append(
                    UnclusteredStructure(
                        input_file=input_file,
                        pdb_id=pdb_id,
                        geometry_quality=score.geometry_quality,
                    )
                )
                continue
            uniprot_accession, uniprot_start, uniprot_end, chain_length, resolution = meta
            clusterable_structures.append(
                QualityStructure(
                    id=f"{pdb_id}:{uniprot_accession}",
                    uniprot_start=uniprot_start,
                    uniprot_end=uniprot_end,
                    resolution_value=resolution,
                    sequence_identity=0.0,
                    chain_length=chain_length,
                    geometry_quality=score.geometry_quality,
                )
            )

    return PrepareQualityStructuresResult(
        clusterable_structures=clusterable_structures,
        unclustered_structures=unclustered_structures,
        no_quality_results=no_quality_results,
        resolution_passed_results=resolution_passed_results,
    )


def process_quality_clusters(
    clusterable_structures: list[QualityStructure],
    cluster_by_uniprot_accession_and_coverage: int,
) -> tuple[set[QualityStructure], set[QualityStructure]]:
    """Process clusterable structures and select top N per UniProt cluster.

    This function performs pure computation without any I/O operations. It groups
    structures by UniProt accession, clusters them by residue-range overlap, and
    selects the top N structures per cluster.

    Args:
        clusterable_structures: List of QualityStructure objects to cluster.
        cluster_by_uniprot_accession_and_coverage: Number of top structures to keep
            per UniProt cluster.

    Returns:
        A tuple of (selected_structures, not_selected_structures) where:
        - selected_structures: Set of QualityStructure objects selected as top N per cluster
        - not_selected_structures: Set of QualityStructure objects not selected
    """
    # Group by UniProt accession first
    accession_groups: dict[str, list[QualityStructure]] = {}
    for qs in clusterable_structures:
        # Extract UniProt accession from id (format: "PDB_ID:UNIPROT_ACC")
        uniprot_accession = qs.id.split(":")[-1]
        accession_groups.setdefault(uniprot_accession, []).append(qs)

    # Cluster by residue-range overlap within each UniProt group, select top N per cluster
    all_clusters: list[list[QualityStructure]] = []
    for structures in accession_groups.values():
        if structures:
            clusters = cluster_structures(structures)
            all_clusters.extend(clusters)

    selected_structures: set[QualityStructure] = set()
    for cluster in all_clusters:
        sorted_cluster = sorted(cluster, key=structure_sort_key)
        for qs in sorted_cluster[:cluster_by_uniprot_accession_and_coverage]:
            selected_structures.add(qs)

    not_selected = set(clusterable_structures) - selected_structures
    return selected_structures, not_selected


def filter_by_pdbe_quality_clustered(
    scores: dict[str, Scores],
    located_ids: LocateStructureFilesByIdResult,
    *,
    minimal_geometry_quality: float = 0.0,
    top: int | None = None,
    pass_given_resolution: bool = False,
    cluster_by_uniprot_accession_and_coverage: int = 1,
) -> list[FilterQualityResult]:
    """Filter structure files by PDBe quality with UniProt-based clustering.

    Structures are grouped by UniProt accession, then clustered by residue-range
    overlap using hierarchical clustering. The top N structures per cluster are
    selected using the geometry-quality-aware sort key.

    This is opt-in via the ``cluster_by_uniprot_accession_and_coverage`` parameter.
    When not set (or 0), falls back to the standard
    [filter_by_pdbe_quality][protein_quest.filters.quality.filter_by_pdbe_quality] path.

    Args:
        scores: A dictionary mapping PDB IDs to their corresponding Scores objects.
        located_ids: A dictionary containing located structure files by PDB ID.
        minimal_geometry_quality: Minimum geometry quality score to pass the filter.
        top: If set, only consider the top N structures based on geometry quality.
        pass_given_resolution: If set, structures with a valid resolution will pass
            regardless of other criteria.
        cluster_by_uniprot_accession_and_coverage: Number of top structures to keep
            per UniProt cluster. If 0, clustering is disabled.

    Returns:
        A list of FilterQualityResult objects.
    """
    if cluster_by_uniprot_accession_and_coverage == 0:
        msg = (
            "Clustering is disabled. Use filter_by_pdbe_quality() instead of "
            "filter_by_pdbe_quality_clustered() when "
            "cluster_by_uniprot_accession_and_coverage is 0."
        )
        raise ValueError(msg)

    # Prepare quality structures (I/O operations)
    quality_result = prepare_quality_structures(
        scores, located_ids, top=len(scores), pass_given_resolution=pass_given_resolution
    )

    results: list[FilterQualityResult] = quality_result.no_quality_results
    results.extend(quality_result.resolution_passed_results)

    # Process clusters (pure computation) - now guaranteed to have non-None geometry_quality
    selected_structures, not_selected = process_quality_clusters(
        quality_result.clusterable_structures, cluster_by_uniprot_accession_and_coverage
    )

    # Build results for selected clustered structures
    results.extend(_build_selected_results(selected_structures, minimal_geometry_quality))

    # Structures not in top N of their cluster get a "not selected" reason
    results.extend(
        FilterQualityResult(
            pdb_id=qs.id,
            geometry_quality=qs.geometry_quality,
            passed=False,
            reason=f"Not selected in UniProt cluster (top {cluster_by_uniprot_accession_and_coverage})",
        )
        for qs in not_selected
    )

    # Unclustered structures (no UniProt accession): apply quality filtering
    # Now guaranteed to have non-None geometry_quality
    nr_from_clusters = len(selected_structures)
    places_left = top - nr_from_clusters if top is not None else len(quality_result.unclustered_structures)
    for us in quality_result.unclustered_structures:
        if places_left <= 0:
            results.append(
                FilterQualityResult(
                    pdb_id=us.pdb_id,
                    input_file=us.input_file,
                    geometry_quality=us.geometry_quality,
                    passed=False,
                    reason=f"Not selected due to top {top} limit reached",
                )
            )
            break
        if us.geometry_quality >= minimal_geometry_quality:
            # Meets quality threshold
            results.append(
                FilterQualityResult(
                    pdb_id=us.pdb_id,
                    input_file=us.input_file,
                    geometry_quality=us.geometry_quality,
                    passed=True,
                    reason="No UniProt accession but meets quality threshold",
                )
            )
        else:
            # Fails quality threshold
            results.append(
                FilterQualityResult(
                    pdb_id=us.pdb_id,
                    input_file=us.input_file,
                    geometry_quality=us.geometry_quality,
                    passed=False,
                    reason=(
                        f"No UniProt accession and geometry quality {us.geometry_quality} < {minimal_geometry_quality}"
                    ),
                )
            )

    # Not-found and extras
    results.extend(
        FilterQualityResult(
            pdb_id=pdb_id,
            geometry_quality=scores[pdb_id].geometry_quality,
            reason="File not found",
        )
        for pdb_id in located_ids.not_found
    )
    results.extend(
        FilterQualityResult(
            input_file=extra,
            reason="File not found in quality scores",
        )
        for extra in located_ids.extras
    )

    return results


def write_quality_stats_csv(
    results: list[FilterQualityResult],
    write_stats: StdioPath,
    output_dir: Path,
):
    """Writes a CSV file containing quality statistics for PDB IDs.

    Args:
        results: A list of QualityResult objects containing the quality results for each PDB ID.
        write_stats: Path to the output CSV file to write statistics.
        output_dir: The directory where the output files are located.

    The CSV file will contain the following columns:
        - pdb_id: The PDB ID of the structure.
        - input_file: The path to the input structure file.
        - geometry_quality: The geometry quality score from the Scores object.
        - passed: A boolean indicating whether the structure passed the filter.
        - output_file: The path to the output structure file if it passed the filter.
        - reason: The reason for discarding or passing the structure.
    """
    with write_stats.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["pdb_id", "input_file", "geometry_quality", "passed", "output_file", "reason"]
        )
        writer.writeheader()
        for result in results:
            output_file = ""
            if result.passed and result.input_file is not None:
                output_file = output_dir / result.input_file.name

            writer.writerow(
                {
                    "pdb_id": result.pdb_id,
                    "input_file": result.input_file,
                    "geometry_quality": result.geometry_quality,
                    "passed": result.passed,
                    "output_file": output_file,
                    "reason": result.reason,
                }
            )
