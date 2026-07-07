import csv
import logging
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

from cyclopts.types import StdioPath
from tqdm.rich import tqdm

from protein_quest.clustering import cluster_structures
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
    """

    id: str
    uniprot_accession: str
    uniprot_start: int
    uniprot_end: int
    sequence_identity: float
    chain_length: int
    geometry_quality: float
    input_file: Path

    def __hash__(self) -> int:
        return hash((self.id, self.uniprot_start, self.uniprot_end))

    @property
    def resolution_value(self) -> float:
        """Return ignored resolultion so cluster member sorting is done on geometry quality

        See [structure_sort_key][protein_quest.clustering.structure_sort_key] for sorting logic.
        """
        return 0.0


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


def _quality_threshold_result(
    *,
    pdb_id: str,
    input_file: Path,
    geometry_quality: float,
    minimal_geometry_quality: float,
    passed_reason: str | None = None,
    failed_reason_prefix: str = "Geometry quality score",
) -> FilterQualityResult:
    passed = geometry_quality >= minimal_geometry_quality
    if passed:
        return FilterQualityResult(
            pdb_id=pdb_id,
            input_file=input_file,
            geometry_quality=geometry_quality,
            passed=True,
            reason=passed_reason,
        )

    return FilterQualityResult(
        pdb_id=pdb_id,
        input_file=input_file,
        geometry_quality=geometry_quality,
        passed=False,
        reason=f"{failed_reason_prefix} {geometry_quality} < {minimal_geometry_quality}",
    )


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
            results.append(
                _quality_threshold_result(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=score.geometry_quality,
                    minimal_geometry_quality=minimal_geometry_quality,
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


@dataclass(frozen=True, slots=True)
class QualityClusteringPartitions:
    """Partitioned structures for clustered PDBe quality filtering.

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


def partition_structures_for_quality_clustering(
    input_files: Iterable[Path],
    scores: dict[str, Scores],
    /,
    *,
    pass_given_resolution: bool = False,
) -> QualityClusteringPartitions:
    """Partition structures for clustered PDBe quality filtering.

    This function handles all I/O operations by reading structure files to extract
    UniProt metadata. It returns both clustered structures (with UniProt accession)
    and unclustered structures (without UniProt accession) with quality data.
    Structures with None geometry_quality are separated out.

    Each input file's PDB ID is obtained from the structure's ``name`` field (the ``id``
    field in the CIF file). Files whose structure name is not found in ``scores`` are
    silently skipped.

    Args:
        input_files: Iterable of structure file paths.
        scores: A dictionary mapping PDB IDs to their corresponding Scores objects.
        pass_given_resolution: If set, structures with a valid resolution will pass
            regardless of other criteria.

    Returns:
        A QualityClusteringPartitions object containing:
        - clusterable_structures: List of QualityStructure objects with UniProt metadata
            and non-None geometry_quality
        - unclustered_structures: List of UnclusteredStructure objects without UniProt accession
            but with non-None geometry_quality for filtering
        - no_quality_results: List of FilterQualityResult objects with None geometry_quality
            or failed to read structure or extract metadata
        - resolution_passed_results: List of FilterQualityResult objects that passed due to
            valid resolution when pass_given_resolution is True
    """
    clusterable_structures: list[QualityStructure] = []
    unclustered_structures: list[UnclusteredStructure] = []
    no_quality_results: list[FilterQualityResult] = []
    resolution_passed_results: list[FilterQualityResult] = []

    for input_file in tqdm(input_files, desc="Building clusters from PDBe quality", unit="file"):
        try:
            structure = read_structure(input_file)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"Failed to read structure {input_file}: {e}")
            no_quality_results.append(
                FilterQualityResult(
                    input_file=input_file,
                    passed=False,
                    reason=f"Failed to read structure: {e}",
                )
            )
            continue
        pdb_id = structure.name
        lowered_pdb_id = pdb_id.lower()
        geometry_quality = scores[lowered_pdb_id].geometry_quality if lowered_pdb_id in scores else None

        if pass_given_resolution and structure.resolution != 0.0:
            resolution_passed_results.append(
                FilterQualityResult(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=geometry_quality,
                    passed=True,
                    reason=f"Passed due to valid resolution {structure.resolution}",
                )
            )
            continue

        try:
            metadata = structure_metadata(structure, path=input_file)
        except Exception as e:  # noqa: BLE001
            logger.warning(f"Failed to extract metadata from {input_file}: {e}")
            no_quality_results.append(
                FilterQualityResult(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=geometry_quality,
                    passed=False,
                    reason=f"Failed to extract metadata: {e}",
                )
            )
            continue
        if metadata.is_alphafold and geometry_quality is None:
            resolution_passed_results.append(
                FilterQualityResult(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=geometry_quality,
                    passed=True,
                    reason="AlphaFold structure passes quality filter",
                )
            )
            continue
        if geometry_quality is None:
            no_quality_results.append(
                FilterQualityResult(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=None,
                    passed=False,
                    reason="Missing geometry quality",
                )
            )
            continue
        if metadata.uniprot_accession is None:
            unclustered_structures.append(
                UnclusteredStructure(
                    input_file=input_file,
                    pdb_id=pdb_id,
                    geometry_quality=geometry_quality,
                )
            )
            continue
        clusterable_structures.append(
            QualityStructure(
                id=pdb_id,
                input_file=input_file,
                uniprot_accession=metadata.uniprot_accession,
                uniprot_start=metadata.uniprot_start,
                uniprot_end=metadata.uniprot_end,
                sequence_identity=metadata.sequence_identity,
                chain_length=metadata.chain_length,
                geometry_quality=geometry_quality,
            )
        )

    return QualityClusteringPartitions(
        clusterable_structures=clusterable_structures,
        unclustered_structures=unclustered_structures,
        no_quality_results=no_quality_results,
        resolution_passed_results=resolution_passed_results,
    )


def cluster_and_select_quality_structures(
    clusterable_structures: list[QualityStructure],
    top_per_cluster: int,
    minimal_geometry_quality: float,
) -> list[FilterQualityResult]:
    """Cluster structures by UniProt coverage and select passing results.

    This function performs pure computation without any I/O operations. It groups
    structures by UniProt accession, clusters them by residue-range overlap, and
    marks structures as passed when they are either in the top N of their cluster
    or meet the geometry quality threshold.

    Args:
        clusterable_structures: List of QualityStructure objects to cluster.
        top_per_cluster: Maximum number of top structures to keep per UniProt cluster
            before applying the geometry quality threshold fallback.
        minimal_geometry_quality: Minimum geometry quality threshold for selection.

    Returns:
        A list of FilterQualityResult objects indicating the selection status of each structure.
    """
    accession_groups: dict[str, list[QualityStructure]] = {}
    for qs in clusterable_structures:
        uniprot_accession = qs.uniprot_accession
        accession_groups.setdefault(uniprot_accession, []).append(qs)

    all_clusters: list[list[QualityStructure]] = []
    for structures in accession_groups.values():
        clusters = cluster_structures(structures)
        all_clusters.extend(clusters)

    results: list[FilterQualityResult] = []
    for cluster in all_clusters:
        for i, qs in enumerate(cluster):
            in_top = i < top_per_cluster
            ok_quality = qs.geometry_quality >= minimal_geometry_quality
            if in_top and ok_quality:
                results.append(
                    FilterQualityResult(
                        pdb_id=qs.id,
                        input_file=qs.input_file,
                        geometry_quality=qs.geometry_quality,
                        passed=True,
                    )
                )
            else:
                results.append(
                    FilterQualityResult(
                        pdb_id=qs.id,
                        input_file=qs.input_file,
                        geometry_quality=qs.geometry_quality,
                        passed=False,
                        reason=(
                            "Not selected in UniProt cluster (top "
                            f"{top_per_cluster}) or geometry quality < {minimal_geometry_quality}"
                        ),
                    )
                )
    return results


def filter_by_pdbe_quality_clustered(
    scores: dict[str, Scores],
    input_files: Iterable[Path],
    /,
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
        input_files: Iterable of structure file paths (e.g. from
            [glob_structure_files][protein_quest.structure.files.glob_structure_files]).
        minimal_geometry_quality: Minimum geometry quality score to pass the filter.
        top: Maximum number of unclustered structures (those without a UniProt accession)
            that can pass the filter. Structures are ranked by geometry quality descending
            and only the top N are allowed to pass. If None, all qualifying structures pass.
        pass_given_resolution: If set, structures with a valid resolution will pass
            regardless of other criteria.
        cluster_by_uniprot_accession_and_coverage: Number of top structures to keep
            per UniProt cluster. If 0, clustering is disabled.

    Returns:
        A list of FilterQualityResult objects.
    """
    partitions = partition_structures_for_quality_clustering(
        input_files,
        scores,
        pass_given_resolution=pass_given_resolution,
    )

    results = cluster_and_select_quality_structures(
        partitions.clusterable_structures,
        top_per_cluster=cluster_by_uniprot_accession_and_coverage,
        minimal_geometry_quality=minimal_geometry_quality,
    )

    sorted_unclustered = sorted(partitions.unclustered_structures, key=lambda us: us.geometry_quality, reverse=True)
    passed_unclustered = 0
    for us in sorted_unclustered:
        result = _quality_threshold_result(
            pdb_id=us.pdb_id,
            input_file=us.input_file,
            geometry_quality=us.geometry_quality,
            minimal_geometry_quality=minimal_geometry_quality,
            passed_reason="No UniProt accession but meets quality threshold",
            failed_reason_prefix="No UniProt accession and geometry quality",
        )
        if result.passed and top is not None and passed_unclustered >= top:
            result = FilterQualityResult(
                pdb_id=us.pdb_id,
                input_file=us.input_file,
                geometry_quality=us.geometry_quality,
                passed=False,
                reason=f"Excluded by top {top} limit for unclustered structures: {result.reason}",
            )
        if result.passed:
            passed_unclustered += 1
        results.append(result)

    results.extend(partitions.resolution_passed_results)
    results.extend(partitions.no_quality_results)
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
