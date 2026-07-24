import csv
import logging
from collections.abc import Iterable
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Annotated

import gemmi
from cyclopts import Parameter
from cyclopts.types import NonNegativeInt, NormFloat, PositiveInt, StdioPath
from cyclopts.validators import Number

from protein_quest.alphafold.confidence import ConfidenceFilterQuery, filter_structure_on_confidence
from protein_quest.clustering import cluster_structures, sort_structures
from protein_quest.parallel import SchedulerAddress, map_with_progress
from protein_quest.pdbe.ws import Scores
from protein_quest.structure.formats import read_structure
from protein_quest.structure.metadata import StructureMetadata, structure_metadata
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


@Parameter(name="*")
@dataclass(frozen=True, slots=True)
class CombinedFilterQuery:
    """Query object to apply combined filtering.

    Parameters:
        min_confidence: Minimal confidence (plDDT) for AlphaFold structures to pass the filter.
        min_residues: Min residues in chain A.
        max_residues: Max residues in chain A.
        min_geometry_quality: Minimal geometry quality score to pass the filter.
        min_sequence_identity: Minimum sequence identity ratio to the Uniprot sequence for a structure to be passed.
            If not set then discards structures that are not fully identical to the Uniprot sequence.
            For example if set to 0.8 then structures that have sequence identity below 0.8 are discarded.
        top_uniprot_cluster: Maximum number of files to keep for structures per cluster per Uniprot accession.
            Alphafold structures are excluded from this limit.
        top_non_uniprot: Maximum number of files to keep for structures without Uniprot accession.

    """

    min_confidence: Annotated[float, Parameter(validator=(Number(lte=100, gte=0)))] = 70.0
    min_residues: NonNegativeInt = 0
    max_residues: PositiveInt = 10_000_000
    min_geometry_quality: float = 50.0
    min_sequence_identity: NormFloat = 1.0
    top_uniprot_cluster: NonNegativeInt = 1_000
    top_non_uniprot: NonNegativeInt = 0

    def confidence_filter_query(self) -> "ConfidenceFilterQuery":
        return ConfidenceFilterQuery(
            confidence=self.min_confidence,
            min_residues=self.min_residues,
            max_residues=self.max_residues,
        )


@dataclass(frozen=True, slots=True)
class CombinedFilterResult:
    """Result of combined filtering.

    Compatible with [SortableStructure][protein_quest.clustering.SortableStructure] protocol.

    Parameters:
        input_file: Path to the input structure file.
        pdb_id: PDB ID of the structure.
        metadata: Structure metadata.
        high_confidence_residues_count: Number of residues with high confidence (plDDT)
            for AlphaFold structures.
        geometry_quality: Geometry quality score.
        passed: Whether the structure passed the filter.
        reason: Reason for failure if the structure did not pass the filter.
        output_file: Path to the output structure file if the structure passed the filter.
    """

    input_file: Path | None = None
    pdb_id: str | None = None
    metadata: StructureMetadata | None = None
    high_confidence_residues_count: int | None = None
    geometry_quality: float | None = None
    passed: bool = False
    reason: str | None = None
    output_file: Path | None = None

    @property
    def id(self):
        if not self.pdb_id:
            msg = "pdb_id is not set"
            raise ValueError(msg)
        return self.pdb_id

    @property
    def resolution_value(self) -> float:
        if self.metadata is None:
            return 0.0
        return self.metadata.resolution

    @property
    def chain_length(self) -> int:
        if self.metadata is None:
            return 0
        return self.metadata.chain_length

    @property
    def sequence_identity(self) -> float:
        if self.metadata is None:
            return 0.0
        return self.metadata.sequence_identity


@dataclass(frozen=True, slots=True)
class ResolutionOrGeometryQualityClusterableStructure:
    """Structure that can be clustered.

    Compatible with [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Used internally by the [combined_filter][protein_quest.filters.combined.combined_filter] function.

    Parameters:
        input_file: Path to the input structure file.
        metadata: Structure metadata.
        geometry_quality: Geometry quality score.

    """

    input_file: Path
    metadata: StructureMetadata
    geometry_quality: float | None = None

    def __hash__(self) -> int:
        return hash(self.input_file)

    @property
    def id(self) -> str:
        return self.metadata.id

    @property
    def uniprot_accession(self) -> str:
        if self.metadata.uniprot_accession is None:
            msg = f"Structure {self.metadata.id} does not have a UniProt accession."
            raise ValueError(msg)
        return self.metadata.uniprot_accession

    @property
    def uniprot_start(self) -> int:
        return self.metadata.uniprot_start

    @property
    def uniprot_end(self) -> int:
        return self.metadata.uniprot_end

    @property
    def sequence_identity(self) -> float:
        return self.metadata.sequence_identity

    @property
    def chain_length(self) -> int:
        return self.metadata.chain_length

    @property
    def resolution_value(self) -> float:
        if self.metadata is None or self.metadata.resolution is None:
            return 0.0
        return self.metadata.resolution


@dataclass(frozen=True, slots=True)
class CombinedPartitions:
    """Partitions of structures and results.

    Used internally by the [combined_filter][protein_quest.filters.combined.combined_filter] function.

    Attributes:
        uniprot_with_resolution: Structures with UniProt accession and resolution.
        uniprot_with_geometry_quality: Structures with UniProt accession and geometry quality.
        resolution_only: Structures with resolution but no UniProt accession.
        geometry_quality_only: Structures with geometry quality but no UniProt accession.
        processed: Results of structures that have been processed.
    """

    uniprot_with_resolution: list[ResolutionOrGeometryQualityClusterableStructure] = field(default_factory=list)
    uniprot_with_geometry_quality: list[ResolutionOrGeometryQualityClusterableStructure] = field(default_factory=list)
    resolution_only: list[CombinedFilterResult] = field(default_factory=list)
    geometry_quality_only: list[CombinedFilterResult] = field(default_factory=list)
    processed: list[CombinedFilterResult] = field(default_factory=list)

    def extend(self, others: Iterable["CombinedPartitions"]) -> "CombinedPartitions":
        for other in others:
            self.uniprot_with_resolution.extend(other.uniprot_with_resolution)
            self.uniprot_with_geometry_quality.extend(other.uniprot_with_geometry_quality)
            self.resolution_only.extend(other.resolution_only)
            self.geometry_quality_only.extend(other.geometry_quality_only)
            self.processed.extend(other.processed)
        return self


def _apply_alphafold_filter(
    input_file: Path,
    structure: gemmi.Structure,
    metadata: StructureMetadata,
    filters: CombinedFilterQuery,
    geometry_quality: float | None,
    output_dir: Path,
    copy_method: CopyMethod = "copy",
) -> CombinedFilterResult:
    # See combined_filter docstring for partitioning logic.
    try:
        raw_result = filter_structure_on_confidence(
            structure=structure,
            file=input_file,
            query=filters.confidence_filter_query(),
            filtered_dir=output_dir,
            copy_method=copy_method,
        )
    except Exception as e:  # noqa: BLE001
        return CombinedFilterResult(
            pdb_id=metadata.id,
            input_file=input_file,
            geometry_quality=geometry_quality,
            metadata=metadata,
            passed=False,
            reason=f"Failed to filter structure on confidence: {e}",
        )

    if raw_result.filtered_file is None:
        return CombinedFilterResult(
            pdb_id=metadata.id,
            input_file=input_file,
            geometry_quality=geometry_quality,
            metadata=metadata,
            high_confidence_residues_count=raw_result.count,
            passed=False,
            reason=f"Low confidence or too few or too many residues ({raw_result.count})",
        )

    return CombinedFilterResult(
        pdb_id=metadata.id,
        input_file=input_file,
        geometry_quality=geometry_quality,
        metadata=metadata,
        high_confidence_residues_count=raw_result.count,
        passed=True,
        output_file=raw_result.filtered_file,
    )


def _apply_non_alphafold_filters(
    input_file: Path,
    pdb_id: str,
    metadata: StructureMetadata,
    geometry_quality: float | None,
    filters: CombinedFilterQuery,
    parts: CombinedPartitions,
) -> CombinedPartitions:
    # See combined_filter docstring for partitioning logic.
    if not (filters.min_residues <= metadata.chain_length <= filters.max_residues):
        rrange = f"[{filters.min_residues}, {filters.max_residues}]"
        parts.processed.append(
            CombinedFilterResult(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                metadata=metadata,
                passed=False,
                reason=f"Chain length {metadata.chain_length} not in range {rrange}",
            )
        )
        return parts

    if metadata.uniprot_accession and metadata.sequence_identity < filters.min_sequence_identity:
        parts.processed.append(
            CombinedFilterResult(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                metadata=metadata,
                passed=False,
                reason=f"Sequence identity {metadata.sequence_identity} < {filters.min_sequence_identity}",
            )
        )
        return parts

    # has uniprot with resolution
    if metadata.uniprot_accession and metadata.resolution != 0.0:
        parts.uniprot_with_resolution.append(
            ResolutionOrGeometryQualityClusterableStructure(input_file=input_file, metadata=metadata)
        )
        return parts

    if metadata.uniprot_accession and geometry_quality is not None:
        parts.uniprot_with_geometry_quality.append(
            ResolutionOrGeometryQualityClusterableStructure(
                input_file=input_file,
                metadata=metadata,
                geometry_quality=geometry_quality,
            )
        )
        return parts

    # has resolution only
    if metadata.resolution != 0.0:
        # TODO add min/max resolution filter? currently sorted and top N selected
        parts.resolution_only.append(
            CombinedFilterResult(
                pdb_id=pdb_id,
                input_file=input_file,
                metadata=metadata,
            )
        )
        return parts

    # has geometry_quality only
    if geometry_quality is not None:
        if geometry_quality < filters.min_geometry_quality:
            parts.processed.append(
                CombinedFilterResult(
                    pdb_id=pdb_id,
                    input_file=input_file,
                    geometry_quality=geometry_quality,
                    metadata=metadata,
                    passed=False,
                    reason=f"Geometry quality {geometry_quality} < {filters.min_geometry_quality}",
                )
            )
            return parts
        parts.geometry_quality_only.append(
            CombinedFilterResult(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                metadata=metadata,
            )
        )
        return parts

    parts.processed.append(
        CombinedFilterResult(
            pdb_id=pdb_id,
            input_file=input_file,
            geometry_quality=geometry_quality,
            metadata=metadata,
            passed=False,
            reason="No UniProt accession, no resolution, no geometry quality",
        )
    )
    return parts


def _partition_structure_file_and_apply_alphafold_filter(
    input_file: Path,
    /,
    *,
    scores: dict[str, Scores],
    filters: CombinedFilterQuery,
    output_dir: Path,
    copy_method: CopyMethod = "hardlink",
) -> CombinedPartitions:
    # See combined_filter docstring for partitioning logic.
    parts = CombinedPartitions()

    try:
        structure = read_structure(input_file)
    except Exception as e:  # noqa: BLE001
        logger.warning(f"Failed to read structure {input_file}: {e}")
        parts.processed.append(
            CombinedFilterResult(
                input_file=input_file,
                passed=False,
                reason=f"Failed to read structure: {e}",
            )
        )
        return parts

    pdb_id = structure.name
    lowered_pdb_id = pdb_id.lower()
    geometry_quality = scores[lowered_pdb_id].geometry_quality if lowered_pdb_id in scores else None

    try:
        metadata = structure_metadata(structure, path=input_file)
    except Exception as e:  # noqa: BLE001
        logger.warning(f"Failed to extract metadata from {input_file}: {e}")
        parts.processed.append(
            CombinedFilterResult(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                passed=False,
                reason=f"Failed to extract metadata: {e}",
            )
        )
        return parts

    if metadata.is_alphafold:
        parts.processed.append(
            _apply_alphafold_filter(
                input_file=input_file,
                structure=structure,
                metadata=metadata,
                filters=filters,
                geometry_quality=geometry_quality,
                output_dir=output_dir,
                copy_method=copy_method,
            )
        )
        return parts

    return _apply_non_alphafold_filters(
        input_file=input_file,
        pdb_id=pdb_id,
        metadata=metadata,
        geometry_quality=geometry_quality,
        filters=filters,
        parts=parts,
    )


def _partition_structure_files_and_apply_alphafold_filter(
    input_files: list[Path],
    scores: dict[str, Scores],
    filters: CombinedFilterQuery,
    /,
    *,
    scheduler_address: SchedulerAddress = None,
    output_dir: Path,
    copy_method: CopyMethod = "hardlink",
) -> CombinedPartitions:
    file_partitions = map_with_progress(
        scheduler_address,
        _partition_structure_file_and_apply_alphafold_filter,
        input_files,
        {"tqdm_desc": "Dividing structure files into partitions", "tqdm_unit": "file"},
        scores=scores,
        filters=filters,
        output_dir=output_dir,
        copy_method=copy_method,
    )
    return CombinedPartitions().extend(file_partitions)


def _top_flat_results(
    flat_results: list[CombinedFilterResult], top_non_uniprot: int | None, output_dir: Path, copy_method: CopyMethod
) -> list[CombinedFilterResult]:
    sorted_flat_results = sort_structures(
        flat_results,
    )
    flat_top = top_non_uniprot if top_non_uniprot is not None else len(sorted_flat_results)
    processed_flat: list[CombinedFilterResult] = []
    for i, result in enumerate(sorted_flat_results):
        if i < flat_top and result.input_file:
            output_file = output_dir / result.input_file.name
            fresult = replace(result, passed=True, output_file=output_file)
            copyfile(result.input_file, output_file, copy_method)
        else:
            fresult = replace(
                result,
                passed=False,
                reason=f"Sorted index {i + 1} > {top_non_uniprot}",
            )
        processed_flat.append(fresult)
    return processed_flat


def _top_clustered_results(
    clusterable_structures: Iterable[ResolutionOrGeometryQualityClusterableStructure],
    /,
    *,
    top_uniprot_cluster: int | None,
    minimal_geometry_quality: float,
    output_dir: Path,
    copy_method: CopyMethod,
) -> list[CombinedFilterResult]:
    accession_groups: dict[str, list[ResolutionOrGeometryQualityClusterableStructure]] = {}
    for qs in clusterable_structures:
        uniprot_accession = qs.uniprot_accession
        accession_groups.setdefault(uniprot_accession, []).append(qs)

    all_clusters: list[list[ResolutionOrGeometryQualityClusterableStructure]] = []
    for structures in accession_groups.values():
        clusters = cluster_structures(structures)
        all_clusters.extend(clusters)

    if not all_clusters:
        return []

    if top_uniprot_cluster is None:
        top_uniprot_cluster = max(len(cluster) for cluster in all_clusters)

    results: list[CombinedFilterResult] = []
    for cluster in all_clusters:
        for i, qs in enumerate(cluster):
            in_top = i < top_uniprot_cluster
            ok_quality = qs.geometry_quality >= minimal_geometry_quality if qs.geometry_quality is not None else True
            if in_top and ok_quality:
                output_file = output_dir / qs.input_file.name
                fresult = CombinedFilterResult(
                    pdb_id=qs.id,
                    input_file=qs.input_file,
                    geometry_quality=qs.geometry_quality,
                    metadata=qs.metadata,
                    passed=True,
                    output_file=output_file,
                )
                copyfile(qs.input_file, output_file, copy_method)
                results.append(fresult)
            else:
                reason = []
                if not in_top:
                    reason.append(f"Sorted index {i + 1} > {top_uniprot_cluster}")
                if not ok_quality:
                    reason.append(f"Geometry quality {qs.geometry_quality} below minimum {minimal_geometry_quality}")
                fresult = CombinedFilterResult(
                    pdb_id=qs.id,
                    input_file=qs.input_file,
                    geometry_quality=qs.geometry_quality,
                    metadata=qs.metadata,
                    passed=False,
                    reason="; ".join(reason),
                )
                results.append(fresult)

    return results


def _cluster_and_take_tops(
    partitions: CombinedPartitions,
    filters: CombinedFilterQuery,
    /,
    *,
    output_dir: Path,
    copy_method: CopyMethod = "hardlink",
) -> list[CombinedFilterResult]:
    processed_uniprot_with_resolution = _top_clustered_results(
        partitions.uniprot_with_resolution,
        top_uniprot_cluster=filters.top_uniprot_cluster,
        minimal_geometry_quality=filters.min_geometry_quality,
        output_dir=output_dir,
        copy_method=copy_method,
    )
    partitions.processed.extend(processed_uniprot_with_resolution)

    processed_uniprot_with_geometry_quality = _top_clustered_results(
        partitions.uniprot_with_geometry_quality,
        top_uniprot_cluster=filters.top_uniprot_cluster,
        minimal_geometry_quality=filters.min_geometry_quality,
        output_dir=output_dir,
        copy_method=copy_method,
    )
    partitions.processed.extend(processed_uniprot_with_geometry_quality)

    processed_resolution_only = _top_flat_results(
        partitions.resolution_only, filters.top_non_uniprot, output_dir=output_dir, copy_method=copy_method
    )
    partitions.processed.extend(processed_resolution_only)

    processed_geometry_quality_only = _top_flat_results(
        partitions.geometry_quality_only, filters.top_non_uniprot, output_dir=output_dir, copy_method=copy_method
    )
    partitions.processed.extend(processed_geometry_quality_only)

    return partitions.processed


def combined_filter(
    input_files: list[Path],
    scores: dict[str, Scores],
    filters: CombinedFilterQuery,
    output_dir: Path,
    copy_method: CopyMethod = "hardlink",
    scheduler_address: SchedulerAddress = None,
) -> list[CombinedFilterResult]:
    """Combined filter for PDB/mmCIF files.

    Rules:

    * All non-AlphaFold structures are filtered by number of residues in chain A.
        See `protein-quest filter residue --help` for details.
    * All non-AlphaFold structures with Uniprot accession are filtered by sequence identity.
    * AlphaFold structures are filtered by confidence (plDDT) and afterwards by number of residues in chain A.
        See `protein-quest filter confidence --help` for details.
    * Structures with uniprot accession and resolution are filtered
        by grouping/clustering and  sorting cluster members by resolution.
        See `protein-quest filter resolution --help` for details.
    * Structures with Uniprot accesion and without resolution are filtered by grouping/clustering
        and sorting cluster members by PDBe quality scores.
    * Structures without Uniprot accession and with resolution are filtered/sorted by resolution.
    * Structures without Uniprot accession and without resolution are filtered by PDBe quality scores.
    * Structures without Uniprot accession, without resolution and without PDBe quality scores are discarded.

    ```mermaid
    flowchart TD
        A[Input PDB/mmCIF files] --> B{AlphaFold structure?}

        B -->|Yes| C[Filter by confidence plDDT]
        C --> D[Filter by residues in chain A]

        B -->|No| F[Filter by residues in chain A]
        F --> U[Filter by sequence identity]
        U --> G{UniProt accession?}

        G -->|Yes| H{Resolution available?}
        H -->|Yes| I[Group by UniProt accession and cluster by residue ranges]
        I --> J[Sort cluster members by resolution]
        J --> K[Keep up to top_uniprot_cluster per cluster]

        H -->|No| L[Group by UniProt accession and cluster by residue ranges]
        L --> M[Sort cluster members by PDBe quality]
        M --> N[Keep up to top_uniprot_cluster per cluster]

        G -->|No| O{Resolution available?}
        O -->|No| P[Sort by PDBe quality]
        P --> R[Select up to top_non_uniprot entries]
        O -->|Yes| S[Sort by resolution]
        S --> T[Select up to top_non_uniprot entries]

        D --> Q[Write output structure files]
        K --> Q
        N --> Q
        R --> Q
        T --> Q
    ```

    Arguments:
        input_files: List of input PDB/mmCIF files.
        scores: Dictionary of PDB IDs to PDBe scores.
        filters: Combined filter query.
        output_dir: Directory to write output structure files.
        copy_method: Method to copy files to output directory.
        scheduler_address: Address of the Dask scheduler for parallel processing.

    Returns:
        List of results.
    """
    # Have to do this in two stages because we need gather all structures before we can cluster and/or sort.
    partitions = _partition_structure_files_and_apply_alphafold_filter(
        input_files,
        scores,
        filters,
        output_dir=output_dir,
        copy_method=copy_method,
        scheduler_address=scheduler_address,
    )

    return _cluster_and_take_tops(
        partitions,
        filters,
        output_dir=output_dir,
        copy_method=copy_method,
    )


def combined_filter_summary(
    results: list[CombinedFilterResult],
) -> list[str]:
    nr_total = len(results)
    nr_total_passed = sum(1 for r in results if r.passed)
    nr_total_discarded = nr_total - nr_total_passed
    alphafold_results = [r for r in results if r.metadata and r.metadata.is_alphafold]
    nr_alphafold = len(alphafold_results)
    nr_alphafold_passed = sum(1 for r in alphafold_results if r.passed)
    nr_alphafold_discarded = nr_alphafold - nr_alphafold_passed
    uniprot_results = {r for r in results if r.metadata and r.metadata.uniprot_accession is not None}
    nr_uniprot = len(uniprot_results)
    nr_uniprot_passed = sum(1 for r in uniprot_results if r.passed)
    nr_uniprot_discarded = nr_uniprot - nr_uniprot_passed
    resolution_results = {
        r for r in results if r.metadata and r.metadata.resolution is not None and r.metadata.resolution != 0.0
    }
    nr_resolution = len(resolution_results)
    nr_resolution_passed = sum(1 for r in resolution_results if r.passed)
    nr_resolution_discarded = nr_resolution - nr_resolution_passed
    geometry_quality_results = {r for r in results if r.geometry_quality is not None}
    nr_geometry_quality = len(geometry_quality_results)
    nr_geometry_quality_passed = sum(1 for r in geometry_quality_results if r.passed)
    nr_geometry_quality_discarded = nr_geometry_quality - nr_geometry_quality_passed
    return [
        f"Total structures: {nr_total}, passed: {nr_total_passed}, discarded: {nr_total_discarded}",
        "",  # rich will render newline
        f"AlphaFold structures: {nr_alphafold}, passed: {nr_alphafold_passed}, discarded: {nr_alphafold_discarded}",
        "",
        (f"Structures with UniProt accession: {nr_uniprot}, "
        f"passed: {nr_uniprot_passed}, "
        f"discarded: {nr_uniprot_discarded}"),
        "",
        (f"Structures with resolution: {nr_resolution}, "
        f"passed: {nr_resolution_passed}, "
        f"discarded: {nr_resolution_discarded}"),
        "",
        (f"Structures with geometry quality: {nr_geometry_quality}, "
        f"passed: {nr_geometry_quality_passed}, "
        f"discarded: {nr_geometry_quality_discarded}"),
    ]


def combined_filter_stats(results: list[CombinedFilterResult], path: StdioPath):
    fieldnames = [
        "input_file",
        "pdb_id",
        "uniprot_accession",
        "resolution",
        "high_confidence_residues_count",
        "total_residue_count",
        "method",
        "is_alphafold",
        "uniprot_start",
        "uniprot_end",
        "sequence_identity",
        "chain_length",
        "geometry_quality",
        "passed",
        "output_file",
        "reason",
    ]
    with path.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row: dict[str, str | Path | int | float | None] = {
                "input_file": result.input_file,
                "pdb_id": result.pdb_id,
                "geometry_quality": result.geometry_quality,
                "high_confidence_residues_count": result.high_confidence_residues_count,
                "resolution": result.resolution_value,
                "sequence_identity": result.sequence_identity,
                "chain_length": result.chain_length,
                "total_residue_count": result.chain_length,
                "passed": result.passed,
                "output_file": result.output_file,
                "reason": result.reason,
            }
            if result.metadata:
                row.update(
                    {
                        "uniprot_accession": result.metadata.uniprot_accession,
                        "method": result.metadata.method,
                        "is_alphafold": result.metadata.is_alphafold,
                        "uniprot_start": result.metadata.uniprot_start,
                        "uniprot_end": result.metadata.uniprot_end,
                    }
                )
            writer.writerow(row)
