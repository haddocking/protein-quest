from asyncio.log import logger
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated

from cyclopts import Parameter
from cyclopts.types import PositiveInt
from cyclopts.validators import Number

from protein_quest.filters.quality import QualityStructure, UnclusteredStructure
from protein_quest.parallel import SchedulerAddress, map_with_progress
from protein_quest.pdbe.ws import Scores
from protein_quest.structure.formats import read_structure
from protein_quest.structure.metadata import structure_metadata
from protein_quest.structure.types import StructureMethod


@Parameter(name="*")
@dataclass(frozen=True, slots=True)
class CombinedFilterQuery:
    """Query object to apply combined filtering.

    Parameters:
        min_confidence: Minimal confidence (plDDT) for AlphaFold structures to pass the filter.
        min_residues: Min residues in chain A.
        max_residues: Max residues in chain A.
        min_geometry_quality: Minimal geometry quality score to pass the filter.
        top_uniprot_cluster: Maximum number of files to keep for structures per cluster per Uniprot accession.
            Alphafold structures are excluded from this limit.
            The top N structures in each cluster of the resolution clusters.
            The top N structures in each cluster of the PDBe quality clusters.
        top_non_uniprot: Maximum number of files to keep for structures without Uniprot accession.

    """

    min_confidence: Annotated[float, Parameter(validator=(Number(lte=100, gte=0)))] = 70.0
    min_residues: PositiveInt = 0
    max_residues: PositiveInt = 10_000_000
    min_geometry_quality: float = 50.0
    top_uniprot_cluster: PositiveInt | None = 1_000
    top_non_uniprot: PositiveInt | None = 0


@dataclass(frozen=True, slots=True)
class FilterCombinedResult:
    input_file: Path | None = None
    pdb_id: str | None = None
    uniprot_accession: str | None = None
    resolution: float | None = None
    high_confidence_residues_count: int | None = None
    total_residue_count: int | None = None
    method: StructureMethod | None = None
    uniprot_start: int | None = None
    uniprot_end: int | None = None
    sequence_identity: float | None = None
    chain_length: int | None = None
    geometry_quality: float | None = None
    passed: bool = False
    reason: str | None = None


@dataclass(frozen=True, slots=True)
class ClusterableResolutionStructure:
    id: str
    uniprot_accession: str
    uniprot_start: int
    uniprot_end: int
    sequence_identity: float
    chain_length: int
    resolution: float
    input_file: Path

    def __hash__(self) -> int:
        return hash(self.input_file)

    @property
    def geometry_quality(self) -> float:
        return 0.0


@dataclass(frozen=True, slots=True)
class UnclusteredStructureWithResolution:
    input_file: Path
    pdb_id: str
    chain_length: int
    sequence_identity: float
    resolution_value: float


@dataclass(frozen=True, slots=True)
class CombinedPartitions:
    alpha_fold: list[Path] = field(default_factory=list)
    uniprot_with_resolution: list[ClusterableResolutionStructure] = field(default_factory=list)
    uniprot_with_geometry_quality: list[QualityStructure] = field(default_factory=list)
    resolution_only: list[UnclusteredStructureWithResolution] = field(default_factory=list)
    geometry_quality_only: list[UnclusteredStructure] = field(default_factory=list)
    discarded: list[FilterCombinedResult] = field(default_factory=list)

    def extend(self, others: Iterable["CombinedPartitions"]) -> "CombinedPartitions":
        for other in others:
            self.alpha_fold.extend(other.alpha_fold)
            self.uniprot_with_resolution.extend(other.uniprot_with_resolution)
            self.uniprot_with_geometry_quality.extend(other.uniprot_with_geometry_quality)
            self.resolution_only.extend(other.resolution_only)
            self.geometry_quality_only.extend(other.geometry_quality_only)
            self.discarded.extend(other.discarded)
        return self


def _structure_file2combinedpartition(
    input_file: Path,
    /,
    *,
    scores: dict[str, Scores],
    min_residues: int,
    max_residues: int,
) -> CombinedPartitions:
    parts = CombinedPartitions()

    try:
        structure = read_structure(input_file)
    except Exception as e:  # noqa: BLE001
        logger.warning(f"Failed to read structure {input_file}: {e}")
        parts.discarded.append(
            FilterCombinedResult(
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
        parts.discarded.append(
            FilterCombinedResult(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                passed=False,
                reason=f"Failed to extract metadata: {e}",
            )
        )
        return parts

    if metadata.is_alphafold:
        # call filter_files_on_confidence later with it
        parts.alpha_fold.append(input_file)
        return parts

    if not (min_residues <= metadata.chain_length <= max_residues):
        parts.discarded.append(
            FilterCombinedResult(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                resolution=metadata.resolution,
                uniprot_accession=metadata.uniprot_accession,
                chain_length=metadata.chain_length,
                sequence_identity=metadata.sequence_identity,
                uniprot_start=metadata.uniprot_start,
                uniprot_end=metadata.uniprot_end,
                method=metadata.method,
                passed=False,
                reason=f"Chain length {metadata.chain_length} not in range [{min_residues}, {max_residues}]",
            )
        )
        return parts

    # has uniprot with geometry_quality
    if metadata.uniprot_accession and geometry_quality is not None:
        parts.uniprot_with_geometry_quality.append(
            QualityStructure(
                id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                uniprot_accession=metadata.uniprot_accession,
                chain_length=metadata.chain_length,
                sequence_identity=metadata.sequence_identity,
                uniprot_end=metadata.uniprot_end,
                uniprot_start=metadata.uniprot_start,
            )
        )
        return parts

    # has uniprot with resolution
    if metadata.uniprot_accession and metadata.resolution != 0.0:
        parts.uniprot_with_resolution.append(
            ClusterableResolutionStructure(
                id=pdb_id,
                input_file=input_file,
                resolution=metadata.resolution,
                uniprot_accession=metadata.uniprot_accession,
                chain_length=metadata.chain_length,
                sequence_identity=metadata.sequence_identity,
                uniprot_end=metadata.uniprot_end,
                uniprot_start=metadata.uniprot_start,
            )
        )
        return parts

    # has resolution only
    if metadata.resolution != 0.0:
        parts.resolution_only.append(
            UnclusteredStructureWithResolution(
                pdb_id=pdb_id,
                input_file=input_file,
                resolution_value=metadata.resolution,
                chain_length=metadata.chain_length,
                sequence_identity=metadata.sequence_identity,
            )
        )
        return parts

    # has geometry_quality only
    if geometry_quality is not None:
        parts.geometry_quality_only.append(
            UnclusteredStructure(
                pdb_id=pdb_id,
                input_file=input_file,
                geometry_quality=geometry_quality,
                chain_length=metadata.chain_length,
                sequence_identity=metadata.sequence_identity,
            )
        )
        return parts

    return parts


def partition_structures_for_combined_filter(
    input_files: list[Path],
    scores: dict[str, Scores],
    /,
    *,
    min_residues: int,
    max_residues: int,
    scheduler_address: SchedulerAddress = None,
) -> CombinedPartitions:
    file_partitions = map_with_progress(
        scheduler_address,
        _structure_file2combinedpartition,
        input_files,
        {"tqdm_desc": "Building clusters from PDBe quality", "tqdm_unit": "file"},
        scores=scores,
        min_residues=min_residues,
        max_residues=max_residues,
    )
    return CombinedPartitions().extend(file_partitions)


def combined_filter(
    input_files: list[Path],
    scores: dict[str, Scores],
    filters: CombinedFilterQuery,
    scheduler_address: SchedulerAddress = None,
) -> list[FilterCombinedResult]:
    """

    Rules:
    * All non-AlphaFold structures are filtered by number of residues in chain A.
        See `protein-quest filter residue --help` for details.
    * AlphaFold structures are filtered by confidence (plDDT) and afterwards by number of residues in chain A.
        See `protein-quest filter confidence --help` for details.
    * Structures with uniprot accession and resolution are filtered
        by grouping/clustering and  sorting cluster members by resolution.
        See `protein-quest filter resolution --help` for details.
    * Structures with Uniprot accesion and without resolution are filtered by grouping/clustering
        and sorting cluster members by PDBe quality scores.
    * Structures without Uniprot accession and without resolution are filtered by PDBe quality scores.
    * Structures without Uniprot accession and with resolution are filtered/sorted by resolution.
    * Structures without Uniprot accession, without resolution and without PDBe quality scores are discarded.

    ```mermaid
    flowchart TD
        A[Input PDB/mmCIF files] --> B{AlphaFold structure?}

        B -->|Yes| C[Filter by confidence plDDT]
        C --> D[Filter by residues in chain A]

        B -->|No| F[Filter by residues in chain A]
        F --> G{UniProt accession?}

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

    """
    partitions = partition_structures_for_combined_filter(
        input_files,
        scores,
        min_residues=filters.min_residues,
        max_residues=filters.max_residues,
        scheduler_address=scheduler_address,
    )

    # TODO apply filters, except residue filter on non-alphafold

    results: list[FilterCombinedResult] = []

    results.extend(partitions.discarded)

    return results
