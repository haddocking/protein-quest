"""Filter structure files by resolution rank."""

import csv
import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from itertools import tee
from pathlib import Path
from typing import Literal

from cyclopts.types import StdioPath
from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.clustering import (
    filter_structures_on_clustered_resolution,
    structure_sort_key,
    top_members_of_clusters,
)
from protein_quest.errors import ResolutionUnsetError
from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.structure.formats import read_structure
from protein_quest.structure.metadata import structure_metadata
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


@dataclass
class ResolutionFilterStatistics:
    """Statistics for filtering files based on ranked structure resolution.

    This class is compatible with
    [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Parameters:
        input_file: The path to the input file.
        id: Identifier of the structure.
        uniprot_accession: UniProt accession used for grouping.
        resolution: Resolution from the structure file.
        total_residue_count: Total residues across the whole structure.
        is_alphafold: Whether the structure was predicted by AlphaFold.
        uniprot_start: Lowest UniProt residue position covered by the mapped chain.
        uniprot_end: Highest UniProt residue position covered by the mapped chain.
        sequence_identity: Sequence identity of the mapped chain to UniProt.
        chain_length: Number of residues in the mapped chain.
        passed: Whether the file passed the ranking filter.
        output_file: The path to the output file, if passed.
        discard_reason: If the file was discarded, the reason for discarding it.
    """

    input_file: Path
    id: str
    uniprot_accession: str | None
    resolution: float
    total_residue_count: int
    is_alphafold: bool
    uniprot_start: int
    uniprot_end: int
    sequence_identity: float
    chain_length: int
    passed: bool
    output_file: Path | None
    discard_reason: Exception | None = None

    @property
    def resolution_value(self) -> float:
        return self.resolution

    def __hash__(self) -> int:
        return hash(
            (
                self.input_file,
                self.id,
                self.uniprot_accession,
                self.resolution,
                self.total_residue_count,
                self.is_alphafold,
                self.uniprot_start,
                self.uniprot_end,
                self.sequence_identity,
                self.chain_length,
                self.passed,
                self.output_file,
                # Make Exceptions hashable by their string
                str(self.discard_reason) if self.discard_reason else None,
            )
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ResolutionFilterStatistics):
            return NotImplemented
        return (
            self.input_file == other.input_file
            and self.id == other.id
            and self.uniprot_accession == other.uniprot_accession
            and round(self.resolution, 3) == round(other.resolution, 3)
            and self.total_residue_count == other.total_residue_count
            and self.is_alphafold == other.is_alphafold
            and self.uniprot_start == other.uniprot_start
            and self.uniprot_end == other.uniprot_end
            and round(self.sequence_identity, 3) == round(other.sequence_identity, 3)
            and self.chain_length == other.chain_length
            and self.passed == other.passed
            and self.output_file == other.output_file
            # Make Exceptions comparable by their string
            and str(self.discard_reason) == str(other.discard_reason)
        )


def _load_resolution_statistics_single(input_file: Path) -> ResolutionFilterStatistics:
    """Load resolution statistics for a single structure file.

    Args:
        input_file: Structure file to read metadata from.

    Returns:
        Statistics object with metadata filled in; ``passed`` is always
        ``False`` and ``output_file`` is always ``None``. If loading fails,
        returns statistics with default values and ``discard_reason`` set.
    """
    try:
        metadata = structure_metadata(read_structure(input_file), path=input_file)
        return ResolutionFilterStatistics(
            input_file=input_file,
            id=metadata.id,
            uniprot_accession=metadata.uniprot_accession,
            resolution=metadata.resolution,
            total_residue_count=metadata.total_residue_count,
            is_alphafold=metadata.is_alphafold,
            uniprot_start=metadata.uniprot_start,
            uniprot_end=metadata.uniprot_end,
            sequence_identity=metadata.sequence_identity,
            chain_length=metadata.chain_length,
            passed=False,
            output_file=None,
        )
    except Exception as e:  # noqa: BLE001 - error is handled downstream
        logger.warning("Failed to load metadata from %s", input_file)
        return ResolutionFilterStatistics(
            input_file=input_file,
            id=input_file.stem,
            uniprot_accession=None,
            resolution=0.0,
            total_residue_count=0,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=0,
            passed=False,
            output_file=None,
            discard_reason=e,
        )


def load_resolution_statistics(
    input_files: list[Path],
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
) -> list[ResolutionFilterStatistics]:
    """Load resolution statistics for structure files, optionally in parallel.

    Args:
        input_files: Structure files to read metadata from.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to ``sequential`` will run tasks sequentially.

    Returns:
        Statistics objects with metadata filled in; ``passed`` is always
        ``False`` and ``output_file`` is always ``None``.
    """
    if scheduler_address == "sequential" or (scheduler_address is None and not input_files):
        return list(yield_resolution_statistics(input_files))
    if scheduler_address is None:
        with configure_dask_scheduler(None, name="load-resolution-statistics") as cluster, Client(cluster) as client:
            client.forward_logging()
            return dask_map_with_progress(client, _load_resolution_statistics_single, input_files)
    with (
        configure_dask_scheduler(scheduler_address, name="load-resolution-statistics") as cluster,
        Client(cluster) as client,
    ):
        client.forward_logging()
        return dask_map_with_progress(client, _load_resolution_statistics_single, input_files)


def yield_resolution_statistics(
    input_files: Iterable[Path],
) -> Generator[ResolutionFilterStatistics]:
    """Load resolution statistics for each structure file.

    Args:
        input_files: Structure files to read metadata from.

    Yields:
        Statistics objects with metadata filled in; ``passed`` is always
        ``False`` and ``output_file`` is always ``None``.
    """
    for input_file in tqdm(input_files, unit="file"):
        yield _load_resolution_statistics_single(input_file)


def _split_on_reason(
    stats: Iterable[ResolutionFilterStatistics],
) -> tuple[list[ResolutionFilterStatistics], list[ResolutionFilterStatistics]]:
    stats_list1, stats_list2 = tee(stats)
    discarded = [s for s in stats_list1 if s.discard_reason is not None]
    valid = [s for s in stats_list2 if s.discard_reason is None]
    return valid, discarded


class OutsideTopError(ValueError):
    """Indicates that a structure was ranked outside the top N."""

    def __init__(self, *, top: int, rank: int) -> None:
        super().__init__(f"Rank {rank} > top {top}")
        self.top = top
        self.rank = rank

    def __hash__(self) -> int:
        return hash((self.top, self.rank))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OutsideTopError):
            return NotImplemented
        return self.top == other.top and self.rank == other.rank


class NoUniProtAccessionError(ValueError):
    """Indicates that a structure has no UniProt accession."""

    def __init__(self, input_file: Path) -> None:
        msg = (
            f"No UniProt accession or multiple UniProt accessions in {input_file}. "
            "Use `protein-quest convert structures --uniprots pdbe.csv` "
            "to inject UniProt accessions into files."
        )
        super().__init__(msg)
        self.input_file = input_file

    def __hash__(self) -> int:
        return hash((self.input_file,))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, NoUniProtAccessionError):
            return NotImplemented
        return self.input_file == other.input_file


def _sort_by_resolution_and_top_per_uniprot(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed.

    Structures with lower resolution are preferred.
    If resolution is the same, structures with more residues are preferred.

    Args:
        stats: Resolution statistics to group and rank.
        top: Maximum number of structures to pass.

    Returns:
        All statistics with ``passed`` updated.
        The entries are sorted by uniprot accession and then sorted by
        [sort_structures][protein_quest.clustering.sort_structures].
    """
    grouped: dict[str, list[ResolutionFilterStatistics]] = {}
    for stat in stats:
        if stat.uniprot_accession is None:
            raise NoUniProtAccessionError(stat.input_file)
        grouped.setdefault(stat.uniprot_accession, []).append(stat)

    output: list[ResolutionFilterStatistics] = []
    for group_results in grouped.values():
        ranked = sorted(group_results, key=structure_sort_key)
        for i, stat in enumerate(ranked, 1):
            if i <= top:
                stat.passed = True
            else:
                stat.passed = False
                stat.discard_reason = OutsideTopError(top=top, rank=i)
        output.extend(ranked)

    return output


def _sort_by_resolution_and_global_top(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed

    Args:
        stats: Resolution statistics to rank.
        top: Maximum number of structures to pass.

    Returns:
        All statistics with ``passed`` updated; the entries are sorted by resolution.
          See [structure_sort_key][protein_quest.clustering.structure_sort_key].

    """
    ranked = sorted(stats, key=structure_sort_key)

    for i, result in enumerate(ranked, 1):
        if i <= top:
            result.passed = True
        else:
            result.passed = False
            result.discard_reason = OutsideTopError(top=top, rank=i)

    return ranked


def _uniprot_group_sort_key(results: list[ResolutionFilterStatistics]) -> tuple[float, int, float, int, str]:
    """Take best best member in best cluster of clusters for a uniprot accession and return its sort key."""
    best_cluster = results[0]
    return structure_sort_key(best_cluster)


def _cluster_resolution_stats_by_accession(
    stats: Iterable[ResolutionFilterStatistics],
) -> list[list[ResolutionFilterStatistics]]:
    per_accession_groups: dict[str, list[ResolutionFilterStatistics]] = {}
    for stat in stats:
        if stat.uniprot_accession is None:
            raise NoUniProtAccessionError(stat.input_file)
        per_accession_groups.setdefault(stat.uniprot_accession, []).append(stat)

    clustered_groups: list[list[ResolutionFilterStatistics]] = [
        filter_structures_on_clustered_resolution(per_accession_group, top=len(per_accession_group))
        for per_accession_group in per_accession_groups.values()
    ]

    # Reorder so clusters of best uniprot accession is first
    return sorted(clustered_groups, key=_uniprot_group_sort_key)


def _sort_by_coverage_and_top_per_uniprot(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Cluster resolution stats by UniProt coverage and mark the best members as passed.

    The input statistics are first grouped by UniProt accession, then clustered
    by overlapping residue ranges. Within each cluster, members are ordered by
    [sort_structures][protein_quest.clustering.sort_structures] and the best
    entries are selected in round-robin order.

    Items with discard_reason set are excluded from clustering but included in output.

    The top ``top`` members of each grouped cluster are marked as passed and all
    results are returned in deterministic order.

    Structures with no UniProt accession are not clustered
    they are appended last in deterministic order,
    and marked as passed until the total number of passed entries reaches ``top``.

    Args:
        stats: Resolution statistics to cluster and rank.
        top: Maximum number of entries to mark as passed.

    Returns:
        Clustered resolution statistics with ``passed`` updated.
    """
    sorted_clustered_groups = _cluster_resolution_stats_by_accession(stats)

    return _pick_top_from_clusters(sorted_clustered_groups, top)


def _sort_by_coverage_and_global_top(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Cluster resolution stats by UniProt coverage and mark the best members as passed and flatten groups.

    The input statistics are clustered by overlapping residue ranges. Within each
    cluster, members are ordered by
    [sort_structures][protein_quest.clustering.sort_structures] and the best
    entries are selected in round-robin order.

    The clustered groups are interleaved globally and the first ``top`` results
    are marked as passed.

    Args:
        stats: Resolution statistics to cluster and rank.
        top: Maximum number of entries to mark as passed.

    Returns:
        Clustered resolution statistics with ``passed`` updated.
    """
    listed_stats = list(stats)
    sorted_clustered_groups = _cluster_resolution_stats_by_accession(listed_stats)

    flattened = top_members_of_clusters(sorted_clustered_groups, top=len(listed_stats))
    for i, result in enumerate(flattened, 1):
        if i <= top:
            result.passed = True
        else:
            result.passed = False
            result.discard_reason = OutsideTopError(top=top, rank=i)
    return flattened


def _pick_top_from_clusters(
    sorted_clustered_groups: list[list[ResolutionFilterStatistics]],
    top: int,
) -> list[ResolutionFilterStatistics]:
    # instead taking best member for each cluster to fill top,
    # TODO Take top N of each cluster, this way we cover each domain of the protein.
    # For example given 2 clusters with 5 and 2 members resp, taking the top 2 should return 4 items.
    # while currently it returns 2 items, one from each cluster.
    output: list[ResolutionFilterStatistics] = []
    for group_results in sorted_clustered_groups:
        for i, result in enumerate(group_results, 1):
            if i <= top:
                result.passed = True
            else:
                result.passed = False
                result.discard_reason = OutsideTopError(top=top, rank=i)
        output.extend(group_results)
    return output


def copy_resolution_statistics(
    stats: Iterable[ResolutionFilterStatistics],
    output_dir: Path,
    copy_method: CopyMethod = "copy",
) -> Generator[ResolutionFilterStatistics]:
    """Copy files for passed statistics and set their ``output_file`` path.

    Args:
        stats: Statistics with ``passed`` already set.
        output_dir: Directory where passed files will be written.
        copy_method: How to copy passed files to output directory.

    Yields:
        Statistics with ``output_file`` set for passed entries.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for result in stats:
        if result.passed:
            result.output_file = output_dir / result.input_file.name
            copyfile(result.input_file, result.output_file, copy_method)
        yield result


class SequenceIdentityBelowThresholdError(ValueError):
    """Indicates that a structure has sequence identity below the specified threshold."""

    def __init__(self, input_file: Path, sequence_identity: float, threshold: float) -> None:
        super().__init__(f"Sequence identity {sequence_identity:.3f} below minimal {threshold:.3f} for {input_file}")
        self.input_file = input_file
        self.sequence_identity = sequence_identity
        self.threshold = threshold

    def __hash__(self) -> int:
        return hash((self.input_file, self.sequence_identity, self.threshold))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SequenceIdentityBelowThresholdError):
            return NotImplemented
        return (
            self.input_file == other.input_file
            and round(self.sequence_identity, 3) == round(other.sequence_identity, 3)
            and round(self.threshold, 3) == round(other.threshold, 3)
        )


def filter_on_sequence_identity(
    min_sequence_identity: float, stats: Iterable[ResolutionFilterStatistics]
) -> Generator[ResolutionFilterStatistics]:
    """Discard statistics with sequence identity below the specified threshold.

    Args:
        min_sequence_identity: Minimum sequence identity ratio to the Uniprot sequence for a structure to be passed.
            If not set then discards structures that are not fully identical to the Uniprot sequence.
            For example if set to 0.8 then structures that have sequence identity below 0.8 are discarded.
        stats: Resolution statistics to filter.

    Yields:
        Statistics with ``passed`` set to ``False`` and ``discard_reason`` set for entries

    """
    for stat in stats:
        if stat.discard_reason is None and stat.sequence_identity < min_sequence_identity:
            logger.warning(
                "Discarding %s due to sequence identity %.3f below minimal sequence identity %.3f",
                stat.input_file,
                stat.sequence_identity,
                min_sequence_identity,
            )
            stat.discard_reason = SequenceIdentityBelowThresholdError(
                stat.input_file, stat.sequence_identity, min_sequence_identity
            )
            stat.passed = False
        yield stat


def sort_resolution_statistics(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
    *,
    coverage: bool = False,
    group_by: bool = True,
) -> list[ResolutionFilterStatistics]:
    """Sort resolution statistics and mark the top N as passed based on the specified criteria.

    Args:
        stats: Resolution statistics to sort.
            Each stat should have resolution and when group_by is True then should also have an uniprot_accession.
            Also the stat.discard_reason should be None.
        top: Maximum number of entries to mark as passed.
        coverage: Whether to cluster by coverage.
            See [cluster_structures][protein_quest.clustering.cluster_structures].
        group_by: ``True`` applies top-N per accession. ``False`` applies top-N globally.
            Structures without uniprot accession are never passed.

    Returns:
        Resolution statistics with ``passed`` updated.
    """
    if coverage:
        if group_by:
            sorted_stats = _sort_by_coverage_and_top_per_uniprot(stats, top)
        else:
            sorted_stats = _sort_by_coverage_and_global_top(stats, top)
    else:
        if group_by:
            sorted_stats = _sort_by_resolution_and_top_per_uniprot(stats, top)
        else:
            sorted_stats = _sort_by_resolution_and_global_top(stats, top)
    return sorted_stats


def _mark_as_resolution_filter_non_applicable(
    stats: Iterable[ResolutionFilterStatistics],
    *,
    coverage: bool,
    group_by: bool,
    lax: bool,
) -> Generator[ResolutionFilterStatistics]:
    needs_uniprot = group_by or coverage
    for stat in stats:
        if stat.resolution == 0.0:
            stat.discard_reason = ResolutionUnsetError(stat.input_file)
        if needs_uniprot and stat.uniprot_accession is None:
            stat.discard_reason = NoUniProtAccessionError(stat.input_file)
        if stat.discard_reason is not None and lax:
            stat.passed = True
        yield stat


def filter_files_on_resolution(
    input_files: list[Path],
    output_dir: Path,
    top: int,
    *,
    coverage: bool = False,
    group_by: bool = True,
    min_sequence_identity: float = 1.0,
    lax: bool = False,
    copy_method: CopyMethod = "hardlink",
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
) -> Generator[ResolutionFilterStatistics]:
    """Filter structure files by resolution rank.

    Args:
        input_files: Structure files to rank and filter.
        output_dir: Directory where passed files will be written.
        top: Maximum number of files to keep.
        coverage: Whether to cluster by coverage.
            See [cluster_structures][protein_quest.clustering.cluster_structures].
        group_by: ``True`` applies top-N per accession. ``False`` applies top-N globally.
        min_sequence_identity: Minimum sequence identity ratio to the Uniprot sequence for a structure to be passed.
            If not set then discards structures that are not fully identical to the Uniprot sequence.
            For example if set to 0.8 then structures that have sequence identity below 0.8 are discarded.
        lax: Whether to be lax in filtering. If ``True``, files that fail to load or without resolution
            will be marked as passed instead of being discarded.
            When ``group_by`` is ``True``, files without uniprot accessions will
            also be marked as passed instead of being discarded.
        copy_method: How to copy passed files to output directory.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to ``sequential`` will run tasks sequentially.

    Yields:
        Objects describing the filtering result for each input file.
    """
    stats = load_resolution_statistics(input_files, scheduler_address)

    stats_with_reason = _mark_as_resolution_filter_non_applicable(stats, coverage=coverage, group_by=group_by, lax=lax)
    applicable_stats, non_applicable_stats = _split_on_reason(stats_with_reason)

    seqind_stats = filter_on_sequence_identity(min_sequence_identity, applicable_stats)
    seqind_applicable_stats, seqind_non_applicable_stats = _split_on_reason(seqind_stats)

    sorted_stats = sort_resolution_statistics(seqind_applicable_stats, top, coverage=coverage, group_by=group_by)

    yield from copy_resolution_statistics(sorted_stats, output_dir, copy_method)
    yield from copy_resolution_statistics(seqind_non_applicable_stats, output_dir, copy_method)
    if lax:
        yield from copy_resolution_statistics(non_applicable_stats, output_dir, copy_method)
    else:
        yield from non_applicable_stats


def write_resolution_stats(stats: Iterable[ResolutionFilterStatistics], output: StdioPath) -> None:
    """Write resolution filter statistics to a CSV file.

    Args:
        stats: Resolution filter statistics to write.
        output: Output file path or "-" for stdout.
    """
    fieldnames = [
        "input_file",
        "id",
        "uniprot_accession",
        "resolution",
        "total_residue_count",
        "is_alphafold",
        "uniprot_start",
        "uniprot_end",
        "sequence_identity",
        "chain_length",
        "passed",
        "output_file",
        "discard_reason",
        "discard_reason_type",
    ]

    with output.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for stat in stats:
            writer.writerow(
                {
                    "input_file": stat.input_file,
                    "id": stat.id,
                    "uniprot_accession": stat.uniprot_accession or "",
                    "resolution": stat.resolution,
                    "total_residue_count": stat.total_residue_count,
                    "is_alphafold": stat.is_alphafold,
                    "uniprot_start": stat.uniprot_start,
                    "uniprot_end": stat.uniprot_end,
                    "sequence_identity": f"{stat.sequence_identity:.3f}",
                    "chain_length": stat.chain_length,
                    "passed": stat.passed,
                    "output_file": stat.output_file or "",
                    "discard_reason": str(stat.discard_reason) if stat.discard_reason else "",
                    "discard_reason_type": type(stat.discard_reason).__name__ if stat.discard_reason else "",
                }
            )
