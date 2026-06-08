"""Filter structure files by resolution rank."""

import csv
import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from cyclopts.types import StdioPath
from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.clustering import (
    filter_structures_on_clustered_resolution,
    sort_structures,
    structure_sort_key,
    top_members_of_clusters,
)
from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.structure import StructureMetadata
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
        metadata = StructureMetadata.from_path(input_file)
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
        return list(iter_resolution_statistics(input_files))
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


def iter_resolution_statistics(
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


def _split_discarded_and_valid(
    stats: Iterable[ResolutionFilterStatistics],
) -> tuple[list[ResolutionFilterStatistics], list[ResolutionFilterStatistics]]:
    stats_list = list(stats)
    discarded = [s for s in stats_list if s.discard_reason is not None]
    valid = [s for s in stats_list if s.discard_reason is None]
    return discarded, valid


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
    """Indicates that a structure has no UniProt accession for grouping."""

    def __init__(self) -> None:
        super().__init__("No UniProt accession for grouping")

    def __hash__(self) -> int:
        return hash(())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, NoUniProtAccessionError):
            return NotImplemented
        return True


def _sort_by_resolution_and_top_per_uniprot(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed.

    Files with no UniProt accession are skipped with a warning
    and appended last.

    Items with discard_reason set are excluded from ranking but included in output.

    AlphaFold structures are preferred over non-AlphaFold.
    Structures with lower resolution are preferred.
    If resolution is the same, structures with more residues are preferred.
    If resolution is missing, those structures are undesirable.

    Args:
        stats: Resolution statistics to group and rank.
        top: Maximum number of structures to pass.

    Returns:
        All statistics with ``passed`` updated.
        The entries are sorted by uniprot accession and then sorted by
        [sort_structures][protein_quest.clustering.sort_structures].
        Entries with no uniprot accession are appended last
        as discarded before entries already given as discarded.
    """
    discarded, valid = _split_discarded_and_valid(stats)

    grouped: dict[str, list[ResolutionFilterStatistics]] = {}
    skipped: list[ResolutionFilterStatistics] = []

    for result in valid:
        if result.uniprot_accession is None:
            result.passed = False
            result.discard_reason = NoUniProtAccessionError()
            skipped.append(result)
            continue
        grouped.setdefault(result.uniprot_accession, []).append(result)

    output: list[ResolutionFilterStatistics] = []
    for group_results in grouped.values():
        ranked = sorted(group_results, key=structure_sort_key)
        for i, result in enumerate(ranked, 1):
            if i <= top:
                result.passed = True
            else:
                result.passed = False
                result.discard_reason = OutsideTopError(top=top, rank=i)
        output.extend(ranked)

    output.extend(sorted(skipped, key=structure_sort_key))
    output.extend(discarded)
    return output


def _sort_by_resolution_and_global_top(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed

    Items with discard_reason set are excluded from ranking but are appended in returned list.

    Args:
        stats: Resolution statistics to rank.
        top: Maximum number of structures to pass.

    Returns:
        All statistics with ``passed`` updated; the entries are sorted by resolution.
          See [structure_sort_key][protein_quest.clustering.structure_sort_key].

    """
    discarded, valid = _split_discarded_and_valid(stats)

    ranked = sorted(valid, key=structure_sort_key)

    for i, result in enumerate(ranked, 1):
        if i <= top:
            result.passed = True
        else:
            result.passed = False
            result.discard_reason = OutsideTopError(top=top, rank=i)

    return ranked + discarded


def _uniprot_group_sort_key(results: list[ResolutionFilterStatistics]) -> tuple[float, int, float, int, str]:
    """Take best best member in best cluster of clusters for a uniprot accession and return its sort key."""
    best_cluster = results[0]
    return structure_sort_key(best_cluster)


def _cluster_resolution_stats_by_accession(
    valid: Iterable[ResolutionFilterStatistics],
) -> tuple[list[ResolutionFilterStatistics], list[list[ResolutionFilterStatistics]]]:
    per_accession_groups: dict[str, list[ResolutionFilterStatistics]] = {}
    accessionless: list[ResolutionFilterStatistics] = []
    for result in valid:
        if result.uniprot_accession is None:
            # It does not make sense to cluster chains from different unknown residues.
            # will put them at bottom of output later
            accessionless.append(result)
            continue
        per_accession_groups.setdefault(result.uniprot_accession, []).append(result)

    clustered_groups: list[list[ResolutionFilterStatistics]] = [
        filter_structures_on_clustered_resolution(per_accession_group, top=len(per_accession_group))
        for per_accession_group in per_accession_groups.values()
    ]

    # Reorder so clusters of best uniprot accession is first
    sorted_clustered_groups = sorted(clustered_groups, key=_uniprot_group_sort_key)
    return accessionless, sorted_clustered_groups


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
    discarded, valid = _split_discarded_and_valid(stats)

    accessionless, sorted_clustered_groups = _cluster_resolution_stats_by_accession(valid)

    return _pick_top_from_clusters(sorted_clustered_groups, top, discarded, accessionless)


def _sort_by_coverage_and_global_top(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
) -> list[ResolutionFilterStatistics]:
    """Cluster resolution stats by UniProt coverage and mark the best members as passed and flatten groups.

    The input statistics are clustered by overlapping residue ranges. Within each
    cluster, members are ordered by
    [sort_structures][protein_quest.clustering.sort_structures] and the best
    entries are selected in round-robin order.

    Items with discard_reason set are excluded from clustering but included in output.

    The clustered groups are interleaved globally and the first ``top`` results
    are marked as passed.

    Structures with no UniProt accession are not clustered and appended last in deterministic order.

    Args:
        stats: Resolution statistics to cluster and rank.
        top: Maximum number of entries to mark as passed.

    Returns:
        Clustered resolution statistics with ``passed`` updated.
    """
    discarded, valid = _split_discarded_and_valid(stats)

    accessionless, sorted_clustered_groups = _cluster_resolution_stats_by_accession(valid)

    flattened = top_members_of_clusters(sorted_clustered_groups, top=len(valid))
    if accessionless:
        flattened.extend(sort_structures(accessionless))
    for i, result in enumerate(flattened, 1):
        if i <= top:
            result.passed = True
        else:
            result.passed = False
            result.discard_reason = OutsideTopError(top=top, rank=i)
    return flattened + discarded


def _pick_top_from_clusters(
    sorted_clustered_groups: list[list[ResolutionFilterStatistics]],
    top: int,
    discarded: list[ResolutionFilterStatistics],
    accessionless: list[ResolutionFilterStatistics],
) -> list[ResolutionFilterStatistics]:
    output: list[ResolutionFilterStatistics] = []
    for group_results in sorted_clustered_groups:
        for i, result in enumerate(group_results, 1):
            if i <= top:
                result.passed = True
            else:
                result.passed = False
                result.discard_reason = OutsideTopError(top=top, rank=i)
        output.extend(group_results)

    if accessionless:
        sorted_accessionless = sort_structures(accessionless)
        for i, result in enumerate(sorted_accessionless, 1):
            if i <= top:
                result.passed = True
            else:
                result.passed = False
                result.discard_reason = OutsideTopError(top=top, rank=i)
        output.extend(sorted_accessionless)

    output.extend(discarded)
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
            stat.discard_reason = ValueError(
                f"Sequence identity {stat.sequence_identity:.3f} below minimal {min_sequence_identity:.3f}"
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


def filter_files_on_resolution(
    input_files: list[Path],
    output_dir: Path,
    top: int,
    *,
    coverage: bool = False,
    group_by: bool = True,
    min_sequence_identity: float = 1.0,
    copy_method: CopyMethod = "hardlink",
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
) -> Generator[ResolutionFilterStatistics]:
    """Filter structure files by resolution rank.

    AlphaFold structures are preferred over non-AlphaFold.
    Structures with lower resolution are preferred.
    If resolution is the same, structures with more residues are preferred.
    If resolution is missing, those structures are undesirable.

    Args:
        input_files: Structure files to rank and filter.
        output_dir: Directory where passed files will be written.
        top: Maximum number of files to keep.
        coverage: Whether to cluster by coverage.
            See [cluster_structures][protein_quest.clustering.cluster_structures].
        group_by: ``True`` applies top-N per accession. Structures without
            uniprot accession are never passed. ``False`` applies top-N globally.
        min_sequence_identity: Minimum sequence identity ratio to the Uniprot sequence for a structure to be passed.
            If not set then discards structures that are not fully identical to the Uniprot sequence.
            For example if set to 0.8 then structures that have sequence identity below 0.8 are discarded.
        copy_method: How to copy passed files to output directory.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to ``sequential`` will run tasks sequentially.

    Yields:
        Objects describing the filtering result for each input file.
    """
    stats = load_resolution_statistics(input_files, scheduler_address)
    stats = filter_on_sequence_identity(min_sequence_identity, stats)
    sorted_stats = sort_resolution_statistics(stats, top, coverage=coverage, group_by=group_by)
    yield from copy_resolution_statistics(sorted_stats, output_dir, copy_method)


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
