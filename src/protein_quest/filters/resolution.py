"""Filter structure files by resolution rank."""

import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.clustering import (
    filter_structures_on_clustered_resolution,
    interleave_longest,
    sort_structures,
    structure_sort_key,
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
        return hash(self.input_file)


def resolution_sort_key(stats: ResolutionFilterStatistics) -> tuple[int, float, int, str]:
    """Sort key for resolution-based filtering.

    AlphaFold structures are preferred over non-AlphaFold.
    Structures with lower resolution are preferred.
    If resolution is the same, structures with more residues are preferred.
    If resolution is missing, those structures are undesirable.

    Output is deterministic and sorted alphabetically by filename.
    """
    if stats.is_alphafold:
        return (0, 0.0, -stats.total_residue_count, stats.input_file.name)
    if stats.resolution != 0.0:
        return (1, stats.resolution, -stats.total_residue_count, stats.input_file.name)
    return (2, 0.0, -stats.total_residue_count, stats.input_file.name)


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


def group_resolution_statistics(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
    group_by: bool = True,
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed.

    In grouped mode, files with no UniProt accession are skipped with a warning
    and appended last. In ungrouped mode, all files are ranked globally and no
    missing-accession warnings are emitted.

    Items with discard_reason set are excluded from ranking but included in output.

    AlphaFold structures are preferred over non-AlphaFold.
    Structures with lower resolution are preferred.
    If resolution is the same, structures with more residues are preferred.
    If resolution is missing, those structures are undesirable.

    Args:
        stats: Resolution statistics to group and rank.
        top: Maximum number of structures to pass.
        group_by: ``True`` applies top-N per uniprot accession. Structures without
            uniprot accession are never passed. ``False`` applies top-N globally.

    Returns:
        All statistics with ``passed`` updated; skipped entries appended last.
        The entries are sorted alphabetically by filename.
    """
    stats_list = list(stats)
    discarded = [s for s in stats_list if s.discard_reason is not None]
    valid = [s for s in stats_list if s.discard_reason is None]

    if not group_by:
        ranked = sorted(valid, key=resolution_sort_key)
        for result in ranked[:top]:
            result.passed = True
        return sorted(ranked + discarded, key=lambda item: item.input_file.name)

    grouped: dict[str, list[ResolutionFilterStatistics]] = {}
    skipped: list[ResolutionFilterStatistics] = []

    for result in valid:
        if result.uniprot_accession is None:
            skipped.append(result)
            continue
        grouped.setdefault(result.uniprot_accession, []).append(result)

    for group_results in grouped.values():
        ranked = sorted(group_results, key=resolution_sort_key)
        for result in ranked[:top]:
            result.passed = True

    output: list[ResolutionFilterStatistics] = []
    for group_results in grouped.values():
        output.extend(group_results)
    output.extend(skipped)
    output.extend(discarded)
    return sorted(output, key=lambda item: item.input_file.name)


def _uniprot_group_sort_key(results: list[ResolutionFilterStatistics]) -> tuple[float, float, int, str]:
    """Take best best member in best cluster of clusters for a uniprot accession and return its sort key."""
    best_cluster = results[0]
    return structure_sort_key(best_cluster)


def coverage_group_resolution_statistics(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
    group_by: bool,
) -> list[ResolutionFilterStatistics]:
    """Cluster resolution stats by UniProt coverage and mark the best members as passed.

    The input statistics are first grouped by UniProt accession, then clustered
    by overlapping residue ranges. Within each cluster, members are ordered by
    [sort_structures][protein_quest.clustering.sort_structures] and the best
    entries are selected in round-robin order.

    Items with discard_reason set are excluded from clustering but included in output.

    When ``group_by`` is ``False``, the clustered groups are interleaved
    globally and the first ``top`` results are marked as passed. Otherwise, the
    top ``top`` members of each grouped cluster are marked as passed and all
    results are returned in deterministic order.

    Structures with no UniProt accession are not clustered.
    When group_by is ``False``, they are appended last in deterministic order.
    When group_by is ``True``, they are appended last in deterministic order,
    and marked as passed until the total number of passed entries reaches ``top``.

    Args:
        stats: Resolution statistics to cluster and rank.
        top: Maximum number of entries to mark as passed.
        group_by: ``False`` interleaves all clustered groups globally;
            ``True`` keeps groups separate per uniprot accession.

    Returns:
        Clustered resolution statistics with ``passed`` updated.
    """
    stats_list = list(stats)
    discarded = [s for s in stats_list if s.discard_reason is not None]
    valid = [s for s in stats_list if s.discard_reason is None]

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

    if not group_by:
        flattened = list(interleave_longest(*sorted_clustered_groups))
        if accessionless:
            flattened.extend(sort_structures(accessionless))
        for result in flattened[:top]:
            result.passed = True
        return flattened + discarded

    return _pick_top_from_clusters(sorted_clustered_groups, top, discarded, accessionless)


def _pick_top_from_clusters(
    sorted_clustered_groups: list[list[ResolutionFilterStatistics]],
    top: int,
    discarded: list[ResolutionFilterStatistics],
    accessionless: list[ResolutionFilterStatistics],
) -> list[ResolutionFilterStatistics]:
    output: list[ResolutionFilterStatistics] = []
    for group_results in sorted_clustered_groups:
        for result in group_results[:top]:
            result.passed = True
        output.extend(group_results)

    if accessionless:
        sorted_accessionless = sort_structures(accessionless)
        passed_count = sum(1 for r in output if r.passed)
        if passed_count < top:
            for result in sorted_accessionless[: top - passed_count]:
                result.passed = True
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


def _filter_on_sequence_identity(
    min_sequence_identity: float, stats: Iterable[ResolutionFilterStatistics]
) -> Generator[ResolutionFilterStatistics]:
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


def filter_files_on_resolution(
    input_files: list[Path],
    output_dir: Path,
    top: int,
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
    stats = _filter_on_sequence_identity(min_sequence_identity, stats)
    if coverage:
        grouped = coverage_group_resolution_statistics(stats, top, group_by=group_by)
    else:
        grouped = group_resolution_statistics(stats, top, group_by=group_by)
    yield from copy_resolution_statistics(grouped, output_dir, copy_method)
