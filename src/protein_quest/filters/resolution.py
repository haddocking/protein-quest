"""Filter structure files by resolution rank."""

import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path

from tqdm.auto import tqdm

from protein_quest.clustering import cluster_structures, interleave_longest, sort_structures, top_members_of_clusters
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

    @property
    def resolution_value(self) -> float:
        return self.resolution


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
        metadata = StructureMetadata.from_path(input_file)
        yield ResolutionFilterStatistics(
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


def group_resolution_statistics(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
    group_by: bool = True,
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed.

    In grouped mode, files with no UniProt accession are skipped with a warning
    and appended last. In ungrouped mode, all files are ranked globally and no
    missing-accession warnings are emitted.

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
    if not group_by:
        ranked = sorted(stats, key=resolution_sort_key)
        for result in ranked[:top]:
            result.passed = True
        return sorted(ranked, key=lambda item: item.input_file.name)

    grouped: dict[str, list[ResolutionFilterStatistics]] = {}
    skipped: list[ResolutionFilterStatistics] = []

    for result in stats:
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
    return sorted(output, key=lambda item: item.input_file.name)


def _coverage_bucket_sort_key(results: list[ResolutionFilterStatistics]) -> tuple[int, float, float, int, str]:
    first = results[0]
    return (0, -first.sequence_identity, first.resolution, -first.chain_length, first.id)


def _coverage_bucket_sort_key_for_accessionless(
    results: list[ResolutionFilterStatistics],
) -> tuple[int, float, float, int, str]:
    first = results[0]
    # Mark accessionless groups as lowest priority
    return (1, -first.sequence_identity, first.resolution, -first.chain_length, first.id)


def _cluster_resolution_bucket(results: list[ResolutionFilterStatistics]) -> list[ResolutionFilterStatistics]:
    if not results:
        return []
    ordered_clusters = [sort_structures(cluster) for cluster in cluster_structures(results)]
    return top_members_of_clusters(ordered_clusters, top=len(results))


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

    When ``group_by`` is ``False``, the clustered groups are interleaved
    globally and the first ``top`` results are marked as passed. Otherwise, the
    top ``top`` members of each grouped cluster are marked as passed and all
    results are returned in deterministic order.

    Args:
        stats: Resolution statistics to cluster and rank.
        top: Maximum number of entries to mark as passed.
        group_by: ``False`` interleaves all clustered groups globally;
            ``True`` keeps groups separate.

    Returns:
        Clustered resolution statistics with ``passed`` updated.
    """
    grouped: dict[str | None, list[ResolutionFilterStatistics]] = {}
    for result in stats:
        grouped.setdefault(result.uniprot_accession, []).append(result)

    clustered_groups: list[list[ResolutionFilterStatistics]] = [
        _cluster_resolution_bucket(group_results) for group_results in grouped.values()
    ]

    accessioned_groups = [group for group in clustered_groups if group and group[0].uniprot_accession is not None]
    accessionless_groups = [group for group in clustered_groups if group and group[0].uniprot_accession is None]
    ordered_groups = sorted(accessioned_groups, key=_coverage_bucket_sort_key)
    ordered_groups.extend(sorted(accessionless_groups, key=_coverage_bucket_sort_key_for_accessionless))

    if not group_by:
        flattened = list(interleave_longest(*ordered_groups))
        for result in flattened[:top]:
            result.passed = True
        return flattened

    output: list[ResolutionFilterStatistics] = []
    for group_results in ordered_groups:
        for result in group_results[:top]:
            result.passed = True
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


def filter_files_on_resolution(
    input_files: list[Path],
    output_dir: Path,
    top: int,
    coverage: bool = False,
    group_by: bool = True,
    copy_method: CopyMethod = "copy",
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
        copy_method: How to copy passed files to output directory.

    Yields:
        Objects describing the filtering result for each input file.
    """
    stats = iter_resolution_statistics(input_files)
    if coverage:
        grouped = coverage_group_resolution_statistics(stats, top, group_by=group_by)
    else:
        grouped = group_resolution_statistics(stats, top, group_by=group_by)
    yield from copy_resolution_statistics(grouped, output_dir, copy_method)
