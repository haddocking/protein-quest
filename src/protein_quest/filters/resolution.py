"""Filter structure files by resolution rank."""

import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from tqdm.auto import tqdm

from protein_quest.structure import structure_metadata
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)

GroupBy = Literal["uniprot_accession"] | None
"""Type for grouping strategy in resolution-based filtering."""


@dataclass
class ResolutionFilterStatistics:
    """Statistics for filtering files based on ranked structure resolution.

    Parameters:
        input_file: The path to the input file.
        uniprot_accession: UniProt accession used for grouping.
        resolution: Resolution from the structure file.
        total_residue_count: Total residues across the whole structure.
        is_alphafold: Whether the structure was predicted by AlphaFold.
        passed: Whether the file passed the ranking filter.
        output_file: The path to the output file, if passed.
    """

    input_file: Path
    uniprot_accession: str | None
    resolution: float
    total_residue_count: int
    is_alphafold: bool
    passed: bool
    output_file: Path | None


def resolution_sort_key(stats: ResolutionFilterStatistics) -> tuple[int, float, int, str]:
    """Rank key for resolution-based filtering.

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
        metadata = structure_metadata(input_file)
        yield ResolutionFilterStatistics(
            input_file=input_file,
            uniprot_accession=metadata.uniprot_accession,
            resolution=metadata.resolution,
            total_residue_count=metadata.total_residue_count,
            is_alphafold=metadata.is_alphafold,
            passed=False,
            output_file=None,
        )


def group_resolution_statistics(
    stats: Iterable[ResolutionFilterStatistics],
    top: int,
    group_by: GroupBy = "uniprot_accession",
) -> list[ResolutionFilterStatistics]:
    """Rank stats by resolution and mark the top N as passed.

    In ``group_by='uniprot_accession'`` mode, files with no UniProt accession
    are skipped with a warning and appended last. In ``group_by=None`` mode,
    all files are ranked globally and no missing-accession warnings are emitted.

    AlphaFold structures are preferred over non-AlphaFold.
    Structures with lower resolution are preferred.
    If resolution is the same, structures with more residues are preferred.
    If resolution is missing, those structures are undesirable.

    Output order is deterministic and sorted alphabetically by filename.

    Args:
        stats: Resolution statistics to group and rank.
        top: Maximum number of structures to pass.
        group_by: Ranking strategy. ``uniprot_accession`` applies top-N per
            accession. Structures without uniprot accession are never passed.
            ``None`` applies top-N globally.

    Returns:
        All statistics with ``passed`` updated; skipped entries appended last.
    """
    if group_by is None:
        ranked = sorted(stats, key=resolution_sort_key)
        for result in ranked[:top]:
            result.passed = True
        return sorted(ranked, key=lambda item: item.input_file.name)

    grouped: dict[str, list[ResolutionFilterStatistics]] = {}
    skipped: list[ResolutionFilterStatistics] = []

    for result in stats:
        if result.uniprot_accession is None:
            logger.warning("No UniProt accession found in %s, skipping.", result.input_file)
            skipped.append(result)
            continue
        grouped.setdefault(result.uniprot_accession, []).append(result)

    for group_results in grouped.values():
        ranked = sorted(group_results, key=resolution_sort_key)
        for result in ranked[:top]:
            result.passed = True

    output: list[ResolutionFilterStatistics] = []
    for group_results in grouped.values():
        output.extend(sorted(group_results, key=lambda item: item.input_file.name))
    output.extend(skipped)
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
    group_by: GroupBy = "uniprot_accession",
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
        group_by: Ranking strategy. ``uniprot_accession`` applies top-N per
            accession. Structures without uniprot accession are never passed.
            ``None`` applies top-N globally.
        copy_method: How to copy passed files to output directory.

    Yields:
        Objects describing the filtering result for each input file.
    """
    stats = iter_resolution_statistics(input_files)
    grouped = group_resolution_statistics(stats, top, group_by=group_by)
    yield from copy_resolution_statistics(grouped, output_dir, copy_method)
