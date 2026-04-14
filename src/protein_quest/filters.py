"""Module for filtering structure files and their contents."""

import logging
from collections.abc import Collection, Generator, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.structure import nr_residues_in_chain, structure_metadata, write_single_chain_structure_file
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)

GroupBy = Literal["uniprot_accession"] | None
"""Type for grouping strategy in resolution-based filtering."""


@dataclass
class ChainFilterStatistics:
    input_file: Path
    chain_id: str
    passed: bool = False
    output_file: Path | None = None
    discard_reason: Exception | None = None


def filter_file_on_chain(
    file_and_chain: tuple[Path, str],
    output_dir: Path,
    out_chain: str = "A",
    copy_method: CopyMethod = "copy",
) -> ChainFilterStatistics:
    input_file, chain_id = file_and_chain
    logger.debug("Filtering %s on chain %s", input_file, chain_id)
    try:
        output_file = write_single_chain_structure_file(
            input_file, chain_id, output_dir, out_chain=out_chain, copy_method=copy_method
        )
        return ChainFilterStatistics(
            input_file=input_file,
            chain_id=chain_id,
            output_file=output_file,
            passed=True,
        )
    except Exception as e:  # noqa: BLE001 - error is handled downstream
        return ChainFilterStatistics(input_file=input_file, chain_id=chain_id, discard_reason=e)


def _filter_files_on_chain_sequentially(
    file2chains: Collection[tuple[Path, str]],
    output_dir: Path,
    out_chain: str = "A",
    copy_method: CopyMethod = "copy",
) -> list[ChainFilterStatistics]:
    results = []
    for file_and_chain in tqdm(file2chains, unit="file"):
        result = filter_file_on_chain(
            file_and_chain,
            output_dir=output_dir,
            out_chain=out_chain,
            copy_method=copy_method,
        )
        results.append(result)
    return results


def filter_files_on_chain(
    file2chains: Collection[tuple[Path, str]],
    output_dir: Path,
    out_chain: str = "A",
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
    copy_method: CopyMethod = "copy",
) -> list[ChainFilterStatistics]:
    """Filter mmcif/PDB files by chain.

    Args:
        file2chains: Which chain to keep for each PDB file.
            First item is the PDB file path, second item is the chain ID.
        output_dir: The directory where the filtered files will be written.
        out_chain: Under what name to write the kept chain.
        scheduler_address: The address of the Dask scheduler.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        copy_method: How to copy when a direct copy is possible.

    Returns:
        Result of the filtering process.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    if scheduler_address == "sequential":
        return _filter_files_on_chain_sequentially(
            file2chains, output_dir, out_chain=out_chain, copy_method=copy_method
        )

    # TODO make logger.debug in filter_file_on_chain show to user when --log
    # GPT-5 generated a fairly difficult setup with a WorkerPlugin, need to find a simpler approach
    with (
        configure_dask_scheduler(
            scheduler_address,
            name="filter-chain",
        ) as cluster,
        Client(cluster) as client,
    ):
        client.forward_logging()
        return dask_map_with_progress(
            client,
            filter_file_on_chain,
            file2chains,
            output_dir=output_dir,
            out_chain=out_chain,
            copy_method=copy_method,
        )


@dataclass
class ResidueFilterStatistics:
    """Statistics for filtering files based on residue count in a specific chain.

    Parameters:
        input_file: The path to the input file.
        residue_count: The number of residues.
        passed: Whether the file passed the filtering criteria.
        output_file: The path to the output file, if passed.
    """

    input_file: Path
    residue_count: int
    passed: bool
    output_file: Path | None


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


def _resolution_rank_key(stats: ResolutionFilterStatistics) -> tuple[int, float, int, str]:
    if stats.resolution != 0.0:
        return (0, stats.resolution, -stats.total_residue_count, stats.input_file.name)
    if stats.is_alphafold:
        return (1, 0.0, -stats.total_residue_count, stats.input_file.name)
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

    If 2 stuctures have the same low resolution the structure with most residues is preferred.

    Output order is deterministic and sorted alphabetically by filename.

    Args:
        stats: Resolution statistics to group and rank.
        top: Maximum number of structures to pass.
        group_by: Ranking strategy. ``uniprot_accession`` applies top-N per
            accession. ``None`` applies top-N globally.

    Returns:
        All statistics with ``passed`` updated; skipped entries appended last.
    """
    if group_by is None:
        ranked = sorted(stats, key=_resolution_rank_key)
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
        ranked = sorted(group_results, key=_resolution_rank_key)
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

    If 2 stuctures have the same low resolution the structure with most residues is preferred.


    Args:
        input_files: Structure files to rank and filter.
        output_dir: Directory where passed files will be written.
        top: Maximum number of files to keep.
        group_by: Ranking strategy. ``uniprot_accession`` applies top-N per
            accession. ``None`` applies top-N globally.
        copy_method: How to copy passed files to output directory.

    Yields:
        Objects describing the filtering result for each input file.
    """
    stats = list(iter_resolution_statistics(input_files))
    grouped = group_resolution_statistics(stats, top, group_by=group_by)
    yield from copy_resolution_statistics(grouped, output_dir, copy_method)


def filter_files_on_residues(
    input_files: list[Path],
    output_dir: Path,
    min_residues: int,
    max_residues: int,
    chain: str = "A",
    copy_method: CopyMethod = "copy",
) -> Generator[ResidueFilterStatistics]:
    """Filter PDB/mmCIF files by number of residues in given chain.

    Args:
        input_files: The list of input PDB/mmCIF files.
        output_dir: The directory where the filtered files will be written.
        min_residues: The minimum number of residues in chain.
        max_residues: The maximum number of residues in chain.
        chain: The chain to count residues of.
        copy_method: How to copy passed files to output directory:

    Yields:
        Objects containing information about the filtering process for each input file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for input_file in tqdm(input_files, unit="file"):
        residue_count = nr_residues_in_chain(input_file, chain=chain)
        passed = min_residues <= residue_count <= max_residues
        if passed:
            output_file = output_dir / input_file.name
            copyfile(input_file, output_file, copy_method)
            yield ResidueFilterStatistics(input_file, residue_count, True, output_file)
        else:
            yield ResidueFilterStatistics(input_file, residue_count, False, None)
