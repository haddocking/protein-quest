"""Module for filtering structure files and their contents."""

import logging
from collections.abc import Collection, Generator
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.structure import nr_residues_in_chain, write_single_chain_structure_file
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


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
