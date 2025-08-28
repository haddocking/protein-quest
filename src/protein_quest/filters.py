"""Module for filtering structure files and their contents."""

import logging
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from typing import cast

from dask.distributed import Client, progress
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.parallel import configure_dask_scheduler
from protein_quest.pdbe.io import (
    locate_structure_file,
    nr_residues_in_chain,
    write_single_chain_pdb_file,
)

logger = logging.getLogger(__name__)


def filter_files_on_chain(
    input_dir: Path,
    # TODO allow to write chain A and B of same pdb
    id2chains: dict[str, str],
    output_dir: Path,
    scheduler_address: str | Cluster | None = None,
    out_chain: str = "A",
) -> list[tuple[str, str, Path | None]]:
    """Filter mmcif/PDB files by chain.

    Args:
        input_dir: The directory containing the input mmcif/PDB files.
        id2chains: Which chain to keep for each PDB ID. Key is the PDB ID, value is the chain ID.
        output_dir: The directory where the filtered files will be written.
        scheduler_address: The address of the Dask scheduler.
        out_chain: Under what name to write the kept chain.

    Returns:
        A list of tuples containing the PDB ID, chain ID, and path to the filtered file.
        Last tuple item is None if something went wrong like chain not present.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    scheduler_address = configure_dask_scheduler(
        scheduler_address,
        name="filter-chain",
    )

    def task(id2chain: tuple[str, str]) -> tuple[str, str, Path | None]:
        pdb_id, chain = id2chain
        input_file = locate_structure_file(input_dir, pdb_id)
        return pdb_id, chain, write_single_chain_pdb_file(input_file, chain, output_dir, out_chain=out_chain)

    with Client(scheduler_address) as client:
        logger.info(f"Follow progress on dask dashboard at: {client.dashboard_link}")

        futures = client.map(task, id2chains.items())

        progress(futures)

        results = client.gather(futures)
        # TODO replace tuple with dataclass
        return cast("list[tuple[str,str, Path | None]]", results)


# TODO rename to be more unique for residue filter or make generic so it can be used for chain filter as well
@dataclass
class FilterStat:
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
    input_files: list[Path], output_dir: Path, min_residues: int, max_residues: int, chain: str = "A"
) -> Generator[FilterStat]:
    """Filter PDB/mmCIF files by number of residues in given chain.

    Args:
        input_files: The list of input PDB/mmCIF files.
        output_dir: The directory where the filtered files will be written.
        min_residues: The minimum number of residues in chain.
        max_residues: The maximum number of residues in chain.
        chain: The chain to count residues of.

    Yields:
        FilterStat objects containing information about the filtering process for each input file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for input_file in tqdm(input_files, unit="file"):
        residue_count = nr_residues_in_chain(input_file, chain=chain)
        passed = min_residues <= residue_count <= max_residues
        if passed:
            output_file = output_dir / input_file.name
            copyfile(input_file, output_file)
            yield FilterStat(input_file, residue_count, True, output_file)
        else:
            yield FilterStat(input_file, residue_count, False, None)
