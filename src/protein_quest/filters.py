import logging
from pathlib import Path
from shutil import copyfile
from typing import cast

from dask.distributed import Client, progress
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.parallel import configure_dask_scheduler
from protein_quest.pdbe.io import (
    glob_structure_files,
    is_chain_in_residues_range,
    locate_structure_file,
    write_single_chain_pdb_file,
)

logger = logging.getLogger(__name__)


def filter_files_on_chain(
    input_dir: Path,
    id2chains: dict[str, str],
    output_dir: Path,
    scheduler_address: str | Cluster | None = None,
    out_chain: str = "A",
) -> list[tuple[str, str, Path | None]]:
    """Filter mmcif/PDB files by chain.

    Args:
        input_dir: The directory containing the input PDB files.
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
        return cast("list[tuple[str,str, Path | None]]", results)


def filter_files_on_residues(
    input_dir: Path,
    output_dir: Path,
    min_residues: int,
    max_residues: int,
) -> tuple[int, int]:
    """Filter PDB/mmCIF files by number of residues in chain A.

    Args:
        input_dir: The directory containing the input PDB files.
        output_dir: The directory where the filtered files will be written.
        min_residues: The minimum number of residues in chain A.
        max_residues: The maximum number of residues in chain A.

    Returns:
        A tuple containing the number of passed and discarded files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    passed_count = 0
    input_files = sorted(glob_structure_files(input_dir))
    for input_file in tqdm(input_files, unit="file"):
        # TODO log the nr of residues in a csv file if --write-stats is given
        if not is_chain_in_residues_range(input_file, min_residues, max_residues, chain="A"):
            continue
        copyfile(input_file, output_dir / input_file.name)
        passed_count += 1

    discarded_count = len(input_files) - passed_count
    return passed_count, discarded_count
