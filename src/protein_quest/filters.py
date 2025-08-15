import logging
from pathlib import Path
from typing import cast

from dask.distributed import Client, progress
from distributed.deploy.cluster import Cluster

from protein_quest.parallel import configure_dask_scheduler
from protein_quest.pdbe.io import write_single_chain_pdb_file

logger = logging.getLogger(__name__)


def _locate_structure_file(root: Path, pdb_id: str) -> Path:
    exts = [".cif.gz", ".cif", ".pdb.gz", ".pdb"]
    # files downloaded from https://www.ebi.ac.uk/pdbe/ website
    # have file names like pdb6t5y.ent or pdb6t5y.ent.gz for a PDB formatted file.
    # TODO support pdb6t5y.ent or pdb6t5y.ent.gz file names
    for ext in exts:
        candidate = root / f"{pdb_id.lower()}{ext}"
        if candidate.exists():
            return candidate
    msg = f"No structure file found for {pdb_id} in {root}"
    raise FileNotFoundError(msg)


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
        input_file = _locate_structure_file(input_dir, pdb_id)
        return pdb_id, chain, write_single_chain_pdb_file(input_file, chain, output_dir, out_chain=out_chain)

    with Client(scheduler_address) as client:
        logger.info(f"Follow progress on dask dashboard at: {client.dashboard_link}")

        futures = client.map(task, id2chains.items())

        progress(futures)

        results = client.gather(futures)
        return cast("list[tuple[str,str, Path | None]]", results)
