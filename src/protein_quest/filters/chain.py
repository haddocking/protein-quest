"""Filter structure files by chain."""

import logging
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.structure.chains import ChainIdSystem, resolve_chain_id_to_label, write_single_chain_structure_file
from protein_quest.structure.formats import read_structure
from protein_quest.utils import CopyMethod

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
    chain_system: ChainIdSystem = "auth",
    copy_method: CopyMethod = "copy",
) -> ChainFilterStatistics:
    input_file, chain_id = file_and_chain
    logger.debug("Filtering %s on chain %s (%s system)", input_file, chain_id, chain_system)
    try:
        structure = read_structure(input_file)
        label_chain_id = resolve_chain_id_to_label(
            structure,
            chain_id,
            chain_system=chain_system,
            source_file=input_file,
        )
        if label_chain_id != chain_id:
            logger.info(
                "Resolved chain id %s -> %s for %s",
                chain_id,
                label_chain_id,
                input_file,
            )
        output_file = write_single_chain_structure_file(
            input_file,
            label_chain_id,
            output_dir,
            out_chain=out_chain,
            copy_method=copy_method,
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
    chain_system: ChainIdSystem = "auth",
    copy_method: CopyMethod = "copy",
) -> list[ChainFilterStatistics]:
    results = []
    for file_and_chain in tqdm(file2chains, unit="file"):
        result = filter_file_on_chain(
            file_and_chain,
            output_dir=output_dir,
            out_chain=out_chain,
            chain_system=chain_system,
            copy_method=copy_method,
        )
        results.append(result)
    return results


def filter_files_on_chain(
    file2chains: Collection[tuple[Path, str]],
    output_dir: Path,
    out_chain: str = "A",
    chain_system: ChainIdSystem = "auth",
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
    copy_method: CopyMethod = "copy",
) -> list[ChainFilterStatistics]:
    """Filter mmcif/PDB files by chain.

    Args:
        file2chains: Which chain to keep for each PDB file.
            First item is the PDB file path, second item is the chain ID.
        output_dir: The directory where the filtered files will be written.
        out_chain: Under what name to write the kept chain.
        chain_system: System for the given chain ids.
            If set to ``auth`` they are translated to label ids before Gemmi processing.
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
            file2chains,
            output_dir,
            out_chain=out_chain,
            chain_system=chain_system,
            copy_method=copy_method,
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
            chain_system=chain_system,
            copy_method=copy_method,
        )
