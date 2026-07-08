"""Filter structure files by chain."""

import logging
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path

from protein_quest.parallel import SchedulerAddress, map_with_progress
from protein_quest.structure.chains import (
    ChainIdSystem,
    get_label2auth_chains,
    write_single_chain_structure_file,
)
from protein_quest.structure.errors import ChainNotFoundError
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
    force: bool = False,
) -> ChainFilterStatistics:
    """Filter a single structure file by chain.

    Args:
        file_and_chain: Tuple of (structure file path, chain ID to keep).
            The chain ID is interpreted according to the `chain_system` argument.
        output_dir: The directory where the filtered file will be written.
        out_chain: Under what name to write the kept chain.
            The chain system id is 'auth'.
        chain_system: System for the given chain ids in `file_and_chain`.
        copy_method: How to copy when a direct copy is possible.
        force: Rewrite existing single-chain outputs and overwrite existing output files.
    """
    input_file, chain_id = file_and_chain
    logger.debug("Filtering %s on chain %s (%s system)", input_file, chain_id, chain_system)
    try:
        auth_chain_id = chain_id
        if chain_system == "label":
            structure = read_structure(input_file)
            l2a = get_label2auth_chains(structure)
            try:
                auth_chain_id = l2a[chain_id]
            except KeyError:
                raise ChainNotFoundError(chain_id, input_file, set(l2a.keys())) from None
            if auth_chain_id != chain_id:
                logger.info(
                    "Resolved label chain id %s -> %s auth for %s",
                    chain_id,
                    auth_chain_id,
                    input_file,
                )

        output_file = write_single_chain_structure_file(
            input_file,
            auth_chain_id,
            output_dir,
            out_chain=out_chain,
            copy_method=copy_method,
            force=force,
        )
        return ChainFilterStatistics(
            input_file=input_file,
            chain_id=chain_id,
            output_file=output_file,
            passed=True,
        )
    except Exception as e:  # noqa: BLE001 - error is handled downstream
        return ChainFilterStatistics(input_file=input_file, chain_id=chain_id, discard_reason=e)


def filter_files_on_chain(
    file2chains: Collection[tuple[Path, str]],
    output_dir: Path,
    out_chain: str = "A",
    chain_system: ChainIdSystem = "auth",
    scheduler_address: SchedulerAddress = None,
    copy_method: CopyMethod = "copy",
    force: bool = False,
) -> list[ChainFilterStatistics]:
    """Filter mmcif/PDB files by chain.

    Args:
        file2chains: Which chain to keep for each PDB file.
            First item is the PDB file path, second item is the chain ID.
        output_dir: The directory where the filtered files will be written.
        out_chain: Under what name to write the kept chain.
        chain_system: System for the given chain ids.
        scheduler_address: The address of the Dask scheduler.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        copy_method: How to copy when a direct copy is possible.
        force: Rewrite existing single-chain outputs and overwrite existing output files.

    Returns:
        Result of the filtering process.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    return map_with_progress(
        scheduler_address,
        filter_file_on_chain,
        file2chains,
        map_with_progress_options={"tqdm_desc": "Filtering on chain", "tqdm_unit": "file"},
        output_dir=output_dir,
        out_chain=out_chain,
        chain_system=chain_system,
        copy_method=copy_method,
        force=force,
    )
