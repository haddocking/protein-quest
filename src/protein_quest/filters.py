"""Module for filtering structure files and their contents."""

import logging
import tarfile
from collections.abc import Collection, Generator
from contextlib import ExitStack
from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Literal

from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.io import locate_structure_file
from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.structure import (
    nr_residues_in_chain,
    write_single_chain_structure_bytes,
    write_single_chain_structure_file,
)
from protein_quest.tarutils import locate_tar_member_name, merge_tar_archives, write_bytes_to_tar
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


def make_chain_file_mapping(
    rows: list[dict[str, str]], input_dir: Path
) -> tuple[set[tuple[Path, str]], list[FileNotFoundError]]:
    """Resolve CSV ``pdb_id,chain`` rows to existing input files.

    Args:
        rows: Parsed CSV rows containing at least ``pdb_id`` and ``chain`` keys.
        input_dir: Directory containing source structure files.

    Returns:
        A tuple with:
            - set of ``(input_file, chain_id)`` entries that were resolved.
            - list of ``FileNotFoundError`` for rows that could not be resolved.
    """
    file2chain: set[tuple[Path, str]] = set()
    errors: list[FileNotFoundError] = []

    for row in rows:
        pdb_id = row["pdb_id"]
        chain = row["chain"]
        try:
            f = locate_structure_file(input_dir, pdb_id)
            file2chain.add((f, chain))
        except FileNotFoundError as e:
            errors.append(e)

    return file2chain, errors


def _read_filter_chain_input_bytes(
    pdb_id: str,
    input_path: Path,
    input_tar: tarfile.TarFile | None,
) -> tuple[str, str, bytes]:
    if input_tar is not None:
        input_name = locate_tar_member_name(input_tar, pdb_id)
        member = input_tar.extractfile(input_name)
        if member is None:
            msg = f"No structure file found for {pdb_id} in {input_tar.name}"
            raise FileNotFoundError(msg)
        return input_name, Path(input_name).name, member.read()

    input_file = locate_structure_file(input_path, pdb_id)
    return str(input_file), input_file.name, input_file.read_bytes()


def _write_filter_chain_output(
    output_name: str,
    output_bytes: bytes,
    output_path: Path,
    output_tar: tarfile.TarFile | None,
) -> bool:
    if output_tar is not None:
        write_bytes_to_tar(output_tar, output_name, output_bytes)
        return True

    output_file = output_path / output_name
    if output_file.exists():
        return False
    output_file.write_bytes(output_bytes)
    return True


@dataclass
class TarChunkResult:
    shard_tar_path: Path
    nr_written: int
    errors: list[Exception]
    discards: list[tuple[str, Exception]]


def _filter_chain_rows_chunk_to_tar(
    task: tuple[list[dict[str, str]], Path],
    input_path: Path,
) -> TarChunkResult:
    row_chunk, shard_tar_path = task
    errors: list[Exception] = []
    discards: list[tuple[str, Exception]] = []
    nr_written = 0
    seen: set[tuple[str, str]] = set()

    input_is_tar = input_path.suffix == ".tar"
    shard_tar_path.parent.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        input_tar = stack.enter_context(tarfile.open(input_path, "r")) if input_is_tar else None
        shard_tar = stack.enter_context(tarfile.open(shard_tar_path, "w"))

        for row in row_chunk:
            pdb_id = row["pdb_id"]
            chain = row["chain"]

            try:
                source_id, input_filename, input_bytes = _read_filter_chain_input_bytes(
                    pdb_id=pdb_id,
                    input_path=input_path,
                    input_tar=input_tar,
                )
            except FileNotFoundError as e:
                errors.append(e)
                continue

            key = (source_id, chain)
            if key in seen:
                continue
            seen.add(key)

            try:
                output_name, output_bytes = write_single_chain_structure_bytes(
                    input_filename=input_filename,
                    input_content=input_bytes,
                    chain2keep=chain,
                )
            except Exception as e:  # noqa: BLE001
                discards.append((input_filename, e))
                continue

            write_bytes_to_tar(shard_tar, output_name, output_bytes)
            nr_written += 1

    if nr_written == 0 and shard_tar_path.exists():
        shard_tar_path.unlink()

    return TarChunkResult(
        shard_tar_path=shard_tar_path,
        nr_written=nr_written,
        errors=errors,
        discards=discards,
    )


def _chunk_rows(rows: list[dict[str, str]], n_chunks: int) -> list[list[dict[str, str]]]:
    if not rows:
        return []
    chunk_size = max(1, ceil(len(rows) / n_chunks))
    return [rows[i : i + chunk_size] for i in range(0, len(rows), chunk_size)]


def _dedup_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[tuple[str, str]] = set()
    deduped_rows: list[dict[str, str]] = []
    for row in rows:
        key = (row["pdb_id"], row["chain"])
        if key in seen:
            continue
        seen.add(key)
        deduped_rows.append(row)
    return deduped_rows


def _filter_chain_rows_to_tar_with_dask(
    rows: list[dict[str, str]],
    input_path: Path,
    output_path: Path,
    scheduler_address: str | Cluster | None,
) -> tuple[int, list[Exception], list[tuple[str, Exception]]]:
    deduped_rows = _dedup_rows(rows)
    if not deduped_rows:
        return 0, [], []

    with configure_dask_scheduler(scheduler_address, name="filter-chain-tar") as cluster, Client(cluster) as client:
        n_workers = max(1, len(client.scheduler_info()["workers"]))
        row_chunks = _chunk_rows(deduped_rows, n_workers)
        shard_paths = [
            output_path.parent / f".{output_path.stem}.worker-{i}{output_path.suffix}" for i in range(len(row_chunks))
        ]
        chunk_tasks = list(zip(row_chunks, shard_paths, strict=False))

        results = dask_map_with_progress(
            client,
            _filter_chain_rows_chunk_to_tar,
            chunk_tasks,
            input_path=input_path,
        )

    total_written = sum(r.nr_written for r in results)
    errors: list[Exception] = [e for r in results for e in r.errors]
    discards: list[tuple[str, Exception]] = [d for r in results for d in r.discards]
    written_shards = [r.shard_tar_path for r in results if r.nr_written > 0]

    try:
        if total_written > 0:
            merge_tar_archives(written_shards, output_path)
    finally:
        for shard in shard_paths:
            if shard.exists():
                shard.unlink()

    return total_written, errors, discards


def filter_chain_rows_mixed_io(
    rows: list[dict[str, str]],
    input_path: Path,
    output_path: Path,
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
) -> tuple[int, list[Exception], list[tuple[str, Exception]]]:
    """Filter chain rows with directory/tar mixed I/O support.

    Supports all combinations of directory and ``.tar`` input/output paths.
    When writing to tar and ``scheduler_address`` is not ``"sequential"``, work
    is parallelized with Dask by writing per-worker shard tars and merging them.

    Args:
        rows: Parsed CSV rows containing at least ``pdb_id`` and ``chain`` keys.
        input_path: Input directory or ``.tar`` archive with structure files.
        output_path: Output directory or ``.tar`` archive path.
        scheduler_address: Dask scheduler address/cluster, ``None`` for local
            cluster, or ``"sequential"`` to disable Dask.

    Returns:
        A tuple with:
            - number of written filtered structures.
            - list of not-found/lookup errors.
            - list of per-input discard errors as ``(input_name, exception)``.
    """
    errors: list[Exception] = []
    discards: list[tuple[str, Exception]] = []
    nr_written = 0
    seen: set[tuple[str, str]] = set()

    input_is_tar = input_path.suffix == ".tar"
    output_is_tar = output_path.suffix == ".tar"

    if output_is_tar and scheduler_address != "sequential":
        return _filter_chain_rows_to_tar_with_dask(
            rows=rows,
            input_path=input_path,
            output_path=output_path,
            scheduler_address=scheduler_address,
        )

    with ExitStack() as stack:
        input_tar = stack.enter_context(tarfile.open(input_path, "r")) if input_is_tar else None
        output_tar = stack.enter_context(tarfile.open(output_path, "w")) if output_is_tar else None
        if not output_is_tar:
            output_path.mkdir(parents=True, exist_ok=True)

        for row in rows:
            pdb_id = row["pdb_id"]
            chain = row["chain"]

            try:
                source_id, input_filename, input_bytes = _read_filter_chain_input_bytes(
                    pdb_id=pdb_id,
                    input_path=input_path,
                    input_tar=input_tar,
                )
            except FileNotFoundError as e:
                errors.append(e)
                continue

            key = (source_id, chain)
            if key in seen:
                continue
            seen.add(key)

            try:
                output_name, output_bytes = write_single_chain_structure_bytes(
                    input_filename=input_filename,
                    input_content=input_bytes,
                    chain2keep=chain,
                )
            except Exception as e:  # noqa: BLE001
                discards.append((input_filename, e))
                continue

            if _write_filter_chain_output(output_name, output_bytes, output_path, output_tar):
                nr_written += 1

    if output_is_tar and nr_written == 0 and output_path.exists():
        output_path.unlink()

    return nr_written, errors, discards


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
