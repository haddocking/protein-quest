"""Helpers for retrieving structure files from 3D Beacons search results."""

import csv
import logging
from collections.abc import Iterable
from dataclasses import dataclass
from io import StringIO, TextIOWrapper
from pathlib import Path

from protein_quest.converter import converter
from protein_quest.io import read_structure, split_name_and_extension, write_structure
from protein_quest.pdbe_3dbeacons.fetch import search_structure_provider_choices
from protein_quest.pdbe_3dbeacons.model import AppUniprotSchemaModelFormat, Provider
from protein_quest.utils import Cacher, PassthroughCacher, retrieve_files

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RetrieveStructureSummary:
    requested: int
    downloaded: int
    skipped: int
    converted: int
    final: int
    cached: int = 0


@dataclass(frozen=True)
class RetrieveStructureRow:
    provider: Provider
    model_identifier: str
    model_url: str
    model_format: AppUniprotSchemaModelFormat


def read_retrieve_structure_rows(file: TextIOWrapper | StringIO) -> list[RetrieveStructureRow]:
    """Reads structures.csv file.

    Args:
        file: A file-like object containing the CSV data.

    Returns:
        A list of `RetrieveStructureRow` objects.

    Throws:
        ClassValidationError: If any row fails validation.
    """
    return [converter.structure(row, RetrieveStructureRow) for row in csv.DictReader(file)]


def _build_structure_filename(
    provider: str,
    model_identifier: str,
    model_format: AppUniprotSchemaModelFormat,
    gzipped: bool,
) -> str:
    # Whenever AppUniprotSchemaModelFormat, changes this lookup should be updated
    format2extension: dict[AppUniprotSchemaModelFormat, str] = {
        "MMCIF": ".cif",
        "PDB": ".pdb",
        "BCIF": ".bcif",
    }
    extension = format2extension.get(model_format)
    if extension is None:
        msg = f"Unsupported model format: {model_format}"
        raise ValueError(msg)

    if gzipped:
        extension += ".gz"

    return f"{provider}~{model_identifier}{extension}"


async def _prepare_structure_downloads(
    rows: Iterable[RetrieveStructureRow],
    output_dir: Path,
    raw: bool,
    cacher: Cacher,
) -> tuple[list[tuple[str, str, bool]], int]:
    urls: list[tuple[str, str, bool]] = []

    # Tested with `uvx --from=httpx[cli] httpx -h 'Accept-Encoding' gzip --download somefile -v <url>`
    # should have Content-Encoding: gzip in the response headers
    gzipped_response_capable_providers: set[Provider] = {"pdbe", "alphafold", "alphafill", "swissmodel"}
    _gzipped_response_incapable_providers: set[Provider] = {"ped", "isoformio"}

    nr_cached = 0
    for row in rows:
        if row.provider not in search_structure_provider_choices:
            msg = f"Unsupported provider: {row.provider}"
            raise ValueError(msg)
        gzip_row = not raw and row.provider in gzipped_response_capable_providers
        filename = _build_structure_filename(row.provider, row.model_identifier, row.model_format, gzip_row)
        if raw:
            raw_file_exists = (output_dir / filename).exists()
            if raw_file_exists:
                logger.debug(f"File {filename} already exists in {output_dir}, skipping download.")
                continue
            if filename in cacher:
                logger.debug(f"File {filename} already exists in cache, skipping download.")
                await cacher.copy_from_cache(output_dir / filename)
                nr_cached += 1
        else:
            # If *.cif.gz exist, skip download
            output_file = output_dir / _build_structure_filename(row.provider, row.model_identifier, "MMCIF", True)
            if output_file.exists():
                logger.debug(f"File {output_file.name} already exists in {output_dir}, skipping download.")
                continue
            if output_file.name in cacher:
                logger.debug(f"File {output_file.name} already exists in cache, skipping download.")
                await cacher.copy_from_cache(output_file)
                nr_cached += 1
                continue
        urls.append((row.model_url, filename, gzip_row))

    return urls, nr_cached


def _finalize_downloaded_structures(
    downloaded: list[Path],
    urls: list[tuple[str, str, bool]],
    output_dir: Path,
    raw: bool,
) -> tuple[int, int]:
    downloaded_by_name = {path.name: path for path in downloaded}
    converted_count = 0
    final_written = 0
    for _, expected_name, _ in urls:
        downloaded_path = downloaded_by_name.get(expected_name)
        if downloaded_path is None:
            continue

        if raw:
            logger.debug(f"Keeping raw downloaded file {downloaded_path} as requested")
            final_written += 1
            continue

        stem, native_ext = split_name_and_extension(downloaded_path.name)
        if native_ext == ".cif.gz":
            logger.debug(f"File {downloaded_path} is already in .cif.gz format, skipping conversion")
            final_written += 1
            continue

        target = output_dir / f"{stem}.cif.gz"
        logger.info(f"Converting {downloaded_path} to {target} format")
        structure_obj = read_structure(downloaded_path)
        write_structure(structure_obj, target)
        # TODO write target to cache, but need target as bytes, not as Path
        downloaded_path.unlink()
        converted_count += 1
        final_written += 1

    return converted_count, final_written


async def retrieve_structures(
    rows: Iterable[RetrieveStructureRow],
    output_dir: Path,
    raw: bool = False,
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
) -> RetrieveStructureSummary:
    if cacher is None:
        cacher = PassthroughCacher()
    urls, nr_cached = await _prepare_structure_downloads(rows=rows, output_dir=output_dir, raw=raw, cacher=cacher)

    downloaded = await retrieve_files(
        urls,
        output_dir,
        max_parallel_downloads=max_parallel_downloads,
        desc="Downloading structure files",
        cacher=cacher,
        raise_for_not_found=False,
    )

    converted_count, final_written = _finalize_downloaded_structures(
        downloaded=downloaded,
        urls=urls,
        output_dir=output_dir,
        raw=raw,
    )

    nr_requested = len(urls)
    nr_downloaded = len(downloaded)
    return RetrieveStructureSummary(
        requested=nr_requested,
        downloaded=nr_downloaded,
        skipped=nr_requested - nr_downloaded,
        converted=converted_count,
        final=final_written,
        cached=nr_cached,
    )
