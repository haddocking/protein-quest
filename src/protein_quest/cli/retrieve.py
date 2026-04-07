"""Retrieve subcommands for protein-quest."""

import asyncio
from typing import Annotated

from cyclopts import App, Parameter

from protein_quest.alphafold.fetch import (
    DownloadableFormat,
    read_af_ids_from_csv,
)
from protein_quest.alphafold.fetch import (
    fetch_many as af_fetch,
)
from protein_quest.cli.common import (
    BatchSize,
    CacheParameter,
    Common,
    InputFile,
    OutputDir,
    console,
    to_cacher,
)
from protein_quest.emdb import fetch as emdb_fetch
from protein_quest.emdb import read_emdb_ids_from_csv
from protein_quest.pdbe import fetch as pdbe_fetch
from protein_quest.pdbe.fetch import read_pdb_ids_from_csv
from protein_quest.pdbe_3dbeacons.retrieve import read_retrieve_structure_rows, retrieve_structures

rprint = console.print


retrieve_app = App(name="retrieve", help="Retrieve structure files")


@retrieve_app.command
def pdbe(
    pdbe_csv: InputFile,
    output_dir: OutputDir,
    /,
    *,
    max_parallel_downloads: BatchSize = 5,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Retrieve mmCIF files from PDBe for PDB IDs in CSV.

    Retrieve mmCIF files from Protein Data Bank in Europe Knowledge Base (PDBe) website
    for unique PDB IDs listed in a CSV file.

    Args:
        pdbe_csv: CSV file with a `pdb_id` column, or with `model_provider` and
            `model_identifier` columns. When using `model_provider`, only rows
            with `model_provider == 'pdbe'` are used. Single-column CSV files
            are also accepted, and the first row is treated as an ID. Use `-` for stdin.
        output_dir: Directory to store downloaded PDBe mmCIF files.
        max_parallel_downloads: Maximum number of parallel downloads.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    pdb_ids = read_pdb_ids_from_csv(pdbe_csv)
    rprint(f"Retrieving {len(pdb_ids)} PDBe entries")

    cacher = to_cacher(cache)

    result = asyncio.run(
        pdbe_fetch.fetch(pdb_ids, output_dir, max_parallel_downloads=max_parallel_downloads, cacher=cacher)
    )
    rprint(f"Retrieved {len(result)} PDBe entries")


@retrieve_app.command
def alphafold(
    alphafold_csv: InputFile,
    output_dir: OutputDir,
    /,
    *,
    format_: Annotated[set[DownloadableFormat], Parameter(name="--format", negative="")] | None = None,
    db_version: str | None = None,
    gzip_files: Annotated[
        bool,
        Parameter(negative=""),
    ] = False,
    all_isoforms: Annotated[
        bool,
        Parameter(negative=""),
    ] = False,
    max_parallel_downloads: BatchSize = 5,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Retrieve AlphaFold files for IDs in CSV.

    Retrieve AlphaFold files from the AlphaFold Protein Structure Database.

    Args:
        alphafold_csv: CSV file with an `af_id` column, or with `model_provider` and
            `model_identifier` columns. When using `model_provider`, only rows
            with `model_provider == 'alphafold'` are used. Single-column CSV
            files are also accepted, and the first row is treated as an ID.
            Use `-` for stdin.
        output_dir: Directory to store downloaded AlphaFold files.
        format_: Formats to retrieve.
        db_version: AlphaFold database version.
        gzip_files: Gzip downloaded files.
        all_isoforms: Return all isoforms.
        max_parallel_downloads: Maximum number of parallel downloads.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    if format_ is None:
        format_ = {"cif"}

    af_ids = read_af_ids_from_csv(alphafold_csv)
    formats = format_

    rprint(f"Retrieving {len(af_ids)} AlphaFold entries with formats {formats}")

    cacher = to_cacher(cache)

    afs = af_fetch(
        af_ids,
        output_dir,
        formats=formats,
        db_version=db_version,
        max_parallel_downloads=max_parallel_downloads,
        cacher=cacher,
        gzip_files=gzip_files,
        all_isoforms=all_isoforms,
    )
    total_nr_files = sum(af.nr_of_files() for af in afs)
    rprint(f"Retrieved {total_nr_files} AlphaFold files and {len(afs)} summaries, written to {output_dir}")


@retrieve_app.command
def emdb(
    emdb_csv: InputFile,
    output_dir: OutputDir,
    /,
    *,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Retrieve EMDB volume files for EMDB IDs in CSV.

    Retrieve volume files from Electron Microscopy Data Bank (EMDB) website
    for unique EMDB IDs listed in a CSV file.

    Args:
        emdb_csv: CSV file with `emdb_id` column. Other columns are ignored.
            Single-column CSV files are also accepted, and the first row is treated as an ID.
            Use `-` for stdin.
        output_dir: Directory to store downloaded EMDB volume files.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    emdb_ids = read_emdb_ids_from_csv(emdb_csv)
    rprint(f"Retrieving {len(emdb_ids)} EMDB entries")

    cacher = to_cacher(cache)

    result = asyncio.run(emdb_fetch(emdb_ids, output_dir, cacher=cacher))
    rprint(f"Retrieved {len(result)} EMDB entries")


@retrieve_app.command
def structure(
    structures_csv: InputFile,
    output_dir: OutputDir,
    /,
    *,
    raw: Annotated[
        bool,
        Parameter(negative=""),
    ] = False,
    max_parallel_downloads: BatchSize = 5,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Retrieve structure files from search structure CSV output.

    Retrieve structure files from model URLs listed in search structure CSV output.

    Args:
        structures_csv: CSV file with `provider`, `model_identifier`, `model_url`, and `model_format` columns.
            Use `-` for stdin.
        output_dir: Directory to store retrieved structure files.
        raw: Download in native format from CSV.
        max_parallel_downloads: Maximum number of parallel downloads.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    rows = read_retrieve_structure_rows(structures_csv)

    cacher = to_cacher(cache)

    summary = asyncio.run(
        retrieve_structures(
            rows,
            output_dir,
            raw=raw,
            max_parallel_downloads=max_parallel_downloads,
            cacher=cacher,
        )
    )
    rprint(
        f"Retrieved structure files "
        f"requested={summary.requested}, downloaded={summary.downloaded}, skipped={summary.skipped}, "
        f"converted={summary.converted}, final={summary.final}, cached={summary.cached}"
    )
