"""Retrieve subcommands for protein-quest."""

import asyncio
import csv
from collections.abc import Mapping
from pathlib import Path
from typing import Annotated

from cyclopts import App, Group, Parameter, validators

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
    OutputFile,
    console,
    to_cacher,
)
from protein_quest.emdb import fetch as emdb_fetch
from protein_quest.emdb import read_emdb_ids_from_csv
from protein_quest.pdbe import fetch as pdbe_fetch
from protein_quest.pdbe.fetch import read_pdb_ids_from_csv
from protein_quest.pdbe_3dbeacons.retrieve import read_retrieve_structure_rows, retrieve_structures

rprint = console.print


def _write_pdbe_stats(result: Mapping[str, Path], output: OutputFile) -> None:
    """Write the PDBe IDs and resulting file paths to CSV."""
    if str(output) != "-":
        output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "output_file"])
        writer.writeheader()
        for pdb_id, output_file in sorted(result.items()):
            writer.writerow({"pdb_id": pdb_id, "output_file": output_file})


retrieve_app = App(name="retrieve", help="Retrieve structure files")

pdbe_archive_mode_group = Group(
    validator=validators.MutuallyExclusive(),
    default_parameter=Parameter(negative="", show_default=False),
)


@retrieve_app.command
def pdbe(
    pdbe_csv: InputFile,
    output_dir: OutputDir,
    /,
    *,
    archived: Annotated[
        bool,
        Parameter(group=pdbe_archive_mode_group),
    ] = False,
    beta_archive: Annotated[
        bool,
        Parameter(group=pdbe_archive_mode_group),
    ] = False,
    max_parallel_downloads: BatchSize = 5,
    cache: CacheParameter | None = None,
    write_stats: OutputFile | None = None,
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
        archived: Retrieve archived versions.
            By default downloads an updated version of the PDB archive mmCIF format file.
            The updated version is generated with standardisation of vocabularies,
            and addition of connectivity information for every chemical compound present in the PDB entry.
        beta_archive: Retrieve files from the wwPDB beta archive.
            Allows 4 character PDB IDs or extended PDB IDs and downloads mmCIF files like ``pdb_00001abc.cif.gz``.
        max_parallel_downloads: Maximum number of parallel downloads.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        write_stats: Write a CSV with `pdb_id` and `output_file` columns. Use `-` for stdout.
        _: Common CLI options.
    """
    pdb_ids = read_pdb_ids_from_csv(pdbe_csv)
    rprint(f"Retrieving {len(pdb_ids)} PDBe entries")

    cacher = to_cacher(cache)

    result = asyncio.run(
        pdbe_fetch.fetch(
            pdb_ids,
            output_dir,
            archived=archived,
            beta_archive=beta_archive,
            max_parallel_downloads=max_parallel_downloads,
            cacher=cacher,
        )
    )
    rprint(f"Retrieved {len(result)} PDBe entries, written to {output_dir}")
    if write_stats:
        _write_pdbe_stats(result, write_stats)
        if str(write_stats) != "-":
            rprint(f"Statistics written to {write_stats}")


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
        format_: Formats to retrieve. Defaults to [`cif`].
            Repeat parameter for multiple formats, for example `--format cif --format pdb`.
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
    total_nr_summaries = sum(1 for af in afs if af.summary is not None)
    rprint(f"Retrieved {total_nr_files} AlphaFold files and {total_nr_summaries} summaries, written to {output_dir}")


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
