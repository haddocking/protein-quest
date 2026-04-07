"""Filter subcommands for protein-quest."""

import csv
import os
from typing import TYPE_CHECKING

from cyclopts import App
from rich.panel import Panel

from protein_quest.alphafold.confidence import ConfidenceFilterQuery, filter_files_on_confidence
from protein_quest.cli.common import (
    CacheParameter,
    Common,
    InputDir,
    InputFile,
    MaxResidues,
    MinResidues,
    OutputDir,
    OutputFile,
    console,
    write_lines,
)
from protein_quest.filters import filter_files_on_chain, filter_files_on_residues
from protein_quest.io import (
    glob_structure_files,
    locate_structure_file,
)
from protein_quest.ss import (
    SecondaryStructureFilterQuery,
    filter_files_on_secondary_structure,
)
from protein_quest.utils import copyfile

if TYPE_CHECKING:
    from pathlib import Path

rprint = console.print


filter_app = App(name="filter", help="Filter files")


@filter_app.command
def confidence(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    filters: ConfidenceFilterQuery | None = None,
    write_stats: OutputFile | None = None,
    scheduler_address: str | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter AlphaFold mmcif/PDB files by confidence (plDDT).

    Filter AlphaFold mmcif/PDB files by confidence (plDDT).
    Passed files are written with residues below threshold removed.

    Args:
        input_dir: Directory with AlphaFold mmcif/PDB files.
        output_dir: Directory to write filtered mmcif/PDB files.
        filters: Confidence filtering criteria.
        write_stats: Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    if filters is None:
        filters = ConfidenceFilterQuery()
    output_dir.mkdir(parents=True, exist_ok=True)
    input_files = sorted(glob_structure_files(input_dir))
    nr_input_files = len(input_files)
    rprint(f"Starting confidence filtering of {nr_input_files} mmcif/PDB files in {input_dir} directory.")

    if write_stats and (write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    cache = cache or CacheParameter()

    passed_count = 0
    results = filter_files_on_confidence(
        input_files, filters, output_dir, copy_method=cache.copy_method, scheduler_address=scheduler_address
    )

    stats_lines = ["input_file,residue_count,passed,output_file"]
    for r in results:
        if r.filtered_file:
            passed_count += 1
        stats_lines.append(f"{r.input_file},{r.count},{r.filtered_file is not None},{r.filtered_file}")

    if write_stats:
        write_lines(write_stats, stats_lines)
    rprint(f"Filtered {passed_count} mmcif/PDB files by confidence, written to {output_dir} directory")
    if str(write_stats) != "-":
        rprint(f"Statistics written to {write_stats}")


@filter_app.command
def chain(
    chains: InputFile,
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    scheduler_address: str | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter on chain.

    For each input PDB/mmCIF and chain combination write a PDB/mmCIF file with just the given chain
    and rename it to chain `A`. Filtering is done in parallel using a Dask cluster.

    Args:
        chains: CSV file with `pdb_id` and `chain` columns. Other columns are ignored.
        input_dir: Directory with PDB/mmCIF files.
            Expected filenames are `{pdb_id}.cif.gz`, `{pdb_id}.cif`, `{pdb_id}.pdb.gz` or `{pdb_id}.pdb`.
        output_dir: Directory to write the single-chain PDB/mmCIF files.
            Output files are in same format as input files.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache or CacheParameter()

    with chains.open("r", encoding="utf-8") as chains_file:
        rows = list(csv.DictReader(chains_file))
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

    if errors:
        msg = f"Some structure files could not be found ({len(errors)} missing), skipping them"
        rprint(Panel(os.linesep.join(map(str, errors)), title=msg, style="red"))

    if not file2chain:
        msg = "No valid structure files found."
        raise ValueError(msg)

    results = filter_files_on_chain(
        file2chain, output_dir, scheduler_address=scheduler_address, copy_method=cache.copy_method
    )

    nr_written = len([r for r in results if r.passed])

    rprint(f"Wrote {nr_written} single-chain PDB/mmCIF files to {output_dir}.")

    for result in results:
        if result.discard_reason:
            rprint(f"[red]Discarding {result.input_file} ({result.discard_reason})[/red]")


@filter_app.command
def residue(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    min_residues: MinResidues = 0,
    max_residues: MaxResidues = 10_000_000,
    write_stats: OutputFile | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter PDB/mmCIF files by number of residues in chain A.

    Filter PDB/mmCIF files by number of residues in chain A.

    Args:
        input_dir: Directory with PDB/mmCIF files (for example from 'filter chain').
        output_dir: Directory to write filtered PDB/mmCIF files. Files are copied without modification.
        min_residues: Min residues in chain A.
        max_residues: Max residues in chain A.
        write_stats: Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache or CacheParameter()

    if write_stats and str(write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    stats_lines = ["input_file,residue_count,passed,output_file"]

    nr_passed = 0
    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)
    rprint(f"Filtering {nr_total} files in {input_dir} directory by number of residues in chain A.")

    for r in filter_files_on_residues(
        input_files,
        output_dir,
        min_residues=min_residues,
        max_residues=max_residues,
        copy_method=cache.copy_method,
    ):
        stats_lines.append(f"{r.input_file},{r.residue_count},{r.passed},{r.output_file or ''}")
        if r.passed:
            nr_passed += 1

    if write_stats:
        write_lines(write_stats, stats_lines)
    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if str(write_stats) != "-":
        rprint(f"Statistics written to {write_stats}")


@filter_app.command
def secondary_structure(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    filters: SecondaryStructureFilterQuery | None = None,
    write_stats: OutputFile | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter PDB/mmCIF files by secondary structure.

    Filter PDB/mmCIF files by secondary structure.

    Args:
        input_dir: Directory with PDB/mmCIF files.
        output_dir: Directory to write filtered PDB/mmCIF files. Files are copied without modification.
        filters: Secondary structure filtering criteria.
        write_stats: Write filter statistics to file.
            In CSV format with columns:
            `<input_file>,<nr_residues>,<nr_helix_residues>,<nr_sheet_residues>,
            <helix_ratio>,<sheet_ratio>,<passed>,<output_file>`.
            Use `-` for stdout.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    if filters is None:
        filters = SecondaryStructureFilterQuery()
    cache = cache or CacheParameter()
    output_dir.mkdir(parents=True, exist_ok=True)

    if write_stats and str(write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)

    stats_lines = [
        "input_file,nr_residues,nr_helix_residues,nr_sheet_residues,helix_ratio,sheet_ratio,passed,output_file"
    ]

    rprint(f"Filtering {nr_total} files in {input_dir} directory by secondary structure.")
    nr_passed = 0
    for input_file, result in filter_files_on_secondary_structure(input_files, query=filters):
        output_file: Path | None = None
        if result.passed:
            output_file = output_dir / input_file.name
            copyfile(input_file, output_file, cache.copy_method)
            nr_passed += 1
        stats_lines.append(
            f"{input_file},{result.stats.nr_residues},{result.stats.nr_helix_residues},"
            f"{result.stats.nr_sheet_residues},{round(result.stats.helix_ratio, 3)},"
            f"{round(result.stats.sheet_ratio, 3)},{result.passed},{output_file or ''}"
        )

    if write_stats:
        write_lines(write_stats, stats_lines)
    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if str(write_stats) != "-":
        rprint(f"Statistics written to {write_stats}")
