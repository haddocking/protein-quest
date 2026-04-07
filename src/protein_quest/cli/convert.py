"""Convert subcommands for protein-quest."""

from pathlib import Path
from typing import Annotated, Literal

from cyclopts import App, Parameter, validators
from rocrate_action_recorder.adapters.cyclopts import INPUT_DIR, OUTPUT_DIR, OUTPUT_FILE
from tqdm import tqdm

from protein_quest.cli.common import CacheParameter, Common, console, write_lines
from protein_quest.io import convert_to_cif_files, glob_structure_files, read_structure
from protein_quest.structure import structure2uniprot_accessions

rprint = console.print


convert_app = App(name="convert", help="Convert files between formats")


@convert_app.command
def uniprot(
    input_dir: Annotated[Path, Parameter(validator=validators.Path(exists=True, file_okay=False)), INPUT_DIR],
    output: Annotated[Path, Parameter(validator=validators.Path(dir_okay=False)), OUTPUT_FILE],
    /,
    *,
    grouped: Annotated[
        bool,
        Parameter(negative=""),
    ] = False,
    _common: Common | None = None,
) -> None:
    """Convert structure files to list of UniProt accessions.

    UniProt accessions are read from database reference of each structure.

    Args:
        input_dir: Directory with structure files. Supported extensions are .cif, .cif.gz, .pdb, .pdb.gz.
        output: Output text file with UniProt accessions (one per line). Use '-' for stdout.
        grouped: Whether to group accessions by structure file.
            If set output changes to `<structure_file1>,<acc1>\\n<structure_file1>,<acc2>` format.
        _common: Common CLI options.
    """
    input_files = sorted(glob_structure_files(input_dir))

    if grouped:
        lines = []
        for input_file in tqdm(input_files, unit="file"):
            s = read_structure(input_file)
            uniprot_accessions = structure2uniprot_accessions(s)
            lines.extend([f"{input_file},{uniprot_accession}" for uniprot_accession in sorted(uniprot_accessions)])
        write_lines(output, lines)
    else:
        uniprot_accessions: set[str] = set()
        for input_file in tqdm(input_files, unit="file"):
            s = read_structure(input_file)
            uniprot_accessions.update(structure2uniprot_accessions(s))
        write_lines(output, sorted(uniprot_accessions))


@convert_app.command
def structures(
    input_dir: Annotated[Path, Parameter(validator=validators.Path(exists=True, file_okay=False)), INPUT_DIR],
    /,
    *,
    output_dir: Annotated[Path | None, Parameter(validator=validators.Path(file_okay=False)), OUTPUT_DIR] = None,
    _format: Annotated[Literal["cif"], Parameter(name="--format")] = "cif",
    cache: CacheParameter | None = None,
    _common: Common | None = None,
) -> None:
    """Convert structure files between formats.

    Convert structure files between formats.

    Args:
        input_dir: Directory with structure files. Supported extensions are .cif, .cif.gz, .pdb, .pdb.gz.
        output_dir: Directory to write converted structure files. If not given, files are written to input_dir.
        _format: Output format to convert to.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _common: Common CLI options.
    """
    cache = cache or CacheParameter()
    output_dir = output_dir if output_dir is not None else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    input_files = sorted(glob_structure_files(input_dir))
    rprint(f"Converting {len(input_files)} files in {input_dir} directory to cif format.")

    for _ in tqdm(
        convert_to_cif_files(
            input_files,
            output_dir,
            copy_method=cache.copy_method,
        ),
        total=len(input_files),
        unit="file",
    ):
        pass

    rprint(f"Converted {len(input_files)} files into {output_dir}.")
