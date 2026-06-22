"""Conversion orchestration for structure files."""

import logging
from collections.abc import Generator, Iterable
from pathlib import Path

from protein_quest.structure.files import split_name_and_extension
from protein_quest.structure.formats import (
    bcif2cif,
    gunzip_file,
    read_structure,
    write_structure,
)
from protein_quest.structure.types import (
    CifOutputFormat,
    Pdb2UniprotMapping,
    cif_output_formats,
    valid_structure_file_extensions,
)
from protein_quest.structure.uniprot import add_uniprot_accessions2structure
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


def convert_to_cif_files(
    input_files: Iterable[Path],
    output_dir: Path,
    copy_method: CopyMethod,
    output_format: CifOutputFormat = ".cif",
    pdb2uniprot: Pdb2UniprotMapping | None = None,
) -> Generator[tuple[Path, Path]]:
    """Convert structure files to CIF format.

    Args:
        input_files: Iterable of structure files to convert.
        output_dir: Directory to save the converted files.
        copy_method: How to copy when no changes are needed to output file.
        output_format: Output file format to write.
        pdb2uniprot: Optional dictionary mapping PDB ID to set of tuples containing chain and UniProt accession.
            If provided, will be used to inject UniProt accessions into structures that lack them.

    Yields:
        A tuple of the input file and the output file."""
    for input_file in input_files:
        output_file = convert_to_cif_file(
            input_file, output_dir, copy_method, output_format=output_format, pdb2uniprot=pdb2uniprot
        )
        yield input_file, output_file


def convert_to_cif_file(
    input_file: Path,
    output_dir: Path,
    copy_method: CopyMethod,
    output_format: CifOutputFormat = ".cif",
    pdb2uniprot: Pdb2UniprotMapping | None = None,
) -> Path:
    """Convert a single structure file to CIF format.

    Args:
        input_file: The structure file to convert.
            See [StructureFileExtensions][protein_quest.io.StructureFileExtensions]
            for supported extensions.
        output_dir: Directory to save the converted file.
        copy_method: How to copy when no changes are needed to output file.
        output_format: Output file format to write.
        pdb2uniprot: Optional dictionary mapping PDB ID to set of tuples containing chain and UniProt accession.
            If provided, will not use any shortcuts for copying files and
            will always read and write the structure to ensure UniProt accessions
            are verified and injected if necessary.

    Returns:
        Path to the converted file.

    Raises:
        ValueError: If the requested output format is not supported."""
    if output_format not in cif_output_formats:
        msg = f"Unsupported output format {output_format}. Supported output formats are: {cif_output_formats}."
        raise ValueError(msg)

    name, extension = split_name_and_extension(input_file.name)
    output_file = output_dir / f"{name}{output_format}"
    if output_file.exists():
        logger.info("Output file %s already exists for input file %s. Skipping.", output_file, input_file)
    elif pdb2uniprot or extension in {".pdb", ".pdb.gz", ".ent", ".ent.gz"}:
        structure = read_structure(input_file)
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)
        write_structure(new_structure, output_file)
    elif extension == ".cif":
        if output_format == ".cif":
            logger.info("File %s is already in .cif format, copying to %s", input_file, output_dir)
            copyfile(input_file, output_file, copy_method)
        else:
            structure = read_structure(input_file)
            write_structure(structure, output_file)
    elif extension == ".cif.gz":
        if output_format == ".cif":
            gunzip_file(input_file, output_file=output_file, keep_original=True)
        else:
            copyfile(input_file, output_file, copy_method)
    elif extension == ".bcif":
        if output_format == ".cif":
            with output_file.open("w") as f:
                f.write(bcif2cif(input_file))
        else:
            structure = read_structure(input_file)
            write_structure(structure, output_file)
    else:
        msg = (
            f"Unsupported file extension {extension} in {input_file}. "
            f"Supported extensions are {valid_structure_file_extensions}."
        )
        raise ValueError(msg)
    return output_file
