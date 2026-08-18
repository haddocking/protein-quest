"""Conversion orchestration for structure files."""

import csv
import logging
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from protein_quest.structure.chains import ChainIdSystem
from protein_quest.structure.files import split_name_and_extension
from protein_quest.structure.formats import (
    bcif2cif,
    gunzip_file,
    read_structure,
    write_structure,
)
from protein_quest.structure.types import (
    CifOutputFormat,
    cif_output_formats,
    valid_structure_file_extensions,
)
from protein_quest.structure.uniprot import add_uniprot_accessions2structure
from protein_quest.uniprot_chains import (
    Pdb2UniprotChainsMapping,
    UniprotChainMappings,
    format_uniprot_chain_mappings,
)
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


@dataclass
class ConversionStatistics:
    """Statistics for a single file conversion.

    Attributes:
        input_file: Path to the input file.
        output_file: Path to the output file.
        uniprot_chain_mappings: Set of UniprotChainMapping that were injected, empty set if none.
    """

    input_file: Path
    output_file: Path
    uniprot_chain_mappings: UniprotChainMappings


def _validate_output_format(output_format: CifOutputFormat) -> None:
    if output_format not in cif_output_formats:
        msg = f"Unsupported output format {output_format}. Supported output formats are: {cif_output_formats}."
        raise ValueError(msg)


def _new_conversion_statistics(
    input_file: Path,
    output_file: Path,
    uniprot_chain_mappings: UniprotChainMappings | None = None,
) -> ConversionStatistics:
    return ConversionStatistics(
        input_file=input_file,
        output_file=output_file,
        uniprot_chain_mappings=set() if uniprot_chain_mappings is None else uniprot_chain_mappings,
    )


def _handle_existing_output_file(input_file: Path, output_file: Path) -> ConversionStatistics:
    logger.info("Output file %s already exists for input file %s. Skipping.", output_file, input_file)
    return _new_conversion_statistics(input_file, output_file)


def _handle_pdb_like_input(
    input_file: Path,
    output_file: Path,
    output_format: CifOutputFormat,
    copy_method: CopyMethod,
    pdb2uniprot: Pdb2UniprotChainsMapping | None,
    chain_system: ChainIdSystem,
    extension: str,
) -> ConversionStatistics:
    structure = read_structure(input_file)
    new_structure, _, uniprot_chain_mappings = add_uniprot_accessions2structure(
        structure, pdb2uniprot, chain_system=chain_system
    )
    if structure is new_structure and extension == output_format:
        msg = "File %s is already in %s format and does not need change, copying to %s"
        logger.info(msg, input_file, output_format, output_file.parent)
        copyfile(input_file, output_file, copy_method)
        return _new_conversion_statistics(input_file, output_file, uniprot_chain_mappings)

    write_structure(new_structure, output_file)
    return _new_conversion_statistics(input_file, output_file, uniprot_chain_mappings)


def _handle_cif_input(
    input_file: Path,
    output_file: Path,
    output_format: CifOutputFormat,
    copy_method: CopyMethod,
) -> ConversionStatistics:
    if output_format == ".cif":
        logger.info("File %s is already in .cif format, copying to %s", input_file, output_file.parent)
        copyfile(input_file, output_file, copy_method)
        return _new_conversion_statistics(input_file, output_file)

    structure = read_structure(input_file)
    write_structure(structure, output_file)
    return _new_conversion_statistics(input_file, output_file)


def _handle_cif_gz_input(
    input_file: Path,
    output_file: Path,
    output_format: CifOutputFormat,
    copy_method: CopyMethod,
) -> ConversionStatistics:
    if output_format == ".cif":
        gunzip_file(input_file, output_file=output_file, keep_original=True)
        return _new_conversion_statistics(input_file, output_file)

    copyfile(input_file, output_file, copy_method)
    return _new_conversion_statistics(input_file, output_file)


def _handle_bcif_input(
    input_file: Path,
    output_file: Path,
    output_format: CifOutputFormat,
) -> ConversionStatistics:
    if output_format == ".cif":
        with output_file.open("w") as f:
            f.write(bcif2cif(input_file))
        return _new_conversion_statistics(input_file, output_file)

    structure = read_structure(input_file)
    write_structure(structure, output_file)
    return _new_conversion_statistics(input_file, output_file)


def convert_to_cif_file(
    input_file: Path,
    output_dir: Path,
    copy_method: CopyMethod,
    output_format: CifOutputFormat = ".cif",
    pdb2uniprot: Pdb2UniprotChainsMapping | None = None,
    chain_system: ChainIdSystem = "auth",
) -> ConversionStatistics:
    """Convert a single structure file to CIF format.

    Args:
        input_file: The structure file to convert.
            See [StructureFileExtensions][protein_quest.structure.types.StructureFileExtensions]
            for supported extensions.
        output_dir: Directory to save the converted file.
        copy_method: How to copy when no changes are needed to output file.
        output_format: Output file format to write.
        pdb2uniprot: Optional dictionary mapping PDB ID to structured UniProt chain mappings.
            If provided, will not use any shortcuts for copying files and
            will always read and write the structure to ensure UniProt accessions
            are verified and injected if necessary.
        chain_system: System of chain ids in ``pdb2uniprot`` mapping.

    Returns:
        ConversionStatistics with details about the conversion.

    Raises:
        ValueError: If the requested output format is not supported."""
    logger.debug("Converting %s", input_file)
    _validate_output_format(output_format)

    name, extension = split_name_and_extension(input_file.name)
    output_file = output_dir / f"{name}{output_format}"

    if output_file.exists():
        return _handle_existing_output_file(input_file, output_file)

    if pdb2uniprot or extension in {".pdb", ".pdb.gz", ".ent", ".ent.gz"}:
        return _handle_pdb_like_input(
            input_file=input_file,
            output_file=output_file,
            output_format=output_format,
            copy_method=copy_method,
            pdb2uniprot=pdb2uniprot,
            chain_system=chain_system,
            extension=extension,
        )

    if extension == ".cif":
        return _handle_cif_input(input_file, output_file, output_format, copy_method)

    if extension == ".cif.gz":
        return _handle_cif_gz_input(input_file, output_file, output_format, copy_method)

    if extension == ".bcif":
        return _handle_bcif_input(input_file, output_file, output_format)

    msg = (
        f"Unsupported file extension {extension} in {input_file}. "
        f"Supported extensions are {valid_structure_file_extensions}."
    )
    raise ValueError(msg)


def write_conversion_stats(
    results: Iterable[ConversionStatistics],
    output_file: Path,
) -> None:
    """Write conversion statistics to a CSV file.

    Args:
        results: Iterable of ConversionStatistics to write.
        output_file: Path to write the CSV file. Use '-' for stdout.
    """

    if str(output_file) != "-":
        output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["input_file", "output_file", "injected", "uniprot_chain_mappings"])
        writer.writeheader()
        writer.writerows(
            {
                "input_file": str(r.input_file),
                "output_file": str(r.output_file),
                "injected": str(bool(r.uniprot_chain_mappings)),
                "uniprot_chain_mappings": format_uniprot_chain_mappings(r.uniprot_chain_mappings),
            }
            for r in results
        )
