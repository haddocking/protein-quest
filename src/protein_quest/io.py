"""Module for structure file input/output."""

import gzip
import logging
import shutil
import tempfile
from collections.abc import Generator, Iterable
from io import StringIO
from pathlib import Path
from typing import Literal, get_args
from urllib.request import urlopen

import gemmi
from mmcif.api.DictionaryApi import DictionaryApi
from mmcif.io.BinaryCifReader import BinaryCifReader
from mmcif.io.BinaryCifWriter import BinaryCifWriter
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter

from protein_quest.utils import CopyMethod, copyfile, user_cache_root_dir

logger = logging.getLogger(__name__)


StructureFileExtensions = Literal[".pdb", ".pdb.gz", ".ent", ".ent.gz", ".cif", ".cif.gz", ".bcif", ".bcif.gz"]
"""Type of supported structure file extensions."""
valid_structure_file_extensions: set[str] = set(get_args(StructureFileExtensions))
"""Set of valid structure file extensions."""

CifOutputFormat = Literal[".cif", ".cif.gz"]
"""Supported output formats for CIF conversion."""
cif_output_formats: set[str] = set(get_args(CifOutputFormat))
"""Set of valid CIF conversion output formats."""


def write_structure(structure: gemmi.Structure, path: Path):
    """Write a gemmi structure to a file.

    Args:
        structure: The gemmi structure to write.
        path: The file path to write the structure to.
            The format depends on the file extension.
            See [StructureFileExtensions][protein_quest.io.StructureFileExtensions]
            for supported extensions.

    Raises:
        ValueError: If the file extension is not supported.
    """
    if path.name.endswith(".pdb") or path.name.endswith(".ent"):
        body: str = structure.make_pdb_string()
        path.write_text(body)
    elif path.name.endswith(".pdb.gz") or path.name.endswith(".ent.gz"):
        body: str = structure.make_pdb_string()
        with gzip.open(path, "wt") as f:
            f.write(body)
    elif path.name.endswith(".cif"):
        # do not write chem_comp so it is viewable by molstar
        # see https://github.com/project-gemmi/gemmi/discussions/362
        doc = structure.make_mmcif_document(gemmi.MmcifOutputGroups(True, chem_comp=False))
        doc.write_file(str(path))
    elif path.name.endswith(".cif.gz"):
        path.write_bytes(structure2cifgz(structure))
    elif path.name.endswith(".bcif"):
        structure2bcif(structure, path)
    elif path.name.endswith(".bcif.gz"):
        structure2bcifgz(structure, path)
    else:
        msg = f"Unsupported file extension in {path.name}. Supported extensions are: {valid_structure_file_extensions}"
        raise ValueError(msg)


def read_structure(file: Path) -> gemmi.Structure:
    """Read a structure from a file.

    Args:
        file: Path to the input structure file.
            See [StructureFileExtensions][protein_quest.io.StructureFileExtensions]
            for supported extensions.

    Returns:
        A gemmi Structure object representing the structure in the file.
    """
    if file.name.endswith(".bcif"):
        return bcif2structure(file)
    if file.name.endswith(".bcif.gz"):
        return bcifgz2structure(file)
    return gemmi.read_structure(str(file))


def bcif2cif(bcif_file: Path) -> str:
    """Convert a binary CIF (bcif) file to a CIF string.

    Args:
        bcif_file: Path to the binary CIF file.

    Returns:
        A string containing the CIF representation of the structure.
    """
    reader = BinaryCifReader()
    container = reader.deserialize(str(bcif_file))
    capture = StringIO()
    writer = PdbxWriter(capture)
    writer.write(container)
    return capture.getvalue()


def bcifgz2structure(bcif_gz_file: Path) -> gemmi.Structure:
    """Read a binary CIF (bcif) gzipped file and return a gemmi Structure object.

    This is slower than other formats because gemmi does not support reading bcif files directly.
    So we first gunzip the file to a temporary location, convert it to a cif string using mmcif package,
    and then read the cif string using gemmi.

    Args:
        bcif_gz_file: Path to the binary CIF gzipped file.

    Returns:
        A gemmi Structure object representing the structure in the bcif.gz file.
    """
    with tempfile.NamedTemporaryFile(suffix=".bcif", delete=True) as tmp_bcif:
        tmp_path = Path(tmp_bcif.name)
        gunzip_file(bcif_gz_file, output_file=tmp_path, keep_original=True)
        return bcif2structure(tmp_path)


def bcif2structure(bcif_file: Path) -> gemmi.Structure:
    """Read a binary CIF (bcif) file and return a gemmi Structure object.

    This is slower than other formats because gemmi does not support reading bcif files directly.
    So we convert it to a cif string first using mmcif package and then read the cif string using gemmi.

    Args:
        bcif_file: Path to the binary CIF file.

    Returns:
        A gemmi Structure object representing the structure in the bcif file.
    """
    cif_content = bcif2cif(bcif_file)
    doc = gemmi.cif.read_string(cif_content)
    block = doc.sole_block()
    return gemmi.make_structure_from_block(block)


def _initialize_dictionary_api(containers) -> DictionaryApi:
    dict_local = user_cache_root_dir() / "mmcif_pdbx_v5_next.dic"
    if not dict_local.exists():
        dict_url = "https://raw.githubusercontent.com/wwpdb-dictionaries/mmcif_pdbx/master/dist/mmcif_pdbx_v5_next.dic"
        logger.info("Downloading mmcif dictionary from %s to %s", dict_url, dict_local)
        dict_local.parent.mkdir(parents=True, exist_ok=True)
        with dict_local.open("wb") as f, urlopen(dict_url) as response:  # noqa: S310 url is hardcoded and https
            f.write(response.read())
    return DictionaryApi(containerList=containers, consolidate=True)


def structure2bcif(structure: gemmi.Structure, bcif_file: Path):
    """Write a gemmi Structure object to a binary CIF (bcif) file.

    This is slower than other formats because gemmi does not support writing bcif files directly.
    So we convert it to a cif string first using gemmi and then convert cif to bcif using mmcif package.

    Args:
        structure: The gemmi Structure object to write.
        bcif_file: Path to the output binary CIF file.
    """
    doc = structure.make_mmcif_document(gemmi.MmcifOutputGroups(True, chem_comp=False))
    containers = []
    with StringIO(doc.as_string()) as sio:
        reader = PdbxReader(sio)
        reader.read(containers)
    dict_api = _initialize_dictionary_api(containers)
    writer = BinaryCifWriter(dictionaryApi=dict_api)
    writer.serialize(str(bcif_file), containers)


def structure2cifgz(structure: gemmi.Structure) -> bytes:
    """Render a gemmi Structure as gzipped mmCIF bytes.

    Args:
        structure: The gemmi Structure object to render.

    Returns:
        Gzipped mmCIF bytes.
    """
    doc = structure.make_mmcif_document(gemmi.MmcifOutputGroups(True, chem_comp=False))
    return gzip.compress(doc.as_string().encode("utf-8"))


def gunzip_file(gz_file: Path, output_file: Path | None = None, keep_original: bool = True) -> Path:
    """Unzip a .gz file.

    Args:
        gz_file: Path to the .gz file.
        output_file: Optional path to the output unzipped file. If None, the .gz suffix is removed from gz_file.
        keep_original: Whether to keep the original .gz file. Default is True.

    Returns:
        Path to the unzipped file.

    Raises:
        ValueError: If output_file is None and gz_file does not end with .gz.
    """
    if output_file is None and not gz_file.name.endswith(".gz"):
        msg = f"If output_file is not provided, {gz_file} must end with .gz"
        raise ValueError(msg)
    out_file = output_file or gz_file.with_suffix("")
    with gzip.open(gz_file, "rb") as f_in, out_file.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_original:
        gz_file.unlink()
    return out_file


def structure2bcifgz(structure: gemmi.Structure, bcif_gz_file: Path):
    """Write a gemmi Structure object to a binary CIF gzipped (bcif.gz) file.

    This is slower than other formats because gemmi does not support writing bcif files directly.
    So we convert it to a cif string first using gemmi and then convert cif to bcif using mmcif package.
    Finally, we gzip the bcif file.

    Args:
        structure: The gemmi Structure object to write.
        bcif_gz_file: Path to the output binary CIF gzipped file.
    """
    with tempfile.NamedTemporaryFile(suffix=".bcif", delete=True) as tmp_bcif:
        tmp_path = Path(tmp_bcif.name)
        structure2bcif(structure, tmp_path)
        with tmp_path.open("rb") as f_in, gzip.open(bcif_gz_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


Pdb2UniprotMapping = dict[str, set[tuple[str, str]]]
"""Dictionary mapping PDB ID to set of tuples containing chain and UniProt accession.

For example, `{'1abc': {('A', 'P12345'), ('B', 'Q67890')}}`.

As mapping can be fetched from Uniprot SPARQL API, which uses uppercase PDB IDs, the keys are lowercase.
"""


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
        A tuple of the input file and the output file.
    """
    for input_file in input_files:
        output_file = convert_to_cif_file(
            input_file, output_dir, copy_method, output_format=output_format, pdb2uniprot=pdb2uniprot
        )
        yield input_file, output_file


def structure_to_uniprot(structure: gemmi.Structure) -> Pdb2UniprotMapping:
    # Via sifts_unp_acc
    chain_uniprots: set[tuple[str, str]] = set()
    for entity in structure.entities:
        uniprots = entity.sifts_unp_acc
        subchains = entity.subchains
        for subchain in subchains:
            for uniprot in uniprots:
                chain_uniprots.add((subchain, uniprot))

    # Gemmi python does not expose entity->dbrefs so we have to go through the mmCIF block
    # Via struct_ref + struct_ref_seq mmcif blocks
    block = structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
    struct_ref = block.get_mmcif_category("_struct_ref.")
    db_name_unp_idx = [i for i, db_name in enumerate(struct_ref["db_name"]) if db_name == "UNP"]
    uniprot_accessions = [struct_ref["pdbx_db_accession"][i] for i in db_name_unp_idx]
    struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
    for i, u in enumerate(struct_ref_seq["pdbx_db_accession"]):
        if u in uniprot_accessions:
            chain_uniprots.add((struct_ref_seq["pdbx_strand_id"][i], u))

    return {structure.name: chain_uniprots}


def _append_uniprot_to_structure(
    structure: gemmi.Structure, chain_uniprot_pairs: set[tuple[str, str]]
) -> gemmi.Structure:
    if not chain_uniprot_pairs:
        return structure
    # sifts_unp_acc is read-only, so we have to inject using mmcif blocks
    block = structure.make_mmcif_block()
    struct_ref = block.get_mmcif_category("_struct_ref.")
    struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
    chain2entity_id: dict[str, str] = {
        chain: entity.name for entity in structure.entities for chain in entity.subchains
    }
    # TODO in chain_uniprot_pairs also include uniprot residues that match aka `A=112-210`
    # so db_align_beg and db_align_end are set so sequence identity can be computed
    fillable_struct_ref_seq_cols = {"align_id", "ref_id", "pdbx_strand_id", "pdbx_db_accession"}
    for chain, uniprot_accession in chain_uniprot_pairs:
        entity_id = chain2entity_id[chain]
        new_id = str(len(struct_ref["id"]) + 1)
        struct_ref["id"].append(new_id)
        struct_ref["entity_id"].append(entity_id)
        struct_ref["db_name"].append("UNP")
        # do not have uniprot id so use None
        struct_ref["db_code"].append(None)
        struct_ref["pdbx_db_accession"].append(uniprot_accession)
        struct_ref["pdbx_db_isoform"].append(None)

        new_seq_id = len(struct_ref_seq["align_id"]) + 1
        struct_ref_seq["align_id"].append(new_seq_id)
        struct_ref_seq["ref_id"].append(new_seq_id)
        struct_ref_seq["pdbx_strand_id"].append(chain)
        struct_ref_seq["pdbx_db_accession"].append(uniprot_accession)
        for col in struct_ref_seq:
            if col not in fillable_struct_ref_seq_cols:
                struct_ref_seq[col].append(None)

    block.set_mmcif_category("_struct_ref.", struct_ref)
    block.set_mmcif_category("_struct_ref_seq.", struct_ref_seq)
    return gemmi.make_structure_from_block(block)


# TODO move to structure.py
def add_uniprot_accessions2structure(
    structure: gemmi.Structure, pdb2uniprot: Pdb2UniprotMapping | None
) -> gemmi.Structure:
    """Add UniProt accessions to a structure if they are missing, based on the provided pdb2uniprot mapping.

    If structure has uniprot accesion that is not in `pdb2uniprot`, it will be left unchanged.

    Args:
        structure: The gemmi Structure object to add UniProt accessions to.
        pdb2uniprot: Dictionary mapping PDB ID to set of tuples containing chain and UniProt accession.
            If provided, will be used to inject UniProt accessions into the structure if they are missing.
            If None, the structure is returned unchanged.

    Returns:
        A gemmi Structure object with UniProt accessions added if they were missing
        or the structure unchanged if all accessions were already present.
    """
    if not pdb2uniprot:
        return structure
    pdb_id = structure.name
    if pdb_id not in pdb2uniprot:
        logger.warning(
            "PDB ID %s not found in pdb2uniprot mapping. Leaving structure unverified and unchanged.", pdb_id
        )
        return structure

    known = structure_to_uniprot(structure)
    missing = pdb2uniprot[pdb_id] - known[pdb_id]
    if not missing:
        # Items in pdb2uniprot are already in structure -> no changes needed
        return structure

    if not missing:
        return structure
    if known[pdb_id]:
        logger.warning(
            "Structure %s has some UniProt accessions that do not match the provided mapping. "
            "Existing: %s, Expected: %s, Missing: %s. Injecting missing accessions.",
            pdb_id,
            known[pdb_id],
            pdb2uniprot[pdb_id],
            missing,
        )
    else:
        logger.info("Injecting UniProt accessions into structure %s: %s", structure.name, missing)
    return _append_uniprot_to_structure(structure, missing)


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
        ValueError: If the requested output format is not supported.
    """
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


def split_name_and_extension(name: str) -> tuple[str, str]:
    """Split a filename into its name and extension.

    `.gz` is considered part of the extension if present.

    Examples:
        Some example usages.

        >>> from protein_quest.pdbe.io import split_name_and_extension
        >>> split_name_and_extension("1234.pdb")
        ('1234', '.pdb')
        >>> split_name_and_extension("1234.pdb.gz")
        ('1234', '.pdb.gz')

    Args:
        name: The filename to split.

    Returns:
        A tuple containing the name and the extension.
    """
    ext = ""
    if name.endswith(".gz"):
        ext = ".gz"
        name = name.removesuffix(".gz")
    i = name.rfind(".")
    if 0 < i < len(name) - 1:
        ext = name[i:] + ext
        name = name[:i]
    return name, ext


def locate_structure_file(root: Path, pdb_id: str) -> Path:
    """Locate a structure file for a given PDB ID in the specified directory.

    Uses [StructureFileExtensions][protein_quest.io.StructureFileExtensions] as potential extensions.
    Also tries different casing of the PDB ID.

    Args:
        root: The root directory to search in.
        pdb_id: The PDB ID to locate.

    Returns:
        The path to the located structure file.

    Raises:
        FileNotFoundError: If no structure file is found for the given PDB ID.
    """
    for ext in valid_structure_file_extensions:
        candidates = (
            root / f"{pdb_id}{ext}",
            root / f"{pdb_id.lower()}{ext}",
            root / f"{pdb_id.upper()}{ext}",
            root / f"pdb{pdb_id.lower()}{ext}",
        )
        for candidate in candidates:
            if candidate.exists():
                return candidate
    msg = f"No structure file found for {pdb_id} in {root}"
    raise FileNotFoundError(msg)


def glob_structure_files(input_dir: Path) -> Generator[Path]:
    """Glob for structure files in a directory.

    Uses [StructureFileExtensions][protein_quest.io.StructureFileExtensions] as valid extensions.
    Does not search recursively.

    Args:
        input_dir: The input directory to search for structure files.

    Yields:
        Paths to the found structure files.
    """
    for ext in valid_structure_file_extensions:
        yield from input_dir.glob(f"*{ext}")
