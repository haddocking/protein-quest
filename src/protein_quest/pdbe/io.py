"""Module for structure file input/output."""

import gzip
import logging
import shutil
from collections.abc import Generator, Iterable
from datetime import UTC, datetime
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

from protein_quest.__version__ import __version__
from protein_quest.utils import CopyMethod, copyfile, user_cache_root_dir

# TODO this modules is used outside protein_quest.pdbe, move to protein_quest.io or split?

logger = logging.getLogger(__name__)

# TODO remove once v0.7.4 of gemmi is released,
# as uv pip install git+https://github.com/project-gemmi/gemmi.git installs 0.7.4.dev0 which does not print leaks
# Swallow gemmi leaked function warnings
gemmi.set_leak_warnings(False)


StructureFileExtensions = Literal[".pdb", ".pdb.gz", ".ent", ".ent.gz", ".cif", ".cif.gz", ".bcif"]
"""Type of supported structure file extensions."""
valid_structure_file_extensions: set[str] = set(get_args(StructureFileExtensions))
"""Set of valid structure file extensions."""


def nr_residues_in_chain(file: Path | str, chain: str = "A") -> int:
    """Returns the number of residues in a specific chain from a structure file.

    Args:
        file: Path to the input structure file.
        chain: Chain to count residues of.

    Returns:
        The number of residues in the specified chain.
    """
    structure = gemmi.read_structure(str(file))
    gchain = find_chain_in_structure(structure, chain)
    if gchain is None:
        logger.warning("Chain %s not found in %s. Returning 0.", chain, file)
        return 0
    return len(gchain)


def find_chain_in_structure(structure: gemmi.Structure, wanted_chain: str) -> gemmi.Chain | None:
    for model in structure:
        chain = find_chain_in_model(model, wanted_chain)
        if chain is not None:
            return chain
    return None


def find_chain_in_model(model: gemmi.Model, wanted_chain: str) -> gemmi.Chain | None:
    chain = model.find_chain(wanted_chain)
    if chain is None:
        # For chain A in 4v92 the find_chain method returns None,
        # however it is prefixed with 'B',
        # so we try again as last char of chain name
        mchains = [c for c in model if c.name.endswith(wanted_chain)]
        if mchains:
            return mchains[0]
    return chain


def write_structure(structure: gemmi.Structure, path: Path):
    """Write a gemmi structure to a file.

    Args:
        structure: The gemmi structure to write.
        path: The file path to write the structure to.
            The format depends on the file extension.
            See [StructureFileExtensions][protein_quest.pdbe.io.StructureFileExtensions]
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
        doc = structure.make_mmcif_document(gemmi.MmcifOutputGroups(True, chem_comp=False))
        cif_str = doc.as_string()
        with gzip.open(path, "wt") as f:
            f.write(cif_str)
    elif path.name.endswith(".bcif"):
        structure2bcif(structure, path)
    else:
        msg = f"Unsupported file extension in {path.name}. Supported extensions are: {valid_structure_file_extensions}"
        raise ValueError(msg)


def read_structure(file: Path) -> gemmi.Structure:
    """Read a structure from a file.

    Args:
        file: Path to the input structure file.
            See [StructureFileExtensions][protein_quest.pdbe.io.StructureFileExtensions]
            for supported extensions.

    Returns:
        A gemmi Structure object representing the structure in the file.
    """
    if file.name.endswith(".bcif"):
        return bcif2structure(file)
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


def gunzip_file(gz_file: Path, output_file: Path | None = None, keep_original: bool = True) -> Path:
    """Unzip a .gz file.

    Args:
        gz_file: Path to the .gz file.
        output_file: Optional path to the output unzipped file. If None, the .gz suffix is removed from gz_file.
        keep_original: Whether to keep the original .gz file. Default is True.

    Returns:
        Path to the unzipped file.
    """
    if not gz_file.name.endswith(".gz"):
        return gz_file
    out_file = output_file or gz_file.with_suffix("")
    with gzip.open(gz_file, "rb") as f_in, out_file.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_original:
        gz_file.unlink()
    return out_file


def convert_to_cif_files(
    input_files: Iterable[Path], output_dir: Path, copy_method: CopyMethod
) -> Generator[tuple[Path, Path]]:
    """Convert structure files to .cif format.

    Args:
        input_files: Iterable of structure files to convert.
        output_dir: Directory to save the converted .cif files.
        copy_method: How to copy when no changes are needed to output file.

    Yields:
        A tuple of the input file and the output file.
    """
    for input_file in input_files:
        output_file = convert_to_cif_file(input_file, output_dir, copy_method)
        yield input_file, output_file


def convert_to_cif_file(input_file: Path, output_dir: Path, copy_method: CopyMethod) -> Path:
    """Convert a single structure file to .cif format.

    Args:
        input_file: The structure file to convert.
            See [StructureFileExtensions][protein_quest.pdbe.io.StructureFileExtensions]
            for supported extensions.
        output_dir: Directory to save the converted .cif file.
        copy_method: How to copy when no changes are needed to output file.

    Returns:
        Path to the converted .cif file.
    """
    name, extension = split_name_and_extension(input_file.name)
    output_file = output_dir / f"{name}.cif"
    if output_file.exists():
        logger.info("Output file %s already exists for input file %s. Skipping.", output_file, input_file)
    elif extension in {".pdb", ".pdb.gz", ".ent", ".ent.gz"}:
        structure = read_structure(input_file)
        write_structure(structure, output_file)
    elif extension == ".cif":
        logger.info("File %s is already in .cif format, copying to %s", input_file, output_dir)
        copyfile(input_file, output_file, copy_method)
    elif extension == ".cif.gz":
        gunzip_file(input_file, output_file=output_file, keep_original=True)
    elif extension == ".bcif":
        with output_file.open("w") as f:
            f.write(bcif2cif(input_file))
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

    Uses [StructureFileExtensions][protein_quest.pdbe.io.StructureFileExtensions] as potential extensions.
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

    Uses [StructureFileExtensions][protein_quest.pdbe.io.StructureFileExtensions] as valid extensions.

    Args:
        input_dir: The input directory to search for structure files.

    Yields:
        Paths to the found structure files.
    """
    for ext in valid_structure_file_extensions:
        yield from input_dir.glob(f"*{ext}")


class ChainNotFoundError(IndexError):
    """Exception raised when a chain is not found in a structure."""

    def __init__(self, chain: str, file: Path | str, available_chains: Iterable[str]):
        super().__init__(f"Chain {chain} not found in {file}. Available chains are: {available_chains}")
        self.chain_id = chain
        self.file = file


def _dedup_helices(structure: gemmi.Structure):
    helix_starts: set[str] = set()
    duplicate_helix_indexes: list[int] = []
    for hindex, helix in enumerate(structure.helices):
        if str(helix.start) in helix_starts:
            logger.debug(f"Duplicate start helix found: {hindex} {helix.start}, removing")
            duplicate_helix_indexes.append(hindex)
        else:
            helix_starts.add(str(helix.start))
    for helix_index in reversed(duplicate_helix_indexes):
        structure.helices.pop(helix_index)


def _dedup_sheets(structure: gemmi.Structure, chain2keep: str):
    duplicate_sheet_indexes: list[int] = []
    for sindex, sheet in enumerate(structure.sheets):
        if sheet.name != chain2keep:
            duplicate_sheet_indexes.append(sindex)
    for sheet_index in reversed(duplicate_sheet_indexes):
        structure.sheets.pop(sheet_index)


def _add_provenance_info(structure: gemmi.Structure, chain2keep: str, out_chain: str):
    old_id = structure.name
    new_id = structure.name + f"{chain2keep}2{out_chain}"
    structure.name = new_id
    structure.info["_entry.id"] = new_id
    new_title = f"From {old_id} chain {chain2keep} to {out_chain}"
    structure.info["_struct.title"] = new_title
    structure.info["_struct_keywords.pdbx_keywords"] = new_title.upper()
    new_si = gemmi.SoftwareItem()
    new_si.classification = gemmi.SoftwareItem.Classification.DataExtraction
    new_si.name = "protein-quest.pdbe.io.write_single_chain_pdb_file"
    new_si.version = str(__version__)
    new_si.date = str(datetime.now(tz=UTC).date())
    structure.meta.software = [*structure.meta.software, new_si]


def chains_in_structure(structure: gemmi.Structure) -> set[gemmi.Chain]:
    """Get a list of chains in a structure."""
    return {c for model in structure for c in model}


# TODO rename to write_single_chain_structure_file
def write_single_chain_pdb_file(
    input_file: Path,
    chain2keep: str,
    output_dir: Path,
    out_chain: str = "A",
    copy_method: CopyMethod = "copy",
) -> Path:
    """Write a single chain from a structure file to a new structure file.

    Also

    - removes ligands and waters
    - renumbers atoms ids
    - removes chem_comp section from cif files
    - adds provenance information to the header like software and input file+chain

    This function is equivalent to the following gemmi commands:

    ```shell
    gemmi convert --remove-lig-wat --select=B --to=cif chain-in/3JRS.cif - | \\
    gemmi convert --from=cif --rename-chain=B:A - chain-out/3JRS_B2A.gemmi.cif
    ```

    Args:
        input_file: Path to the input structure file.
        chain2keep: The chain to keep.
        output_dir: Directory to save the output file.
        out_chain: The chain identifier for the output file.
        copy_method: How to copy when no changes are needed to output file.

    Returns:
        Path to the output structure file

    Raises:
        FileNotFoundError: If the input file does not exist.
        ChainNotFoundError: If the specified chain is not found in the input file.
    """

    logger.debug(f"chain2keep: {chain2keep}, out_chain: {out_chain}")
    structure = read_structure(input_file)
    structure.setup_entities()

    chain = find_chain_in_structure(structure, chain2keep)
    chainnames_in_structure = {c.name for c in chains_in_structure(structure)}
    if chain is None:
        raise ChainNotFoundError(chain2keep, input_file, chainnames_in_structure)
    chain_name = chain.name
    name, extension = split_name_and_extension(input_file.name)
    output_file = output_dir / f"{name}_{chain_name}2{out_chain}{extension}"

    if output_file.exists():
        logger.info("Output file %s already exists for input file %s. Skipping.", output_file, input_file)
        return output_file

    if chain_name == out_chain and len(chainnames_in_structure) == 1:
        logger.info(
            "%s only has chain %s and out_chain is also %s. Copying file to %s.",
            input_file,
            chain_name,
            out_chain,
            output_file,
        )
        copyfile(input_file, output_file, copy_method)
        return output_file

    gemmi.Selection(chain_name).remove_not_selected(structure)
    for m in structure:
        m.remove_ligands_and_waters()
    structure.setup_entities()
    structure.rename_chain(chain_name, out_chain)
    _dedup_helices(structure)
    _dedup_sheets(structure, out_chain)
    _add_provenance_info(structure, chain_name, out_chain)

    write_structure(structure, output_file)

    return output_file
