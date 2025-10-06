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

# TODO remove once v0.7.4 of gemmi is released,
# as uv pip install git+https://github.com/project-gemmi/gemmi.git installs 0.7.4.dev0 which does not print leaks
# Swallow gemmi leaked function warnings
gemmi.set_leak_warnings(False)


StructureFileExtensions = Literal[".pdb", ".pdb.gz", ".ent", ".ent.gz", ".cif", ".cif.gz", ".bcif", ".bcif.gz"]
"""Type of supported structure file extensions."""
valid_structure_file_extensions: set[str] = set(get_args(StructureFileExtensions))
"""Set of valid structure file extensions."""


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
        doc = structure.make_mmcif_document(gemmi.MmcifOutputGroups(True, chem_comp=False))
        cif_str = doc.as_string()
        with gzip.open(path, "wt") as f:
            f.write(cif_str)
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
            See [StructureFileExtensions][protein_quest.io.StructureFileExtensions]
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
