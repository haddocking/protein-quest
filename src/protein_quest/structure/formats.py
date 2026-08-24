"""Structure format read/write and conversion helpers."""

import gzip
import logging
import shutil
import tempfile
from io import StringIO
from pathlib import Path
from typing import Any
from urllib.request import urlopen

import gemmi
from mmcif.api.DictionaryApi import DictionaryApi
from mmcif.io.BinaryCifReader import BinaryCifReader
from mmcif.io.BinaryCifWriter import BinaryCifWriter
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter

from protein_quest.structure.sifts import uniprot_chain_mappings_to_cif
from protein_quest.uniprot_chains import UniprotChainMappings
from protein_quest.utils import user_cache_root_dir

logger = logging.getLogger(__name__)


def _add_em_3d_reconstruction(structure: gemmi.Structure, block: gemmi.cif.Block):
    """Append EM reconstruction resolution metadata to an mmCIF block.

    This is needed because Gemmi reads EM resolution from
    ``_em_3d_reconstruction.resolution`` into ``Structure.resolution``, but the
    generated mmCIF document does not recreate that category automatically.

    Args:
        structure: Structure whose EM resolution metadata may need to be
            restored to the output mmCIF block.
        block: mmCIF block produced for ``structure`` that should receive the
            ``_em_3d_reconstruction`` fields when they are missing.
    """
    if structure.resolution <= 0 or not _is_em_method(structure):
        return
    if block.get_mmcif_category("_em_3d_reconstruction."):
        msg = (
            "Gemmi has support for _em_3d_reconstruction, please create "
            "issue to remove _add_em_3d_reconstruction function"
        )
        raise ValueError(msg)
    block.set_pair("_em_3d_reconstruction.entry_id", structure.name)
    block.set_pair("_em_3d_reconstruction.id", "1")
    block.set_pair("_em_3d_reconstruction.resolution", str(structure.resolution))


def _is_em_method(structure: gemmi.Structure) -> bool:
    # Could have used ./metadata.py::_structure_method, but do not to prevent circular import
    try:
        experimental_method = structure.info["_exptl.method"].lower()
    except KeyError:
        experimental_method = ""
    return "electron microscopy" in experimental_method


def _make_mmcif_document(
    structure: gemmi.Structure,
    *,
    uniprot_chain_mappings: UniprotChainMappings | None = None,
) -> gemmi.cif.Document:
    """Create an mmCIF document and preserve EM resolution metadata when needed."""
    doc = structure.make_mmcif_document()
    block = doc.sole_block()

    _add_em_3d_reconstruction(structure, block)
    if uniprot_chain_mappings:
        uniprot_chain_mappings_to_cif(uniprot_chain_mappings, block)
    return doc


def write_structure(
    structure: gemmi.Structure,
    path: Path,
    uniprot_chain_mappings: UniprotChainMappings | None = None,
):
    """Write a gemmi structure to a file.

    Args:
        structure: The gemmi structure to write.
        path: The file path to write the structure to.
            The format depends on the file extension.
            See [StructureFileExtensions][protein_quest.structure.types.StructureFileExtensions]
            for supported extensions.
        uniprot_chain_mappings: Optional UniProt chain mappings in label system to
            emit as ``_pdbx_sifts_unp_segments`` for mmCIF outputs.

    Raises:
        ValueError: If the file extension is not supported."""
    if path.name.endswith(".pdb") or path.name.endswith(".ent"):
        body: str = structure.make_pdb_string()
        path.write_text(body)
    elif path.name.endswith(".pdb.gz") or path.name.endswith(".ent.gz"):
        body = structure.make_pdb_string()
        with gzip.open(path, "wt") as f:
            f.write(body)
    elif path.name.endswith(".cif"):
        doc = _make_mmcif_document(structure, uniprot_chain_mappings=uniprot_chain_mappings)
        doc.write_file(str(path))
    elif path.name.endswith(".cif.gz"):
        path.write_bytes(structure2cifgz(structure, uniprot_chain_mappings=uniprot_chain_mappings))
    elif path.name.endswith(".bcif"):
        structure2bcif(structure, path, uniprot_chain_mappings=uniprot_chain_mappings)
    elif path.name.endswith(".bcif.gz"):
        structure2bcifgz(structure, path, uniprot_chain_mappings=uniprot_chain_mappings)
    else:
        msg = f"Unsupported file extension in {path.name}."
        raise ValueError(msg)


def read_structure(file: Path) -> gemmi.Structure:
    """Read a structure from a file.

    Args:
        file: Path to the input structure file.
            See [StructureFileExtensions][protein_quest.structure.types.StructureFileExtensions]
            for supported extensions.

    Returns:
        A gemmi Structure object representing the structure in the file."""
    if file.name.endswith(".bcif"):
        return bcif2structure(file)
    if file.name.endswith(".bcif.gz"):
        return bcifgz2structure(file)
    # TODO reading 3c5f_updated.cif.gz gives `RuntimeError: _pdbx_sifts_xref_db: seq_id not found: 230`
    return gemmi.read_structure(str(file))


def _bcif_to_cif_block(file: Path, *, source_file: Path | None = None) -> gemmi.cif.Block | None:
    source = source_file or file
    cif_content = bcif2cif(file)
    if not cif_content:
        logger.info("Unable to read mmCIF block from %s: empty CIF content after BCIF conversion", source)
        return None
    try:
        doc = gemmi.cif.read_string(cif_content)
    except ValueError as exc:
        logger.info("Unable to read mmCIF block from %s: %s", source, exc)
        return None
    return doc.sole_block() if len(doc) else None


def _bcif_gz_to_cif_block(file: Path) -> gemmi.cif.Block | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".bcif", delete=True) as tmp_bcif:
            tmp_path = Path(tmp_bcif.name)
            gunzip_file(file, output_file=tmp_path, keep_original=True)
            return _bcif_to_cif_block(tmp_path, source_file=file)
    except (OSError, EOFError) as exc:
        logger.info("Unable to read mmCIF block from %s: %s", file, exc)
        return None


def read_structure_as_cif_block(file: Path) -> gemmi.cif.Block | None:
    """Read a structure file as a raw mmCIF block when possible.

    Useful because Gemmi structure objects do not preserve all mmCIF categories.

    Args:
        file: Path to an input file intended to be readable by
            ``gemmi.cif.read_file``.

    Returns:
        The sole mmCIF block from the file, or ``None`` when the file cannot be
        parsed as mmCIF or contains no blocks.
    """
    if file.name.endswith(".bcif"):
        return _bcif_to_cif_block(file)
    if file.name.endswith(".bcif.gz"):
        return _bcif_gz_to_cif_block(file)
    try:
        doc = gemmi.cif.read_file(str(file))
    except ValueError as exc:
        logger.info("Unable to read mmCIF block from %s: %s", file, exc)
        return None
    if len(doc) == 0:
        return None
    return doc.sole_block()


def bcif2cif(bcif_file: Path) -> str:
    """Convert a binary CIF (bcif) file to a CIF string.

    Args:
        bcif_file: Path to the binary CIF file.

    Returns:
        A string containing the CIF representation of the structure."""
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
        A gemmi Structure object representing the structure in the bcif.gz file."""
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
        A gemmi Structure object representing the structure in the bcif file."""
    cif_content = bcif2cif(bcif_file)
    doc = gemmi.cif.read_string(cif_content)
    block = doc.sole_block()
    return gemmi.make_structure_from_block(block)


def _initialize_dictionary_api(containers: list[Any]) -> DictionaryApi:
    dict_local = user_cache_root_dir() / "mmcif_pdbx_v5_next.dic"
    if not dict_local.exists():
        dict_url = "https://raw.githubusercontent.com/wwpdb-dictionaries/mmcif_pdbx/master/dist/mmcif_pdbx_v5_next.dic"
        logger.info("Downloading mmcif dictionary from %s to %s", dict_url, dict_local)
        dict_local.parent.mkdir(parents=True, exist_ok=True)
        with dict_local.open("wb") as f, urlopen(dict_url) as response:
            f.write(response.read())
    return DictionaryApi(containerList=containers, consolidate=True)


def structure2bcif(
    structure: gemmi.Structure,
    bcif_file: Path,
    uniprot_chain_mappings: UniprotChainMappings | None = None,
):
    """Write a gemmi Structure object to a binary CIF (bcif) file.

    This is slower than other formats because gemmi does not support writing bcif files directly.
    So we convert it to a cif string first using gemmi and then convert cif to bcif using mmcif package.

    Args:
        structure: The gemmi Structure object to write.
        bcif_file: Path to the output binary CIF file.
        uniprot_chain_mappings: Optional UniProt chain mappings in label system to
            emit as ``_pdbx_sifts_unp_segments`` for mmCIF outputs."""
    doc = _make_mmcif_document(structure, uniprot_chain_mappings=uniprot_chain_mappings)
    containers: list[Any] = []
    with StringIO(doc.as_string()) as sio:
        reader = PdbxReader(sio)
        reader.read(containers)
    dict_api = _initialize_dictionary_api(containers)
    writer = BinaryCifWriter(dictionaryApi=dict_api)
    writer.serialize(str(bcif_file), containers)


def structure2cifgz(
    structure: gemmi.Structure,
    uniprot_chain_mappings: UniprotChainMappings | None = None,
) -> bytes:
    """Render a gemmi Structure as gzipped mmCIF bytes.

    Args:
        structure: The gemmi Structure object to render.
        uniprot_chain_mappings: Optional UniProt chain mappings in label system to
            emit as ``_pdbx_sifts_unp_segments`` for mmCIF outputs.

    Returns:
        Gzipped mmCIF bytes."""
    doc = _make_mmcif_document(structure, uniprot_chain_mappings=uniprot_chain_mappings)
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
        ValueError: If output_file is None and gz_file does not end with .gz."""
    if output_file is None and not gz_file.name.endswith(".gz"):
        msg = f"If output_file is not provided, {gz_file} must end with .gz"
        raise ValueError(msg)
    out_file = output_file or gz_file.with_suffix("")
    with gzip.open(gz_file, "rb") as f_in, out_file.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_original:
        gz_file.unlink()
    return out_file


def structure2bcifgz(
    structure: gemmi.Structure,
    bcif_gz_file: Path,
    uniprot_chain_mappings: UniprotChainMappings | None = None,
):
    """Write a gemmi Structure object to a binary CIF gzipped (bcif.gz) file.

    This is slower than other formats because gemmi does not support writing bcif files directly.
    So we convert it to a cif string first using gemmi and then convert cif to bcif using mmcif package.
    Finally, we gzip the bcif file.

    Args:
        structure: The gemmi Structure object to write.
        bcif_gz_file: Path to the output binary CIF gzipped file.
        uniprot_chain_mappings: Optional UniProt chain mappings in label system to
            emit as ``_pdbx_sifts_unp_segments`` for mmCIF outputs."""
    with tempfile.NamedTemporaryFile(suffix=".bcif", delete=True) as tmp_bcif:
        tmp_path = Path(tmp_bcif.name)
        structure2bcif(structure, tmp_path, uniprot_chain_mappings=uniprot_chain_mappings)
        with tmp_path.open("rb") as f_in, gzip.open(bcif_gz_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
