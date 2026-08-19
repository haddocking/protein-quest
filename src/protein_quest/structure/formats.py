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

from protein_quest.utils import user_cache_root_dir

logger = logging.getLogger(__name__)


def _sifts_polymer_rows(structure: gemmi.Structure) -> dict[str, list[str]]:
    """Build minimal SIFTS residue rows for mmCIF output.

    Gemmi reconstructs ``Entity.sifts_unp_acc`` and ``Residue.sifts_unp`` from
    ``_pdbx_sifts_xref_db`` when reading mmCIF. This helper only emits the
    subset of columns needed for that roundtrip.

    Args:
        structure: Structure whose polymer residues may contain SIFTS UniProt
            annotations.

    Returns:
        Mapping of ``_pdbx_sifts_xref_db`` column names to row values. Returns
        empty column lists when the structure has no SIFTS residue annotations
        that can be written safely.
    """
    rows: dict[str, list[str]] = {
        "entity_id": [],
        "asym_id": [],
        "seq_id_ordinal": [],
        "seq_id": [],
        "observed": [],
        "unp_res": [],
        "unp_num": [],
        "unp_acc": [],
    }
    subchain_to_entity = {
        subchain: entity for entity in structure.entities for subchain in entity.subchains if entity.sifts_unp_acc
    }

    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            if not polymer:
                continue
            subchain_id = polymer.subchain_id()
            entity = subchain_to_entity.get(subchain_id)
            if entity is None:
                continue

            for residue in polymer:
                unp_res, acc_index, unp_num = residue.sifts_unp
                if unp_num <= 0 or not unp_res:
                    continue
                rows["entity_id"].append(entity.name)
                rows["asym_id"].append(subchain_id)
                # Gemmi only consumes rows with seq_id_ordinal == 1.
                rows["seq_id_ordinal"].append("1")
                rows["seq_id"].append(str(residue.label_seq))
                rows["observed"].append("y")
                rows["unp_res"].append(unp_res)
                rows["unp_num"].append(str(unp_num))
                rows["unp_acc"].append(entity.sifts_unp_acc[acc_index])

    # TODO also write _pdbx_sifts_unp_segments as uniprot_chain_mappings_from_cif reads those
    return rows


def _add_sifts_xref_db(structure: gemmi.Structure, block: gemmi.cif.Block):
    """Append minimal SIFTS residue annotations to an mmCIF block.

    This is needed because the mmCIF document produced by Gemmi is missing the
    ``_pdbx_sifts_xref_db`` category, even when the input structure still has
    SIFTS UniProt annotations in memory.

    Args:
        structure: Structure that may contain in-memory SIFTS annotations on
            entities and residues.
        block: mmCIF block produced for ``structure`` that should receive the
            ``_pdbx_sifts_xref_db`` category when SIFTS rows are available.
    """
    rows = _sifts_polymer_rows(structure)
    if not rows["entity_id"]:
        return
    block.set_mmcif_category("_pdbx_sifts_xref_db.", rows)


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


def _make_mmcif_document(structure: gemmi.Structure) -> gemmi.cif.Document:
    """Create an mmCIF document and preserve EM resolution metadata when needed."""
    doc = structure.make_mmcif_document()
    block = doc.sole_block()

    _add_em_3d_reconstruction(structure, block)
    _add_sifts_xref_db(structure, block)
    # TODO set _chem_comp.type from current always False to correct value, see
    # https://github.com/project-gemmi/gemmi/discussions/362
    # and
    # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.type.html
    return doc


def write_structure(structure: gemmi.Structure, path: Path):
    """Write a gemmi structure to a file.

    Args:
        structure: The gemmi structure to write.
        path: The file path to write the structure to.
            The format depends on the file extension.
            See [StructureFileExtensions][protein_quest.structure.types.StructureFileExtensions]
            for supported extensions.

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
        doc = _make_mmcif_document(structure)
        doc.write_file(str(path))
    elif path.name.endswith(".cif.gz"):
        path.write_bytes(structure2cifgz(structure))
    elif path.name.endswith(".bcif"):
        structure2bcif(structure, path)
    elif path.name.endswith(".bcif.gz"):
        structure2bcifgz(structure, path)
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


def structure2bcif(structure: gemmi.Structure, bcif_file: Path):
    """Write a gemmi Structure object to a binary CIF (bcif) file.

    This is slower than other formats because gemmi does not support writing bcif files directly.
    So we convert it to a cif string first using gemmi and then convert cif to bcif using mmcif package.

    Args:
        structure: The gemmi Structure object to write.
        bcif_file: Path to the output binary CIF file."""
    doc = _make_mmcif_document(structure)
    containers: list[Any] = []
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
        Gzipped mmCIF bytes."""
    doc = _make_mmcif_document(structure)
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


def structure2bcifgz(structure: gemmi.Structure, bcif_gz_file: Path):
    """Write a gemmi Structure object to a binary CIF gzipped (bcif.gz) file.

    This is slower than other formats because gemmi does not support writing bcif files directly.
    So we convert it to a cif string first using gemmi and then convert cif to bcif using mmcif package.
    Finally, we gzip the bcif file.

    Args:
        structure: The gemmi Structure object to write.
        bcif_gz_file: Path to the output binary CIF gzipped file."""
    with tempfile.NamedTemporaryFile(suffix=".bcif", delete=True) as tmp_bcif:
        tmp_path = Path(tmp_bcif.name)
        structure2bcif(structure, tmp_path)
        with tmp_path.open("rb") as f_in, gzip.open(bcif_gz_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
