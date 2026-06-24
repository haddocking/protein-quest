import logging
from pathlib import Path

import gemmi
import pytest
from cattrs import ClassValidationError
from orjson import JSONDecodeError
from platformdirs import user_cache_dir

from protein_quest.__version__ import __version__
from protein_quest.pdbe.fetch import sync_fetch
from protein_quest.structure.chains import (
    CHAIN_PROVENANCE_SOFTWARE_NAME,
    nr_of_residues_in_total,
    nr_residues_in_chain,
    retrieve_chain_extraction_provenance,
    write_single_chain_structure_file,
)
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.formats import read_structure

logger = logging.getLogger(__name__)


def test_write_single_chain_structure_file_happypath(sample2_cif: Path, tmp_path: Path):
    input_structure = read_structure(sample2_cif)

    output_file = write_single_chain_structure_file(
        input_file=sample2_cif,
        chain2keep="A",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is not None
    assert output_file.name == "2Y29_A2Z.cif.gz"
    assert output_file.exists()

    # Chain checks
    structure = read_structure(output_file)
    assert len(structure) == 1
    model = structure[0]
    assert len(model) == 1
    chain = model[0]
    assert chain.name == "Z"
    assert len(chain) == 6

    # Unchanged main/info
    assert structure.name == input_structure.name
    assert structure.info["_entry.id"] == input_structure.info["_entry.id"]
    assert structure.info["_struct.title"] == input_structure.info["_struct.title"]

    # Added software item
    assert len(structure.meta.software) == len(input_structure.meta.software) + 1
    extracted_provenance = retrieve_chain_extraction_provenance(structure)
    assert extracted_provenance is not None
    software_item, provenance = extracted_provenance
    assert software_item.name == CHAIN_PROVENANCE_SOFTWARE_NAME
    assert software_item.version == __version__
    assert software_item.classification == gemmi.SoftwareItem.Classification.DataExtraction
    assert provenance.chain2keep == "A"
    assert provenance.out_chain == "Z"


def test_write_single_chain_structure_file_with_secondary_structure(sample_cif: Path, tmp_path: Path):
    output_file = write_single_chain_structure_file(
        input_file=sample_cif,
        chain2keep="A",
        output_dir=tmp_path,
    )
    structure = read_structure(output_file)
    assert len(structure.helices) == 4
    assert len(structure.sheets) == 1


def test_write_single_chain_structure_file_already_exists(
    sample2_cif: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    fake_output_file = tmp_path / "2Y29_A2Z.cif.gz"
    fake_output_file.write_text("fake content")
    caplog.set_level(logging.INFO)

    output_file = write_single_chain_structure_file(
        input_file=sample2_cif,
        chain2keep="A",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file == fake_output_file
    assert "Skipping" in caplog.text


def test_write_single_chain_structure_file_unknown_chain(sample2_cif: Path, tmp_path: Path):
    with pytest.raises(ChainNotFoundError):
        write_single_chain_structure_file(
            sample2_cif,
            chain2keep="B",
            output_dir=tmp_path,
            out_chain="Z",
        )


def test_write_single_chain_structure_file_unknown_format(tmp_path: Path):
    with pytest.raises(RuntimeError, match="Unknown format"):
        write_single_chain_structure_file(
            tmp_path / "nonexistent_file.xyz",
            chain2keep="B",
            output_dir=tmp_path,
            out_chain="Z",
        )


@pytest.fixture
def download_cache_dir() -> Path:
    cache_dir = Path(user_cache_dir("protein-quest-tests"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def fetch_pdb_file_cached(pdb_id: str, download_cache_dir: Path) -> Path:
    cache_fn = download_cache_dir / f"{pdb_id.lower()}.cif.gz"
    if cache_fn.exists():
        logger.info(f"[cache hit] Using cached file: {cache_fn}")
        return cache_fn

    logger.info(f"[cache miss] Downloading file for {pdb_id} to {cache_fn}")
    fetched_files = sync_fetch([pdb_id], download_cache_dir, max_parallel_downloads=1)
    fetched_file = fetched_files[pdb_id]
    assert cache_fn == fetched_file
    return fetched_file


@pytest.fixture
def lowercase_chain_cif(download_cache_dir: Path) -> Path:
    """Big model (ribosome complex) with >36 chains"""
    return fetch_pdb_file_cached("5KCS", download_cache_dir)


def test_write_single_chain_structure_file_lowercase_chain(lowercase_chain_cif: Path, tmp_path: Path):
    output_file = write_single_chain_structure_file(
        input_file=lowercase_chain_cif,
        chain2keep="l",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is not None
    structure = read_structure(output_file)
    assert len(structure) == 1  # One model
    model = structure[0]
    assert len(model) == 1  # One chain
    chain = model[0]
    assert chain.name == "Z"
    assert len(chain) == 123


@pytest.fixture
def multi_model_cif(download_cache_dir: Path) -> Path:
    """A multi model structure"""
    return fetch_pdb_file_cached("2lmc", download_cache_dir)


def test_write_single_chain_structure_file_multi_model_cif(multi_model_cif: Path, tmp_path: Path):
    output_file = write_single_chain_structure_file(
        input_file=multi_model_cif,
        chain2keep="B",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is not None
    structure = read_structure(output_file)
    assert len(structure) == 1  # One model
    model = structure[0]
    assert len(model) == 1  # One chain
    chain = model[0]
    assert chain.name == "Z"
    assert len(chain) == 61
    subchains = {chain: entity.name for entity in structure.entities for chain in entity.subchains}
    expected_subchains = {"Z": "1"}
    assert subchains == expected_subchains


def test_nr_residues_in_chain(sample2_cif: Path):
    residue_count = nr_residues_in_chain(sample2_cif, chain="A")

    assert residue_count == 8


def test_nr_residues_in_chain_wrongchain(sample2_cif: Path, caplog: pytest.LogCaptureFixture):
    residue_count = nr_residues_in_chain(sample2_cif, chain="Z")

    assert residue_count == 0
    assert "Chain Z not found in" in caplog.text


def test_nr_of_residues_in_total(sample2_cif: Path):
    structure = read_structure(sample2_cif)
    total_residues = nr_of_residues_in_total(structure)

    assert total_residues == 8


@pytest.mark.parametrize(
    "contact_author, exception_type", [("not a valid json", JSONDecodeError), ('{"a:":42}', ClassValidationError)]
)
def test_retrieve_chain_extraction_provenance_badjson(
    contact_author: str, exception_type: type, caplog: pytest.LogCaptureFixture
):
    caplog.set_level(logging.WARNING)
    structure = gemmi.Structure()
    software_item = gemmi.SoftwareItem()
    software_item.name = CHAIN_PROVENANCE_SOFTWARE_NAME
    software_item.contact_author = contact_author
    structure.meta.software = [software_item]

    with pytest.raises(exception_type):
        retrieve_chain_extraction_provenance(structure)
