import logging
from pathlib import Path

import gemmi
import pytest
from cattrs import ClassValidationError
from orjson import JSONDecodeError

from protein_quest.__version__ import __version__
from protein_quest.pdbe.fetch import sync_fetch
from protein_quest.structure.chains import (
    CHAIN_PROVENANCE_SOFTWARE_NAME,
    ChainIdSystem,
    find_chain_in_structure,
    get_label2auth_chains,
    label_auth_mismatch,
    nr_of_residues_in_total,
    nr_residues_in_chain,
    resolve_chain_id_to_label,
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
    ("wanted", "system"),
    [
        ("B", "auth"),
        ("A", "label"),
    ],
)
def test_find_chain_in_structure(cif_8rw8: Path, wanted: str, system: ChainIdSystem):
    structure = read_structure(cif_8rw8)

    found_chain = find_chain_in_structure(structure, wanted, chain_system=system)

    assert found_chain is not None
    # Gemmi chain names are auth system for 8rw8, while the lookup input is label.
    assert found_chain.name == "B"


@pytest.mark.parametrize(
    ("cif_fixture", "expected"),
    [
        # TODO write_single_chain_structure_file should update label_auth_chains, which it does not now, expected {A:A} when write works correctly
        pytest.param("sample_cif", {"B": "A"}, id="3JRS_B2A"),
        pytest.param("sample2_cif", {"A": "A"}, id="2Y29"),
        pytest.param("af_cif", {"A": "A"}, id="AF-A0A0C5B5G6-F1"),
        pytest.param("nmr_cif", {"A": "A"}, id="1AMB"),
        pytest.param("em_cif", {"A": "A"}, id="8W77"),
        pytest.param("sample_multispan_cif", {"A": "A"}, id="6O5I"),
        pytest.param(
            "multi_accession_cif",
            {"A": "A", "B": "B", "C": "N", "D": "F", "E": "J"},
            id="1A02",
        ),
        pytest.param("multi_accession_chain_cif", {"A": "A"}, id="1UN5"),
        pytest.param(
            "multi_entity_cif",
            {
                "A": "I",
                "B": "J",
                "C": "A",
                "D": "B",
                "E": "C",
                "F": "D",
                "G": "E",
                "H": "F",
                "I": "G",
                "J": "H",
            },
            id="1F66",
        ),
        pytest.param("cif_8rw8", {"A": "B"}, id="8rw8"),
    ],
)
def test_get_label2auth_chains(cif_fixture: str, expected: dict[str, str], request: pytest.FixtureRequest):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)

    label2auth_chains = get_label2auth_chains(structure)

    assert label2auth_chains == expected


@pytest.mark.parametrize(
    ("input_chain", "chain_system", "expected"),
    [
        pytest.param("B", "auth", "A", id="auth-to-label"),
        pytest.param("A", "label", "A", id="label-passthrough"),
    ],
)
def test_resolve_chain_id_to_label(
    cif_8rw8: Path,
    input_chain: str,
    chain_system: ChainIdSystem,
    expected: str,
):
    structure = read_structure(cif_8rw8)

    resolved_chain = resolve_chain_id_to_label(
        structure,
        input_chain,
        chain_system=chain_system,
        source_file=cif_8rw8,
    )

    assert resolved_chain == expected


@pytest.mark.parametrize(
    ("chains_map", "expected"),
    [
        pytest.param({"A": "A", "B": "B"}, False, id="no-mismatch"),
        pytest.param({"A": "I", "B": "B"}, True, id="has-mismatch"),
        pytest.param({}, False, id="empty-map"),
    ],
)
def test_label_auth_mismatch(chains_map: dict[str, str], expected: bool):
    assert label_auth_mismatch(chains_map) is expected


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
