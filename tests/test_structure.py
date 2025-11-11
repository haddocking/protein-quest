import logging
from pathlib import Path

import gemmi
import pytest
from platformdirs import user_cache_dir

from protein_quest.io import read_structure
from protein_quest.pdbe.fetch import sync_fetch
from protein_quest.structure import (
    ChainNotFoundError,
    nr_residues_in_chain,
    structure2uniprot_accessions,
    write_single_chain_structure_file,
)

logger = logging.getLogger(__name__)


def test_write_single_chain_structure_file_happypath(sample2_cif: Path, tmp_path: Path):
    output_file = write_single_chain_structure_file(
        input_file=sample2_cif,
        chain2keep="A",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is not None
    assert output_file.name == "2Y29_A2Z.cif.gz"
    assert output_file.exists()
    structure = read_structure(output_file)
    assert len(structure) == 1  # One model
    model = structure[0]
    assert len(model) == 1  # One chain
    chain = model[0]
    assert chain.name == "Z"
    assert len(chain) == 6  # 6 residues in chain Z


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
        chain2keep="A",
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
    assert len(chain) == 59


def test_nr_residues_in_chain(sample2_cif: Path):
    residue_count = nr_residues_in_chain(sample2_cif, chain="A")

    assert residue_count == 8


def test_nr_residues_in_chain_wrongchain(sample2_cif: Path, caplog):
    residue_count = nr_residues_in_chain(sample2_cif, chain="Z")

    assert residue_count == 0
    assert "Chain Z not found in" in caplog.text


def test_structure2uniprot_accessions_present(sample2_cif: Path):
    structure = read_structure(sample2_cif)
    accessions = structure2uniprot_accessions(structure)

    assert accessions == {"P05067"}


def test_structure2uniprot_accessions_missing(sample_cif: Path, caplog):
    # Empty struct_ref category to simulate missing UniProt accessions
    structure_with_unp = read_structure(sample_cif)
    block_without_struct_ref = structure_with_unp.make_mmcif_block(
        gemmi.MmcifOutputGroups(True, chem_comp=False, struct_ref=False)
    )
    structure = gemmi.make_structure_from_block(block_without_struct_ref)

    accessions = structure2uniprot_accessions(structure)

    assert accessions == set()
    assert "No UniProt accessions found in structure 3JRSB2A" in caplog.text
