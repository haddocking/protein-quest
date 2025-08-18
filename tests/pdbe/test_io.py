from pathlib import Path

import gemmi
import pytest

from protein_quest.pdbe.io import nr_residues_in_chain, write_single_chain_pdb_file


@pytest.fixture
def cif_path() -> Path:
    return Path(__file__).parent / "fixtures" / "2y29.cif"


def test_write_single_chain_pdb_file_happypath(cif_path: Path, tmp_path: Path):
    output_file = write_single_chain_pdb_file(
        input_file=cif_path,
        chain2keep="A",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is not None
    assert output_file.name == "2y29_A2Z.cif"
    assert output_file.exists()
    structure = gemmi.read_structure(str(output_file))
    assert len(structure) == 1  # One model
    model = structure[0]
    assert len(model) == 1  # One chain
    chain = model[0]
    assert chain.name == "Z"
    assert len(chain) == 6  # 6 residues in chain Z


def test_write_single_chain_pdb_file_unknown_chain(cif_path: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture):
    output_file = write_single_chain_pdb_file(
        cif_path,
        chain2keep="B",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is None
    assert "Chain B not found in" in caplog.text


def test_nr_residues_in_chain(cif_path: Path):
    residue_count = nr_residues_in_chain(cif_path, chain="A")

    assert residue_count == 8


def test_nr_residues_in_chain_wrongchain(cif_path: Path, caplog):
    residue_count = nr_residues_in_chain(cif_path, chain="Z")

    assert residue_count == 0
    assert "Chain Z not found in" in caplog.text

