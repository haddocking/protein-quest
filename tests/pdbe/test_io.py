import logging
from pathlib import Path

import gemmi
import pytest

from protein_quest.pdbe.io import (
    ChainNotFoundError,
    glob_structure_files,
    locate_structure_file,
    nr_residues_in_chain,
    write_single_chain_pdb_file,
    write_structure,
)


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


def test_write_single_chain_pdb_file_already_exists(cif_path: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture):
    fake_output_file = tmp_path / "2y29_A2Z.cif"
    fake_output_file.write_text("fake content")
    caplog.set_level(logging.INFO)

    output_file = write_single_chain_pdb_file(
        input_file=cif_path,
        chain2keep="A",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file == fake_output_file
    assert "Skipping" in caplog.text


def test_write_single_chain_pdb_file_unknown_chain(cif_path: Path, tmp_path: Path):
    with pytest.raises(ChainNotFoundError):
        write_single_chain_pdb_file(
            cif_path,
            chain2keep="B",
            output_dir=tmp_path,
            out_chain="Z",
        )


def test_write_single_chain_pdb_file_unknown_format(tmp_path: Path):
    with pytest.raises(RuntimeError, match="Unknown format"):
        write_single_chain_pdb_file(
            tmp_path / "nonexistent_file.xyz",
            chain2keep="B",
            output_dir=tmp_path,
            out_chain="Z",
        )


def test_nr_residues_in_chain(cif_path: Path):
    residue_count = nr_residues_in_chain(cif_path, chain="A")

    assert residue_count == 8


def test_nr_residues_in_chain_wrongchain(cif_path: Path, caplog):
    residue_count = nr_residues_in_chain(cif_path, chain="Z")

    assert residue_count == 0
    assert "Chain Z not found in" in caplog.text


@pytest.mark.parametrize("extension", [".pdb", ".pdb.gz", ".cif", ".cif.gz"])
def test_write_structure(cif_path: Path, tmp_path: Path, extension: str):
    structure = gemmi.read_structure(str(cif_path))
    output_file = tmp_path / f"bla{extension}"

    write_structure(structure, output_file)

    found_files = list(glob_structure_files(tmp_path))
    assert len(found_files) == 1
    assert found_files[0].name == output_file.name


def test_locate_structure_file(cif_path: Path):
    root = cif_path.parent
    result = locate_structure_file(root, "2y29")

    assert result == cif_path
