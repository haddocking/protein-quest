from pathlib import Path

import gemmi
import pytest

from protein_quest.alphafold.confidence import filter_out_low_confidence_residues, find_high_confidence_residues


@pytest.fixture
def sample_pdb_file() -> Path:
    return Path(__file__).parent / "AF-A1YPR0-F1-model_v4.pdb"


@pytest.fixture
def sample_pdb(sample_pdb_file: Path) -> gemmi.Structure:
    return gemmi.read_structure(str(sample_pdb_file))


def test_find_high_confidence_residues(sample_pdb: gemmi.Structure):
    residues = list(find_high_confidence_residues(sample_pdb, 90))

    assert len(residues) == 22


def test_filter_out_low_confidence_residues(sample_pdb: gemmi.Structure):
    # Make sure we start with >22 residues
    assert len(sample_pdb[0][0]) == 619

    residues = set(find_high_confidence_residues(sample_pdb, 90))
    new_structure = filter_out_low_confidence_residues(sample_pdb, residues)

    assert len(new_structure[0][0]) == 22
