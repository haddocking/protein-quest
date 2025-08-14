from pathlib import Path

import pytest

from protein_quest.alphafold.confidence import filter_out_low_confidence_residues, find_high_confidence_residues


@pytest.fixture
def sample_pdb() -> Path:
    return Path(__file__).parent / "AF-A1YPR0-F1-model_v4.pdb"


def test_find_high_confidence_residues(sample_pdb: Path):
    residues = list(find_high_confidence_residues(sample_pdb, 90))

    assert len(residues) == 22


def test_filter_out_low_confidence_residues(sample_pdb: Path, tmp_path: Path):
    residues = set(find_high_confidence_residues(sample_pdb, 90))
    out_pdb_file = tmp_path / "filtered.pdb"

    filter_out_low_confidence_residues(sample_pdb, residues, out_pdb_file)

    assert out_pdb_file.stat().st_size < sample_pdb.stat().st_size
