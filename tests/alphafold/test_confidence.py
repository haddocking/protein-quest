from pathlib import Path

import gemmi
import pytest

from protein_quest.alphafold.confidence import (
    ConfidenceFilterQuery,
    ConfidenceFilterResult,
    filter_files_on_confidence,
    filter_out_low_confidence_residues,
    find_high_confidence_residues,
)
from protein_quest.pdbe.io import nr_residues_in_chain


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


def test_filter_files_on_confidence(sample_pdb_file: Path, tmp_path: Path):
    input_files = [sample_pdb_file]
    query = ConfidenceFilterQuery(
        confidence=90,
        max_threshold=40,
        min_threshold=10,
    )

    results = list(filter_files_on_confidence(input_files, query, tmp_path))

    expected = [
        ConfidenceFilterResult(
            input_file=sample_pdb_file.name,
            count=22,
            filtered_file=tmp_path / sample_pdb_file.name,
        )
    ]

    assert results == expected
    assert results[0].filtered_file is not None
    assert results[0].filtered_file.exists()
    assert nr_residues_in_chain(results[0].filtered_file) == 22
