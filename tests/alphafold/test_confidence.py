from pathlib import Path

import gemmi
import pytest

from protein_quest.alphafold.confidence import (
    ConfidenceFilterQuery,
    ConfidenceFilterResult,
    converter,
    filter_files_on_confidence,
    filter_out_low_confidence_residues,
    find_high_confidence_residues,
)
from protein_quest.structure.chains import nr_residues_in_chain
from protein_quest.structure.formats import read_structure


@pytest.fixture
def sample_pdb(af_pdb: Path) -> gemmi.Structure:
    return read_structure(af_pdb)


def test_find_high_confidence_residues(sample_pdb: gemmi.Structure):
    residues = list(find_high_confidence_residues(sample_pdb, 90))

    assert len(residues) == 22


def test_filter_out_low_confidence_residues(sample_pdb: gemmi.Structure):
    # Make sure we start with >22 residues
    assert len(sample_pdb[0][0]) == 619

    residues = set(find_high_confidence_residues(sample_pdb, 90))
    new_structure = filter_out_low_confidence_residues(sample_pdb, residues)

    assert len(new_structure[0][0]) == 22


def test_filter_files_on_confidence(af_pdb: Path, tmp_path: Path):
    input_files = [af_pdb]
    query = ConfidenceFilterQuery(
        confidence=90,
        max_residues=40,
        min_residues=10,
    )

    results = filter_files_on_confidence(input_files, query, tmp_path, scheduler_address="sequential")

    expected = [
        ConfidenceFilterResult(
            input_file=af_pdb.name,
            count=22,
            filtered_file=tmp_path / af_pdb.name,
        )
    ]

    assert results == expected
    assert results[0].filtered_file is not None
    assert results[0].filtered_file.exists()
    assert nr_residues_in_chain(results[0].filtered_file) == 22


def test_query_converter():
    result = converter.structure(
        {
            "confidence": 90.0,
            "min_residues": 10,
            "max_residues": 100,
        },
        ConfidenceFilterQuery,
    )
    expected = ConfidenceFilterQuery(
        confidence=90.0,
        min_residues=10,
        max_residues=100,
    )
    assert result == expected


@pytest.mark.parametrize(
    "raw,match",
    [
        (
            {
                "confidence": 42,
                "min_residues": -10,
                "max_residues": 100,
            },
            "is not a valid positive integer",
        ),
        (
            {
                "confidence": 1234,
                "min_residues": 10,
                "max_residues": 100,
            },
            "is not a valid percentage",
        ),
        (
            {
                "confidence": -10,
                "min_residues": 10,
                "max_residues": 100,
            },
            "is not a valid percentage",
        ),
        (
            {
                "confidence": "not a number",
                "min_residues": 10,
                "max_residues": 100,
            },
            "could not convert string to float",
        ),
        (
            {
                "confidence": 42,
                "min_residues": "not a number",
                "max_residues": 100,
            },
            "invalid literal for int",
        ),
    ],
)
def test_query_converter_bad_confidence(raw: dict[str, str | int], match: str):
    with pytest.RaisesGroup(pytest.RaisesExc(ValueError, match=match)):
        converter.structure(
            raw,
            ConfidenceFilterQuery,
        )


@pytest.mark.parametrize(
    "raw,match",
    [
        (
            {
                "confidence": 1,
                "min_residues": 80,
                "max_residues": 20,
            },
            "cannot be larger than",
        ),
    ],
)
def test_query_converter_bad_residues(raw: dict[str, str | int], match: str):
    with pytest.raises(ValueError, match=match):
        converter.structure(
            raw,
            ConfidenceFilterQuery,
        )
