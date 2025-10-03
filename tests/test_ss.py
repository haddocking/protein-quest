from pathlib import Path

import gemmi
import pytest

from protein_quest.converter import converter
from protein_quest.io import read_structure
from protein_quest.ss import (
    SecondaryStructureFilterQuery,
    SecondaryStructureFilterResult,
    SecondaryStructureStats,
    filter_file_on_secondary_structure,
    filter_files_on_secondary_structure,
    filter_on_secondary_structure,
    nr_of_residues_in_helix,
    nr_of_residues_in_sheet,
    nr_of_residues_in_total,
)


@pytest.fixture
def sample_structure(sample_cif: Path) -> gemmi.Structure:
    return read_structure(sample_cif)


@pytest.fixture(scope="module")
def sample_stats() -> SecondaryStructureStats:
    return SecondaryStructureStats(
        nr_residues=173,
        nr_helix_residues=58,
        nr_sheet_residues=59,
        helix_ratio=58 / 173,
        sheet_ratio=59 / 173,
    )


def test_nr_of_residues_in_helix(sample_structure: gemmi.Structure):
    expected_count = 58
    assert nr_of_residues_in_helix(sample_structure) == expected_count


def test_nr_of_residues_in_sheet(sample_structure: gemmi.Structure):
    expected_count = 59
    assert nr_of_residues_in_sheet(sample_structure) == expected_count


def test_nr_of_residues_in_total(sample_structure: gemmi.Structure):
    expected_count = 173
    assert nr_of_residues_in_total(sample_structure) == expected_count


@pytest.mark.parametrize(
    "query, expected_passed",
    [
        # abs_min_helix_residues
        (SecondaryStructureFilterQuery(abs_min_helix_residues=58), True),
        (SecondaryStructureFilterQuery(abs_min_helix_residues=59), False),
        # abs_max_helix_residues
        (SecondaryStructureFilterQuery(abs_max_helix_residues=58), True),
        (SecondaryStructureFilterQuery(abs_max_helix_residues=57), False),
        # abs_min_sheet_residues
        (SecondaryStructureFilterQuery(abs_min_sheet_residues=59), True),
        (SecondaryStructureFilterQuery(abs_min_sheet_residues=60), False),
        # abs_max_sheet_residues
        (SecondaryStructureFilterQuery(abs_max_sheet_residues=59), True),
        (SecondaryStructureFilterQuery(abs_max_sheet_residues=58), False),
        # ratio_min_helix_residues
        (SecondaryStructureFilterQuery(ratio_min_helix_residues=58 / 173), True),
        (SecondaryStructureFilterQuery(ratio_min_helix_residues=0.4), False),
        # ratio_max_helix_residues
        (SecondaryStructureFilterQuery(ratio_max_helix_residues=58 / 173), True),
        (SecondaryStructureFilterQuery(ratio_max_helix_residues=0.3), False),
        # ratio_min_sheet_residues
        (SecondaryStructureFilterQuery(ratio_min_sheet_residues=59 / 173), True),
        (SecondaryStructureFilterQuery(ratio_min_sheet_residues=0.4), False),
        # ratio_max_sheet_residues
        (SecondaryStructureFilterQuery(ratio_max_sheet_residues=59 / 173), True),
        (SecondaryStructureFilterQuery(ratio_max_sheet_residues=0.3), False),
        # multiple
        (
            SecondaryStructureFilterQuery(
                ratio_min_helix_residues=0.1,
                ratio_min_sheet_residues=0.1,
            ),
            True,
        ),
        (
            SecondaryStructureFilterQuery(
                ratio_min_helix_residues=0.9,
                ratio_min_sheet_residues=0.9,
            ),
            False,
        ),
    ],
)
def test_filter_on_secondary_structure(
    sample_structure: gemmi.Structure,
    sample_stats: SecondaryStructureStats,
    query: SecondaryStructureFilterQuery,
    expected_passed: bool,
):
    assert query.is_actionable()
    result = filter_on_secondary_structure(sample_structure, query)
    expected = SecondaryStructureFilterResult(
        stats=sample_stats,
        passed=expected_passed,
    )
    assert result == expected


def test_SecondaryStructureFilterQuery_non_actionable():
    query = SecondaryStructureFilterQuery()
    assert not query.is_actionable()


def test_filter_on_secondary_structure_raises_on_zero_conditions(sample_structure: gemmi.Structure):
    query = SecondaryStructureFilterQuery()
    with pytest.raises(ValueError, match="No filtering conditions provided"):
        filter_on_secondary_structure(sample_structure, query)


def test_filter_on_secondary_structure_raise_on_zero_residues():
    structure = gemmi.Structure()
    query = SecondaryStructureFilterQuery(abs_min_helix_residues=1)
    with pytest.raises(ValueError, match="Structure has zero residues"):
        filter_on_secondary_structure(structure, query)


def test_filter_file_on_secondary_structure(sample_cif: Path, sample_stats: SecondaryStructureStats):
    query = SecondaryStructureFilterQuery(abs_min_helix_residues=1)
    result = filter_file_on_secondary_structure(sample_cif, query)
    expected = SecondaryStructureFilterResult(stats=sample_stats, passed=True)
    assert result == expected


def test_filter_files_on_secondary_structure(sample_cif: Path, sample_stats: SecondaryStructureStats):
    query = SecondaryStructureFilterQuery(abs_min_helix_residues=1)
    result = filter_files_on_secondary_structure([sample_cif], query)
    expected = {sample_cif: SecondaryStructureFilterResult(stats=sample_stats, passed=True)}
    assert dict(result) == expected


def test_converter():
    raw = {
        "abs_min_helix_residues": 10,
        "abs_max_helix_residues": 20,
        "abs_min_sheet_residues": 5,
        "abs_max_sheet_residues": 15,
        "ratio_min_helix_residues": 0.1,
        "ratio_max_helix_residues": 0.3,
        "ratio_min_sheet_residues": 0.2,
        "ratio_max_sheet_residues": 0.4,
    }
    result = converter.structure(raw, SecondaryStructureFilterQuery)
    expected = SecondaryStructureFilterQuery(
        abs_min_helix_residues=10,
        abs_max_helix_residues=20,
        abs_min_sheet_residues=5,
        abs_max_sheet_residues=15,
        ratio_min_helix_residues=0.1,
        ratio_max_helix_residues=0.3,
        ratio_min_sheet_residues=0.2,
        ratio_max_sheet_residues=0.4,
    )
    assert result == expected


@pytest.mark.parametrize(
    "raw, match",
    [
        (
            {
                "abs_min_helix_residues": 20,
                "abs_max_helix_residues": 10,
            },
            "must be smaller than max",
        ),
        (
            {
                "abs_min_sheet_residues": 15,
                "abs_max_sheet_residues": 5,
            },
            "must be smaller than max",
        ),
        (
            {
                "ratio_min_helix_residues": 0.4,
                "ratio_max_helix_residues": 0.3,
            },
            "must be smaller than max",
        ),
        (
            {
                "ratio_min_sheet_residues": 0.5,
                "ratio_max_sheet_residues": 0.4,
            },
            "must be smaller than max",
        ),
    ],
)
def test_converter_raises_on_invalid_range(raw, match):
    with pytest.raises(ValueError, match=match):
        converter.structure(raw, SecondaryStructureFilterQuery)


@pytest.mark.parametrize(
    "raw, match",
    [
        (
            {
                "ratio_min_sheet_residues": 4.2,
            },
            "is not a valid ratio",
        ),
        (
            {
                "ratio_min_sheet_residues": "not a number",
            },
            "could not convert string to float",
        ),
    ],
)
def test_converter_raises_on_invalid_ratio(raw, match):
    with pytest.RaisesGroup(pytest.RaisesExc(ValueError, match=match)):
        converter.structure(raw, SecondaryStructureFilterQuery)
