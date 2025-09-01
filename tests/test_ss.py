from pathlib import Path

import gemmi
import pytest

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
def sample_cif() -> Path:
    """Downloaded from https://www.rcsb.org/structure/3JRS
    and filtered with
    `write_single_chain_pdb_file(Path('tests/fixtures/3JRS.cif.gz'), 'A', Path('tests/fixtures/'))`
    """
    return Path(__file__).parent / "fixtures" / "3JRS_A2A.cif.gz"


@pytest.fixture
def sample_structure(sample_cif: Path) -> gemmi.Structure:
    return gemmi.read_structure(str(sample_cif))


@pytest.fixture(scope="module")
def sample_stats() -> SecondaryStructureStats:
    return SecondaryStructureStats(
        nr_residues=174,
        nr_helix_residues=58,
        nr_sheet_residues=59,
        helix_ratio=58 / 174,
        sheet_ratio=59 / 174,
    )


def test_nr_of_residues_in_helix(sample_structure: gemmi.Structure):
    expected_count = 58
    assert nr_of_residues_in_helix(sample_structure) == expected_count


def test_nr_of_residues_in_sheet(sample_structure: gemmi.Structure):
    expected_count = 59
    assert nr_of_residues_in_sheet(sample_structure) == expected_count


def test_nr_of_residues_in_total(sample_structure: gemmi.Structure):
    expected_count = 174
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
        (SecondaryStructureFilterQuery(ratio_min_helix_residues=58 / 174), True),
        (SecondaryStructureFilterQuery(ratio_min_helix_residues=0.4), False),
        # ratio_max_helix_residues
        (SecondaryStructureFilterQuery(ratio_max_helix_residues=58 / 174), True),
        (SecondaryStructureFilterQuery(ratio_max_helix_residues=0.3), False),
        # ratio_min_sheet_residues
        (SecondaryStructureFilterQuery(ratio_min_sheet_residues=59 / 174), True),
        (SecondaryStructureFilterQuery(ratio_min_sheet_residues=0.4), False),
        # ratio_max_sheet_residues
        (SecondaryStructureFilterQuery(ratio_max_sheet_residues=59 / 174), True),
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
    result = filter_on_secondary_structure(sample_structure, query)
    expected = SecondaryStructureFilterResult(
        stats=sample_stats,
        passed=expected_passed,
    )
    assert result == expected


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
