from pathlib import Path

from protein_quest.filters.residues import (
    ResidueFilterStatistics,
    filter_files_on_residues,
)


def test_filter_files_on_residues_min_max_range(sample_cif: Path, sample2_cif: Path, tmp_path: Path):
    results = list(
        filter_files_on_residues(
            input_files=[sample_cif, sample2_cif],
            output_dir=tmp_path,
            min_residues=100,
            max_residues=200,
        )
    )
    expected_passed = ResidueFilterStatistics(
        input_file=sample_cif,
        residue_count=173,
        passed=True,
        output_file=tmp_path / sample_cif.name,
    )
    assert expected_passed.output_file and expected_passed.output_file.exists()
    expected_discarded = ResidueFilterStatistics(
        input_file=sample2_cif,
        residue_count=8,
        passed=False,
        output_file=None,
    )

    assert results == [expected_passed, expected_discarded]
