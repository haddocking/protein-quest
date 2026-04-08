"""Filter CLI tests for protein-quest."""

import csv
from pathlib import Path

import pytest

from protein_quest.cli import main


def test_filter_chain_happy_path(sample2_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test filter chain command with valid input."""
    chains_fn = tmp_path / "chains.csv"
    chains_fn.write_text("pdb_id,chain\n2Y29,A\n")

    argv = [
        "filter",
        "chain",
        str(chains_fn),
        str(sample2_cif.parent),
        str(tmp_path),
        "--copy-method",
        "copy",
    ]

    main(argv)

    output_file = tmp_path / "2Y29_A2A.cif.gz"
    assert output_file.exists()

    captured = capsys.readouterr()
    assert "Wrote 1 single-chain PDB/mmCIF files to" in captured.err


def test_filter_chain_input_file_notfound(tmp_path: Path):
    """Test filter chain command with missing input file."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    chains_fn = tmp_path / "chains.csv"
    chains_fn.write_text("pdb_id,chain\n2Y29,A\n")

    argv = [
        "filter",
        "chain",
        str(chains_fn),
        str(input_dir),
        str(output_dir),
    ]

    with pytest.raises(ValueError, match="No valid structure files found"):
        main(argv)


def test_filter_residue(sample_cif: Path, sample2_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test filter residue command."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    local_sample = input_dir / sample_cif.name
    local_sample.symlink_to(sample_cif)
    local_sample2 = input_dir / sample2_cif.name
    local_sample2.symlink_to(sample2_cif)
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    stats_fn = tmp_path / "stats.csv"

    argv = [
        "filter",
        "residue",
        str(input_dir),
        str(output_dir),
        "--min-residues",
        "100",
        "--max-residues",
        "200",
        "--copy-method",
        "symlink",
        "--write-stats",
        str(stats_fn),
    ]

    main(argv)

    # Check output files
    output_files = list(output_dir.iterdir())
    assert len(output_files) == 1
    expected_passed_file = output_dir / sample_cif.name
    assert expected_passed_file in output_files

    # Check stats file
    with stats_fn.open() as f:
        rows = list(csv.DictReader(f))
    # Input files processed in alphabetical order
    expected_stats = [
        {
            "input_file": str(local_sample2),
            "residue_count": "8",
            "passed": "False",
            "output_file": "",
        },
        {
            "input_file": str(local_sample),
            "residue_count": "173",
            "passed": "True",
            "output_file": str(expected_passed_file),
        },
    ]
    assert rows == expected_stats

    # Check captured output
    captured = capsys.readouterr()
    assert "by number of residues in chain A" in captured.err
    assert "Wrote 1 files to" in captured.err
    assert "Statistics written to" in captured.err


def test_filter_secondary_structure(
    sample_cif: Path, sample2_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    """Test filter secondary-structure command."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    local_sample = input_dir / sample_cif.name
    local_sample.symlink_to(sample_cif)
    local_sample2 = input_dir / sample2_cif.name
    local_sample2.symlink_to(sample2_cif)
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    stats_fn = tmp_path / "ss_stats.csv"

    argv = [
        "filter",
        "secondary-structure",
        str(input_dir),
        str(output_dir),
        "--abs-min-helix-residues",
        "10",
        "--copy-method",
        "symlink",
        "--write-stats",
        str(stats_fn),
    ]

    main(argv)

    # Check output files
    output_files = list(output_dir.iterdir())
    assert len(output_files) == 1
    expected_passed_file = output_dir / sample_cif.name
    assert expected_passed_file in output_files

    # Check stats file
    with stats_fn.open() as f:
        rows = list(csv.DictReader(f))
    expected_stats = [
        {
            "helix_ratio": "0.0",
            "input_file": str(local_sample2),
            "nr_helix_residues": "0",
            "nr_residues": "8",
            "nr_sheet_residues": "0",
            "output_file": "",
            "passed": "False",
            "sheet_ratio": "0.0",
        },
        {
            "input_file": str(local_sample),
            "nr_residues": "173",
            "nr_helix_residues": "58",
            "nr_sheet_residues": "59",
            "helix_ratio": f"{58 / 173:.3f}",
            "sheet_ratio": f"{59 / 173:.3f}",
            "passed": "True",
            "output_file": str(expected_passed_file),
        },
    ]
    assert rows == expected_stats

    # Check captured output
    captured = capsys.readouterr()
    assert "by secondary structure" in captured.err
    assert "Wrote 1 files to" in captured.err
    assert "Statistics written to" in captured.err
