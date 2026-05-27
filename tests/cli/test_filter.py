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


class TestResolution:
    def test_default_grouping(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        """Test filter resolution ranks structures per UniProt accession."""
        fixtures_dir = Path(__file__).resolve().parents[1] / "fixtures"
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in ("2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz", "1amb_updated.cif.gz"):
            (input_dir / fixture_name).symlink_to(fixtures_dir / fixture_name)

        output_dir = tmp_path / "output"

        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--top",
            "1",
            "--copy-method",
            "symlink",
        ]

        main(argv)

        output_files = sorted(path.name for path in output_dir.iterdir())
        assert output_files == ["2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz"]

        captured = capsys.readouterr()
        assert "Filtering 3 files" in captured.err
        assert "Wrote 2 files to" in captured.err

    def test_uses_default_top(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        """Test filter resolution defaults to top 1000 and stays quiet about stats when unused."""
        fixtures_dir = Path(__file__).resolve().parents[1] / "fixtures"
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in ("2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz", "1amb_updated.cif.gz"):
            (input_dir / fixture_name).symlink_to(fixtures_dir / fixture_name)

        output_dir = tmp_path / "output"

        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--copy-method",
            "symlink",
        ]

        main(argv)

        output_files = sorted(path.name for path in output_dir.iterdir())
        assert output_files == ["1amb_updated.cif.gz", "2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz"]

        captured = capsys.readouterr()
        assert "Wrote 3 files to" in captured.err
        assert "Statistics written to" not in captured.err

    def test_write_stats(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        """Test filter resolution writes CSV stats even when output file uses .log extension."""
        fixtures_dir = Path(__file__).resolve().parents[1] / "fixtures"
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in ("2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz", "1amb_updated.cif.gz"):
            (input_dir / fixture_name).symlink_to(fixtures_dir / fixture_name)

        output_dir = tmp_path / "output"
        stats_fn = tmp_path / "some.log"

        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--copy-method",
            "symlink",
            "--write-stats",
            str(stats_fn),
        ]

        main(argv)

        with stats_fn.open() as handle:
            rows = list(csv.DictReader(handle))

        assert len(rows) == 3
        assert rows[0] == {
            "input_file": str(input_dir / "2Y29.cif.gz"),
            "id": "2Y29",
            "uniprot_accession": "P05067",
            "resolution": "2.3",
            "total_residue_count": "8",
            "is_alphafold": "False",
            "uniprot_start": "687",
            "uniprot_end": "692",
            "sequence_identity": "1.333",
            "chain_length": "8",
            "passed": "True",
            "output_file": str(output_dir / "2Y29.cif.gz"),
        }
        assert rows[1] == {
            "input_file": str(input_dir / "1amb_updated.cif.gz"),
            "id": "1AMB",
            "uniprot_accession": "P05067",
            "resolution": "0.0",
            "total_residue_count": "28",
            "is_alphafold": "False",
            "uniprot_start": "672",
            "uniprot_end": "699",
            "sequence_identity": "1.000",
            "chain_length": "28",
            "passed": "True",
            "output_file": str(output_dir / "1amb_updated.cif.gz"),
        }
        assert rows[2] == {
            "input_file": str(input_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
            "id": "AF-A0A0C5B5G6-F1",
            "uniprot_accession": "A0A0C5B5G6",
            "resolution": "0.0",
            "total_residue_count": "16",
            "is_alphafold": "True",
            "uniprot_start": "1",
            "uniprot_end": "16",
            "sequence_identity": "1.000",
            "chain_length": "16",
            "passed": "True",
            "output_file": str(output_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
        }

        captured = capsys.readouterr()
        assert "Statistics written to" in captured.err
        assert str(stats_fn) in captured.err

    def test_no_groupby_and_write_stats(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        fixtures_dir = Path(__file__).resolve().parents[1] / "fixtures"
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in ("2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz", "1amb_updated.cif.gz"):
            (input_dir / fixture_name).symlink_to(fixtures_dir / fixture_name)

        output_dir = tmp_path / "output"
        stats_fn = tmp_path / "stats.csv"

        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--no-group-by",
            "--top",
            "2",
            "--copy-method",
            "symlink",
            "--write-stats",
            str(stats_fn),
        ]

        main(argv)

        output_files = sorted(path.name for path in output_dir.iterdir())
        assert output_files == ["2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz"]

        with stats_fn.open() as handle:
            rows = list(csv.DictReader(handle))

        expected_rows = [
            {
                "input_file": str(input_dir / "2Y29.cif.gz"),
                "id": "2Y29",
                "uniprot_accession": "P05067",
                "resolution": "2.3",
                "total_residue_count": "8",
                "is_alphafold": "False",
                "uniprot_start": "687",
                "uniprot_end": "692",
                "sequence_identity": "1.333",
                "chain_length": "8",
                "passed": "True",
                "output_file": str(output_dir / "2Y29.cif.gz"),
            },
            {
                "input_file": str(input_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
                "id": "AF-A0A0C5B5G6-F1",
                "uniprot_accession": "A0A0C5B5G6",
                "resolution": "0.0",
                "total_residue_count": "16",
                "is_alphafold": "True",
                "uniprot_start": "1",
                "uniprot_end": "16",
                "sequence_identity": "1.000",
                "chain_length": "16",
                "passed": "True",
                "output_file": str(output_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
            },
            {
                "input_file": str(input_dir / "1amb_updated.cif.gz"),
                "id": "1AMB",
                "uniprot_accession": "P05067",
                "resolution": "0.0",
                "total_residue_count": "28",
                "is_alphafold": "False",
                "uniprot_start": "672",
                "uniprot_end": "699",
                "sequence_identity": "1.000",
                "chain_length": "28",
                "passed": "False",
                "output_file": "",
            },
        ]
        assert rows == expected_rows

        captured = capsys.readouterr()
        assert "global resolution ranking (no grouping)" in captured.err

    def test_no_group(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        fixtures_dir = Path(__file__).resolve().parents[1] / "fixtures"
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in ("2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz", "1amb_updated.cif.gz"):
            (input_dir / fixture_name).symlink_to(fixtures_dir / fixture_name)

        output_dir = tmp_path / "output"

        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--no-group-by",
            "--top",
            "2",
            "--copy-method",
            "symlink",
        ]

        main(argv)

        output_files = sorted(path.name for path in output_dir.iterdir())
        assert output_files == ["2Y29.cif.gz", "AF-A0A0C5B5G6-F1-model_v6.cif.gz"]

        captured = capsys.readouterr()
        assert "global resolution ranking (no grouping)" in captured.err

    def test_rejects_zero_top(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        argv = [
            "filter",
            "resolution",
            str(tmp_path),
            str(tmp_path),
            "--top",
            "0",
        ]

        with pytest.raises(SystemExit) as exc_info:
            main(argv)

        assert exc_info.value.code == 1
        captured = capsys.readouterr()
        assert "Invalid value" in captured.err
