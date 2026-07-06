"""Filter CLI tests for protein-quest."""

import csv
import json
import textwrap
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
        "--scheduler-address",
        "sequential",
    ]

    main(argv)

    output_file = tmp_path / "2Y29_A2A.cif.gz"
    assert output_file.exists()

    captured = capsys.readouterr()
    assert "Wrote 1 single-chain PDB/mmCIF files to" in captured.err


def test_filter_chain_multi_accession_happy_path(
    multi_accession_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    """Test filter chain command handles multi-accession structures."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    local_multi = input_dir / multi_accession_cif.name
    local_multi.symlink_to(multi_accession_cif)

    chains_fn = tmp_path / "chains.csv"
    # 1A02 maps the three UniProt accessions to chains N, F, and J.
    chains_fn.write_text(
        textwrap.dedent(
            """\
            pdb_id,chain
            1A02,N
            1A02,F
            1A02,J
            """
        )
    )

    argv = [
        "filter",
        "chain",
        str(chains_fn),
        str(input_dir),
        str(tmp_path),
        "--scheduler-address",
        "sequential",
        "--copy-method",
        "copy",
    ]

    main(argv)

    output_files = {path.name for path in tmp_path.glob("1a02_*2A.cif.gz")}
    assert output_files == {"1a02_N2A.cif.gz", "1a02_F2A.cif.gz", "1a02_J2A.cif.gz"}

    captured = capsys.readouterr()
    assert "Wrote 3 single-chain PDB/mmCIF files to" in captured.err


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


@pytest.mark.parametrize(
    ("chain_id", "extra_args"),
    [
        pytest.param("B", [], id="defaults-to-auth"),
        pytest.param("A", ["--chain-system", "label"], id="explicit-label-system"),
    ],
)
def test_filter_chain_system_modes(cif_8rw8: Path, tmp_path: Path, chain_id: str, extra_args: list[str]):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    local_cif = input_dir / cif_8rw8.name
    local_cif.symlink_to(cif_8rw8)

    chains_fn = tmp_path / "chains.csv"
    chains_fn.write_text(f"pdb_id,chain\n8rw8,{chain_id}\n")

    argv = [
        "filter",
        "chain",
        str(chains_fn),
        str(input_dir),
        str(tmp_path),
        *extra_args,
        "--scheduler-address",
        "sequential",
    ]

    main(argv)

    output_files = {path.name for path in tmp_path.glob("*8rw8*_B2A.cif.gz")}
    assert len(output_files) == 1


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
    assert "by number of residues in chain A" in captured.err.replace("\n", " ")
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
            "nr_helix_residues": "54",
            "nr_sheet_residues": "61",
            "helix_ratio": f"{54 / 173:.3f}",
            "sheet_ratio": f"{61 / 173:.3f}",
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
    def test_with_defaults_and_stats(
        self, sample_cif: Path, sample2_cif: Path, af_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in (sample_cif, sample2_cif, af_cif):
            (input_dir / fixture_name.name).symlink_to(fixture_name)

        output_dir = tmp_path / "output"
        stats_fn = tmp_path / "stats.csv"
        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--copy-method",
            "symlink",
            "--scheduler-address",
            "sequential",
            "--write-stats",
            str(stats_fn),
        ]

        main(argv)

        output_files = {path.name for path in output_dir.iterdir()}
        assert output_files == {"2Y29.cif.gz", "3jrs_updated_B2A.cif.gz"}

        captured = capsys.readouterr()
        assert "Filtering 3 files" in captured.err
        assert "Wrote 2 files to" in captured.err

        stats = list(csv.DictReader(stats_fn.open()))
        expected_stats = [
            {
                "input_file": (input_dir / sample_cif.name).as_posix(),
                "id": "3JRS",
                "uniprot_accession": "Q8VZS8",
                "resolution": "2.05",
                "total_residue_count": "173",
                "is_alphafold": "False",
                "uniprot_start": "8",
                "uniprot_end": "211",
                "sequence_identity": "1.000",
                "chain_length": "173",
                "passed": "True",
                "output_file": (output_dir / sample_cif.name).as_posix(),
                "discard_reason": "",
                "discard_reason_type": "",
            },
            {
                "input_file": (input_dir / sample2_cif.name).as_posix(),
                "id": "2Y29",
                "uniprot_accession": "P05067",
                "resolution": "2.3",
                "total_residue_count": "8",
                "is_alphafold": "False",
                "uniprot_start": "687",
                "uniprot_end": "692",
                "sequence_identity": "1.000",
                "chain_length": "8",
                "passed": "True",
                "output_file": (output_dir / sample2_cif.name).as_posix(),
                "discard_reason": "",
                "discard_reason_type": "",
            },
            {
                "input_file": (input_dir / af_cif.name).as_posix(),
                "id": "AF-A0A0C5B5G6-F1",
                "uniprot_accession": "A0A0C5B5G6",
                "resolution": "0.0",
                "total_residue_count": "16",
                "is_alphafold": "True",
                "uniprot_start": "1",
                "uniprot_end": "16",
                "sequence_identity": "1.000",
                "chain_length": "16",
                "passed": "False",
                "output_file": "",
                "discard_reason": f"Resolution is unset for {(input_dir / af_cif.name).as_posix()}",
                "discard_reason_type": "ResolutionUnsetError",
            },
        ]
        assert stats == expected_stats

    def test_with_all_options(self, all_cifs: list[Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in all_cifs:
            (input_dir / fixture_name.name).symlink_to(fixture_name)

        output_dir = tmp_path / "output"
        stats_fn = tmp_path / "stats.csv"
        argv = [
            "filter",
            "resolution",
            str(input_dir),
            str(output_dir),
            "--copy-method",
            "symlink",
            "--scheduler-address",
            "sequential",
            "--write-stats",
            str(stats_fn),
            "--top",
            "2",
            "--no-group-by-uniprot-accession",
            "--no-coverage",
            "--min-sequence-identity",
            "0.9",
            "--lax",
        ]

        main(argv)

        output_files = {path.name for path in output_dir.iterdir()}
        assert output_files == {
            "1amb_updated.cif.gz",
            "2Y29.cif.gz",
            "3jrs_updated_B2A.cif.gz",
            "AF-A0A0C5B5G6-F1-model_v6.cif.gz",
        }

        captured = capsys.readouterr()
        assert "Filtering 7 files" in captured.err
        assert "Wrote 4 files to" in captured.err
        assert "Additionally wrote 2 files to" in captured.err

        stats = list(csv.DictReader(stats_fn.open()))
        expected_stats = [
            {
                "chain_length": "173",
                "discard_reason": "",
                "discard_reason_type": "",
                "id": "3JRS",
                "input_file": str(input_dir / "3jrs_updated_B2A.cif.gz"),
                "is_alphafold": "False",
                "output_file": str(output_dir / "3jrs_updated_B2A.cif.gz"),
                "passed": "True",
                "resolution": "2.05",
                "sequence_identity": "1.000",
                "total_residue_count": "173",
                "uniprot_accession": "Q8VZS8",
                "uniprot_end": "211",
                "uniprot_start": "8",
            },
            {
                "chain_length": "8",
                "discard_reason": "",
                "discard_reason_type": "",
                "id": "2Y29",
                "input_file": str(input_dir / "2Y29.cif.gz"),
                "is_alphafold": "False",
                "output_file": str(output_dir / "2Y29.cif.gz"),
                "passed": "True",
                "resolution": "2.3",
                "sequence_identity": "1.000",
                "total_residue_count": "8",
                "uniprot_accession": "P05067",
                "uniprot_end": "692",
                "uniprot_start": "687",
            },
            {
                "chain_length": "260",
                "discard_reason": "Rank 3 > top 2",
                "discard_reason_type": "OutsideTopError",
                "id": "8W77",
                "input_file": str(input_dir / "8w77_updated.cif.gz"),
                "is_alphafold": "False",
                "output_file": "",
                "passed": "False",
                "resolution": "3.61",
                "sequence_identity": "1.000",
                "total_residue_count": "260",
                "uniprot_accession": "P0ABE7",
                "uniprot_end": "127",
                "uniprot_start": "23",
            },
            {
                "chain_length": "131",
                "discard_reason": "Rank 4 > top 2",
                "discard_reason_type": "OutsideTopError",
                "id": "1UN5",
                "input_file": str(input_dir / "1un5.cif.gz"),
                "is_alphafold": "False",
                "output_file": "",
                "passed": "False",
                "resolution": "2.6",
                "sequence_identity": "0.967",
                "total_residue_count": "131",
                "uniprot_accession": "P03950",
                "uniprot_end": "147",
                "uniprot_start": "25",
            },
            {
                "chain_length": "1346",
                "discard_reason": f"Sequence identity 0.816 below minimal 0.900 for {input_dir / '6O5I.cif.gz'}",
                "discard_reason_type": "SequenceIdentityBelowThresholdError",
                "id": "6O5I",
                "input_file": str(input_dir / "6O5I.cif.gz"),
                "is_alphafold": "False",
                "output_file": "",
                "passed": "False",
                "resolution": "1.24",
                "sequence_identity": "0.816",
                "total_residue_count": "1346",
                "uniprot_accession": "O00255",
                "uniprot_end": "593",
                "uniprot_start": "1",
            },
            {
                "chain_length": "28",
                "discard_reason": f"Resolution is unset for {input_dir / '1amb_updated.cif.gz'}",
                "discard_reason_type": "ResolutionUnsetError",
                "id": "1AMB",
                "input_file": str(input_dir / "1amb_updated.cif.gz"),
                "is_alphafold": "False",
                "output_file": str(output_dir / "1amb_updated.cif.gz"),
                "passed": "True",
                "resolution": "0.0",
                "sequence_identity": "1.000",
                "total_residue_count": "28",
                "uniprot_accession": "P05067",
                "uniprot_end": "699",
                "uniprot_start": "672",
            },
            {
                "chain_length": "16",
                "discard_reason": f"Resolution is unset for {input_dir / 'AF-A0A0C5B5G6-F1-model_v6.cif.gz'}",
                "discard_reason_type": "ResolutionUnsetError",
                "id": "AF-A0A0C5B5G6-F1",
                "input_file": str(input_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
                "is_alphafold": "True",
                "output_file": str(output_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
                "passed": "True",
                "resolution": "0.0",
                "sequence_identity": "1.000",
                "total_residue_count": "16",
                "uniprot_accession": "A0A0C5B5G6",
                "uniprot_end": "16",
                "uniprot_start": "1",
            },
        ]
        assert stats == expected_stats

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


def test_pdbe_quality(
    tmp_path: Path,
    all_cifs: list[Path],
):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for fixture_name in all_cifs:
        (input_dir / fixture_name.name).symlink_to(fixture_name)

    # Fetched using `protein-quest search pdbe-quality ./cif_fixture.ids ./cif_fixture.qual.json`
    quality_json = tmp_path / "quality.json"
    quality_json.write_text(
        json.dumps(
            {
                "1un5": {
                    "geometry_quality": 45.23,
                    "data_quality": None,
                    "overall_quality": 45.23,
                    "experiment_data_available": False,
                },
                "1a02": {
                    "geometry_quality": 6.15,
                    "data_quality": 7.84,
                    "overall_quality": 6.73,
                    "experiment_data_available": True,
                },
                "2y29": {
                    "geometry_quality": 55.9,
                    "data_quality": 81.09,
                    "overall_quality": 63.83,
                    "experiment_data_available": True,
                },
                "3jrs": {
                    "geometry_quality": 31.58,
                    "data_quality": 22.4,
                    "overall_quality": 27.13,
                    "experiment_data_available": True,
                },
                "1amb": {
                    "geometry_quality": None,
                    "data_quality": None,
                    "overall_quality": None,
                    "experiment_data_available": "unknown",
                },
                "6o5i": {
                    "geometry_quality": 64.38,
                    "data_quality": 61.45,
                    "overall_quality": 63.17,
                    "experiment_data_available": True,
                },
                "8w77": {
                    "geometry_quality": 75.0,
                    "data_quality": 75.0,
                    "overall_quality": 75.0,
                    "experiment_data_available": True,
                },
            }
        )
    )
    output_dir = tmp_path / "output"
    stats_csv = tmp_path / "stats.csv"

    argv = [
        "filter",
        "pdbe-quality",
        str(input_dir),
        str(quality_json),
        str(output_dir),
        "--write-stats",
        str(stats_csv),
    ]
    main(argv)

    assert output_dir.exists()
    assert stats_csv.exists()

    stats = list(csv.DictReader(stats_csv.open()))
    expected_stats = [
        {
            "geometry_quality": "75.0",
            "input_file": str(input_dir / "8w77_updated.cif.gz"),
            "output_file": str(output_dir / "8w77_updated.cif.gz"),
            "passed": "True",
            "pdb_id": "8w77",
            "reason": "",
        },
        {
            "pdb_id": "6o5i",
            "input_file": str(input_dir / "6O5I.cif.gz"),
            "geometry_quality": "64.38",
            "passed": "True",
            "output_file": str(output_dir / "6O5I.cif.gz"),
            "reason": "",
        },
        {
            "pdb_id": "2y29",
            "input_file": str(input_dir / "2Y29.cif.gz"),
            "geometry_quality": "55.9",
            "passed": "True",
            "output_file": str(output_dir / "2Y29.cif.gz"),
            "reason": "",
        },
        {
            "pdb_id": "1un5",
            "input_file": str(input_dir / "1un5.cif.gz"),
            "geometry_quality": "45.23",
            "passed": "False",
            "output_file": "",
            "reason": "Geometry quality score 45.23 < 50.0",
        },
        {
            "pdb_id": "3jrs",
            "input_file": str(input_dir / "3jrs_updated_B2A.cif.gz"),
            "geometry_quality": "31.58",
            "passed": "False",
            "output_file": "",
            "reason": "Geometry quality score 31.58 < 50.0",
        },
        {
            "pdb_id": "1amb",
            "input_file": str(input_dir / "1amb_updated.cif.gz"),
            "geometry_quality": "",
            "passed": "False",
            "output_file": "",
            "reason": "No geometry quality score",
        },
        {
            "pdb_id": "1a02",
            "input_file": "",
            "geometry_quality": "6.15",
            "passed": "False",
            "output_file": "",
            "reason": "File not found",
        },
        {
            "pdb_id": "",
            "input_file": str(input_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
            "geometry_quality": "",
            "passed": "False",
            "output_file": "",
            "reason": "File not found in quality scores",
        },
    ]
    assert stats == expected_stats


def test_pdbe_quality_with_clustering(
    tmp_path: Path,
    all_cifs: list[Path],
):
    """Test filter pdbe-quality command with clustering by UniProt accession and coverage."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for fixture_name in all_cifs:
        (input_dir / fixture_name.name).symlink_to(fixture_name)

    # Fetched using `protein-quest search pdbe-quality ./cif_fixture.ids ./cif_fixture.qual.json`
    quality_json = tmp_path / "quality.json"
    quality_json.write_text(
        json.dumps(
            {
                "1un5": {
                    "geometry_quality": 45.23,
                    "data_quality": None,
                    "overall_quality": 45.23,
                    "experiment_data_available": False,
                },
                "1a02": {
                    "geometry_quality": 6.15,
                    "data_quality": 7.84,
                    "overall_quality": 6.73,
                    "experiment_data_available": True,
                },
                "2y29": {
                    "geometry_quality": 55.9,
                    "data_quality": 81.09,
                    "overall_quality": 63.83,
                    "experiment_data_available": True,
                },
                "3jrs": {
                    "geometry_quality": 31.58,
                    "data_quality": 22.4,
                    "overall_quality": 27.13,
                    "experiment_data_available": True,
                },
                "1amb": {
                    "geometry_quality": None,
                    "data_quality": None,
                    "overall_quality": None,
                    "experiment_data_available": "unknown",
                },
                "6o5i": {
                    "geometry_quality": 64.38,
                    "data_quality": 61.45,
                    "overall_quality": 63.17,
                    "experiment_data_available": True,
                },
                "8w77": {
                    "geometry_quality": 75.0,
                    "data_quality": 75.0,
                    "overall_quality": 75.0,
                    "experiment_data_available": True,
                },
            }
        )
    )
    output_dir = tmp_path / "output"
    stats_csv = tmp_path / "stats.csv"

    argv = [
        "filter",
        "pdbe-quality",
        str(input_dir),
        str(quality_json),
        str(output_dir),
        "--write-stats",
        str(stats_csv),
        "--minimal-geometry-quality",
        "30.0",
        "--cluster-by-uniprot-accession-and-coverage",
        "1",
    ]
    main(argv)

    assert output_dir.exists()
    assert stats_csv.exists()

    stats = list(csv.DictReader(stats_csv.open()))
    expected_stats = [
        {
            "geometry_quality": "75.0",
            "input_file": str(tmp_path / "input" / "8w77_updated.cif.gz"),
            "output_file": str(tmp_path / "output" / "8w77_updated.cif.gz"),
            "passed": "True",
            "pdb_id": "8w77",
            "reason": "",
        },
        {
            "geometry_quality": "64.38",
            "input_file": str(tmp_path / "input" / "6O5I.cif.gz"),
            "output_file": str(tmp_path / "output" / "6O5I.cif.gz"),
            "passed": "True",
            "pdb_id": "6o5i",
            "reason": "",
        },
        {
            "geometry_quality": "55.9",
            "input_file": str(tmp_path / "input" / "2Y29.cif.gz"),
            "output_file": str(tmp_path / "output" / "2Y29.cif.gz"),
            "passed": "True",
            "pdb_id": "2y29",
            "reason": "",
        },
        {
            "geometry_quality": "45.23",
            "input_file": str(tmp_path / "input" / "1un5.cif.gz"),
            "output_file": str(tmp_path / "output" / "1un5.cif.gz"),
            "passed": "True",
            "pdb_id": "1un5",
            "reason": "",
        },
        {
            "geometry_quality": "31.58",
            "input_file": str(tmp_path / "input" / "3jrs_updated_B2A.cif.gz"),
            "output_file": str(tmp_path / "output" / "3jrs_updated_B2A.cif.gz"),
            "passed": "True",
            "pdb_id": "3jrs",
            "reason": "",
        },
        {
            "geometry_quality": "",
            "input_file": str(tmp_path / "input" / "1amb_updated.cif.gz"),
            "output_file": "",
            "passed": "False",
            "pdb_id": "1amb",
            "reason": "No geometry quality score",
        },
        {
            "geometry_quality": "",
            "input_file": str(tmp_path / "input" / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"),
            "output_file": "",
            "passed": "False",
            "pdb_id": "af-a0a0c5b5g6-f1",
            "reason": "No quality score found for PDB ID",
        },
    ]
    assert stats == expected_stats
