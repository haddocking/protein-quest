import csv
from pathlib import Path
from textwrap import dedent

import pytest

from protein_quest.cli import main, make_parser


def test_make_parser_help(capsys: pytest.CaptureFixture[str]):
    in_args = ["--help"]
    parser = make_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(in_args)

    captured = capsys.readouterr()
    assert "Protein Quest CLI" in captured.out


@pytest.mark.vcr
def test_search_uniprot(capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture):
    argv = [
        "search",
        "uniprot",
        "--taxon-id",
        "9606",
        "--reviewed",
        "--limit",
        "1",
        "-",
    ]

    main(argv)

    captured = capsys.readouterr()
    expected = "A0A024R1R8\n"
    assert captured.out == expected
    assert "Searching for UniProt accessions" in captured.err
    assert "Found 1 UniProt accessions, written to <stdout>" in captured.err
    assert "There may be more results available" in caplog.text


@pytest.mark.vcr
def test_search_pdbe(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("P00811\n")
    output_file = tmp_path / "pdbe_results.csv"
    argv = [
        "search",
        "pdbe",
        "--limit",
        "150",
        "--min-residues",
        "360",  # P00811 has 377 residues and 5 full PDB entries
        str(input_text),
        str(output_file),
    ]

    main(argv)

    result = output_file.read_text()
    expected = dedent("""\
        uniprot_accession,pdb_id,method,resolution,uniprot_chains,chain,chain_length
        P00811,9C6P,X-Ray_Crystallography,1.66,A/B=1-377,A,377
        P00811,9C81,X-Ray_Crystallography,1.7,A/B=1-377,A,377
        P00811,9C83,X-Ray_Crystallography,2.9,A/B=1-377,A,377
        P00811,9C84,X-Ray_Crystallography,1.7,A/B=1-377,A,377
        P00811,9DHL,X-Ray_Crystallography,1.88,A/B=1-377,A,377
        """)
    assert result == expected

    captured = capsys.readouterr()
    assert "Finding PDB entries for 1 uniprot accessions" in captured.err
    assert "Before filtering found 120 PDB entries for 1 uniprot accessions." in captured.err
    assert "After filtering on chain length (360, None) remained 5 PDB entries for 1 uniprot" in captured.err
    assert "Written to " in captured.err


@pytest.mark.vcr
def test_search_uniprot_details(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("P05067\nA0A0B5AC95\n")
    output_csv = tmp_path / "uniprot_details.csv"
    argv = [
        "search",
        "uniprot-details",
        str(input_text),
        str(output_csv),
    ]

    main(argv)

    result = output_csv.read_text()
    expected = dedent("""\
        uniprot_accession,uniprot_id,sequence_length,reviewed,protein_name,taxon_id,taxon_name
        A0A0B5AC95,INS1A_CONGE,115,True,Con-Ins G1a,6491,Conus geographus
        P05067,A4_HUMAN,770,True,Amyloid-beta precursor protein,9606,Homo sapiens
        """)
    assert result == expected
    captured = capsys.readouterr()
    assert "Retrieving UniProt entry details for 2 uniprot accessions" in captured.err
    assert "Retrieved details for 2 UniProt entries, written to " in captured.err


@pytest.mark.vcr
def test_search_alphafold(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("P00811\n")
    output_file = tmp_path / "af_results.csv"

    argv = [
        "search",
        "alphafold",
        str(input_text),
        str(output_file),
    ]

    main(argv)

    result = output_file.read_text()

    expected = dedent("""\
        uniprot_accession,af_id
        P00811,P00811
        """)
    assert result == expected

    captured = capsys.readouterr()
    assert "Finding AlphaFold entries for 1 uniprot accessions" in captured.err
    assert "Found 1 AlphaFold entries, written to " in captured.err


def test_filter_chain_happy_path(sample2_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    chains_fn = tmp_path / "chains.csv"
    chains_fn.write_text("pdb_id,chain\n2Y29,A\n")

    argv = [
        "filter",
        "chain",
        str(chains_fn),
        str(sample2_cif.parent),
        str(tmp_path),
    ]

    main(argv)

    output_file = tmp_path / "2Y29_A2A.cif.gz"
    assert output_file.exists()

    captured = capsys.readouterr()
    assert "Wrote 1 single-chain PDB/mmCIF files to" in captured.err


def test_filter_chain_input_file_notfound(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
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

    with pytest.raises(SystemExit):
        main(argv)

    assert not any(output_dir.iterdir())

    captured = capsys.readouterr()
    assert "No structure file found for 2Y29" in captured.err


def test_filter_residue(sample_cif: Path, sample2_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
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
