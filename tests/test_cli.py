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
