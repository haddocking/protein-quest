from pathlib import Path
from textwrap import dedent

import pytest

from protein_quest.cli import main


@pytest.mark.vcr
def test_search_uniprot(capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture):
    """Test search uniprot command."""
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
def test_search_uniprot_with_provenance(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Test search uniprot with provenance recording."""
    output_file = Path("uniprot_accessions.txt")
    argv = [
        "search",
        "uniprot",
        "--taxon-id",
        "9606",
        "--reviewed",
        "--limit",
        "1",
        str(output_file),
        "--prov",
    ]
    monkeypatch.chdir(tmp_path)

    main(argv)

    assert output_file.exists()
    prov_file = tmp_path / "ro-crate-metadata.json"
    assert prov_file.exists()
    body = prov_file.read_text()
    assert '"@type": "CreateAction"' in body
    assert '"name": "protein-quest"' in body
    assert body.count("uniprot_accessions.txt") >= 4


@pytest.mark.vcr
def test_search_pdbe(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test search pdbe command."""
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


@pytest.mark.default_cassette("test_search_pdbe.yaml")
@pytest.mark.vcr
def test_search_pdbe_top_resolution_per_accession(tmp_path: Path):
    """Test search pdbe command with top-resolution filtering."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("P00811\n")
    output_file = tmp_path / "pdbe_results.csv"
    argv = [
        "search",
        "pdbe",
        "--limit",
        "150",
        "--min-residues",
        "360",
        "--top-resolution-per-uniprot-accession",
        "2",
        str(input_text),
        str(output_file),
    ]

    main(argv)

    result = output_file.read_text()
    expected = dedent("""\
        uniprot_accession,pdb_id,method,resolution,uniprot_chains,chain,chain_length
        P00811,9C6P,X-Ray_Crystallography,1.66,A/B=1-377,A,377
        P00811,9C81,X-Ray_Crystallography,1.7,A/B=1-377,A,377
        """)
    assert result == expected


@pytest.mark.default_cassette("test_search_pdbe_bad_chain_length.yaml")
@pytest.mark.vcr
def test_search_pdbe_bad_chain_length(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    """Test search pdbe with bad chain length."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("Q9NTW7\n")
    output_file = tmp_path / "pdbe_results.csv"
    argv = [
        "search",
        "pdbe",
        str(input_text),
        str(output_file),
    ]

    main(argv)

    assert len(output_file.read_text()) == 159

    assert "Could not determine chain length for " in caplog.text
    assert "Q9NTW7 / 1X5W chain A from 'A=-'" in caplog.text
    assert "No chain length for this entry." in caplog.text


@pytest.mark.default_cassette("test_search_pdbe_bad_chain_length.yaml")
@pytest.mark.vcr
def test_search_pdbe_bad_chain_length_with_min_nokeep(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    """Test search pdbe with bad chain length and min residues without keep."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("Q9NTW7\n")
    output_file = tmp_path / "pdbe_results.csv"
    argv = [
        "search",
        "pdbe",
        "--min-residues",
        "42",
        str(input_text),
        str(output_file),
    ]

    main(argv)

    assert len(output_file.read_text()) == 122

    log = caplog.text
    assert (
        "Filtering out PDB entry '1X5W' belonging to uniprot accession 'Q9NTW7' due to invalid chain length from 'A=-'"
        in log
    )


@pytest.mark.default_cassette("test_search_pdbe_bad_chain_length.yaml")
@pytest.mark.vcr
def test_search_pdbe_bad_chain_length_with_min_keep(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    """Test search pdbe with bad chain length and min residues with keep."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("Q9NTW7\n")
    output_file = tmp_path / "pdbe_results.csv"
    argv = [
        "search",
        "pdbe",
        "--keep-invalid",
        "--min-residues",
        "42",
        str(input_text),
        str(output_file),
    ]

    main(argv)

    assert len(output_file.read_text()) == 159

    log = caplog.text
    assert "for completeness not filtering it out" in log
    assert "Could not determine chain length for " in log
    assert "Q9NTW7 / 1X5W chain A from 'A=-'" in log
    assert "No chain length for this entry." in log


@pytest.mark.vcr
def test_search_structure(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test search structure command."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("Q9NTW7\n")
    output_file = tmp_path / "structure_results.csv"
    raw_output_file = tmp_path / "structure.results.json"
    argv = [
        "search",
        "structure",
        "--limit",
        "1",
        "--source",
        "alphafold",
        "--source",
        "alphafill",
        "--raw",
        str(raw_output_file),
        str(input_text),
        str(output_file),
    ]

    main(argv)

    assert len(output_file.read_text()) == 260
    assert len(raw_output_file.read_text()) == 3522

    captured = capsys.readouterr()
    assert "Finding structures for 1 uniprot accessions" in captured.err
    assert "Written raw results to" in captured.err
    assert "Found 2 structures, written to" in captured.err


@pytest.mark.vcr
def test_search_structure_all_sources(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test search structure with all sources."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("Q9NTW7\n")
    output_file = tmp_path / "structure_results.csv"
    raw_output_file = tmp_path / "structure.results.json"
    argv = [
        "search",
        "structure",
        "--limit",
        "1",
        "--source",
        "all",
        "--raw",
        str(raw_output_file),
        str(input_text),
        str(output_file),
    ]

    main(argv)

    output_content = output_file.read_text()
    raw_content = raw_output_file.read_text()
    # Check that output contains expected data (length may vary due to API changes)
    assert "uniprot_accession" in output_content
    assert "provider" in output_content
    assert "Q9NTW7" in output_content
    # Check that raw output exists and contains expected data
    assert '"uniprot_entry"' in raw_content
    assert '"structures"' in raw_content

    captured = capsys.readouterr()
    assert "Finding structures for 1 uniprot accessions" in captured.err
    assert "Written raw results to" in captured.err
    assert "Found 5 structures, written to" in captured.err


@pytest.mark.vcr
def test_search_emdb(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test search emdb command."""
    input_text = tmp_path / "uniprot_accessions.txt"
    input_text.write_text("O14646\n")
    output_file = tmp_path / "emdbs.csv"
    argv = [
        "search",
        "emdb",
        "--limit",
        "150",
        str(input_text),
        str(output_file),
    ]

    main(argv)

    result = output_file.read_text()
    expected = dedent("""\
        uniprot_accession,emdb_id
        O14646,EMD-47841
        O14646,EMD-49406
        """)
    assert result == expected
    captured = capsys.readouterr()
    assert "Finding EMDB entries for 1 uniprot accessions" in captured.err
    assert "Found 2 EMDB entries, written to" in captured.err
    assert str(output_file) in captured.err


@pytest.mark.vcr
def test_search_uniprot_details(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test search uniprot-details command."""
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
    """Test search alphafold command."""
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
