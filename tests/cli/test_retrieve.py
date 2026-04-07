"""Retrieve CLI tests for protein-quest."""

from pathlib import Path

import pytest

from protein_quest.cli import main


@pytest.mark.vcr
def test_retrieve_structure_happy_path(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Test retrieve structure command."""
    input_csv = tmp_path / "structures.csv"
    input_csv.write_text(
        "provider,model_identifier,model_url,model_format\n"
        "swissmodel,Q9NTW7_329-603:5v3m.1.C,https://swissmodel.expasy.org/3d-beacons/uniprot/Q9NTW7.cif?range=329-603&template=5v3m.1.C&provider=swissmodel,MMCIF\n"
    )
    output_dir = tmp_path / "downloads"

    main(["retrieve", "structure", "--no-cache", str(input_csv), str(output_dir)])

    assert (output_dir / "swissmodel~Q9NTW7_329-603:5v3m.1.C.cif.gz").exists()
    captured = capsys.readouterr()
    assert "downloaded=1" in captured.err
    assert "converted=0" in captured.err
    assert "cached=0" in captured.err
