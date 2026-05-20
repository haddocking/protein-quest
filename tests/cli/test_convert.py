"""Convert CLI tests for protein-quest."""

from pathlib import Path

import pytest

from protein_quest.cli import main


def test_convert_structures_to_cifgz(sample_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    input_file = input_dir / sample_cif.name
    input_file.write_bytes(sample_cif.read_bytes())
    output_dir = tmp_path / "output"

    main(
        [
            "convert",
            "structures",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--output-format",
            ".cif.gz",
        ]
    )

    assert (output_dir / "3JRS_B2A.cif.gz").exists()
    captured = capsys.readouterr()
    assert ".cif.gz" in captured.err
