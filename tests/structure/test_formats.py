import gzip
from pathlib import Path

import pytest

from protein_quest.structure.formats import (
    gunzip_file,
    read_structure,
    structure2bcif,
    structure2cifgz,
    write_structure,
)
from protein_quest.structure.types import valid_structure_file_extensions


@pytest.mark.parametrize("extension", valid_structure_file_extensions)
def test_write_structure(sample2_cif: Path, tmp_path: Path, extension: str):
    structure = read_structure(sample2_cif)
    output_file = tmp_path / f"bla{extension}"

    write_structure(structure, output_file)

    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_write_structure_invalid_extension(sample2_cif: Path, tmp_path: Path):
    structure = read_structure(sample2_cif)
    output_file = tmp_path / "bla.txt"

    with pytest.raises(ValueError, match="Unsupported file extension"):
        write_structure(structure, output_file)


@pytest.mark.parametrize("extension", valid_structure_file_extensions)
def test_read_structure(sample2_cif: Path, tmp_path: Path, extension: str):
    # We only have cif as fixture, so convert to other formats first
    structure_from_cif = read_structure(sample2_cif)
    thefile = tmp_path / f"foo{extension}"
    write_structure(structure_from_cif, thefile)

    structure_from_extension = read_structure(thefile)

    assert structure_from_extension.make_pdb_string() == structure_from_cif.make_pdb_string()


def test_structure2cifgz(sample2_cif: Path):
    structure = read_structure(sample2_cif)

    cif_gz_bytes = structure2cifgz(structure)

    assert isinstance(cif_gz_bytes, bytes)
    assert len(cif_gz_bytes) > 0
    observed_cif_text = gzip.decompress(cif_gz_bytes).decode("utf-8")
    assert observed_cif_text.startswith("data_")


def test_gunzip_of_non_gz_file(tmp_path: Path):
    input_file = tmp_path / "bla.zip"
    input_file.write_text("some text")

    with pytest.raises(ValueError, match="must end with \\.gz"):
        gunzip_file(input_file)

@pytest.mark.parametrize(
    ["output_fn"],
    [
        ("em_structure.cif",),
        ("em_structure.cif.gz",),
    ],
)
def test_em_structure_retains_resolution(em_cif: Path, tmp_path: Path, output_fn: str):
    structure = read_structure(em_cif)
    assert structure.resolution == 3.61
    output_file = tmp_path / output_fn

    write_structure(structure, output_file)

    written_structure = read_structure(output_file)

    assert written_structure.resolution == 3.61
