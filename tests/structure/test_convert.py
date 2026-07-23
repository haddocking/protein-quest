from pathlib import Path

import pytest

from protein_quest.structure.convert import convert_to_cif_file, convert_to_cif_files, read_structure
from protein_quest.structure.formats import gunzip_file, structure2bcif, write_structure


def test_convert_cifgz_to_cif(sample_cif: Path, tmp_path: Path):
    stats = convert_to_cif_file(sample_cif, tmp_path, copy_method="symlink")

    assert stats.output_file.exists()
    assert stats.output_file.name == "3jrs_updated_B2A.cif"
    assert stats.output_file.stat().st_size > sample_cif.stat().st_size  # uncompressed file is larger
    assert not stats.output_file.is_symlink()
    assert stats.output_file.stat().st_nlink == 1  # not a hard link


def test_convert_many_cifgz_to_cif(sample_cif: Path, tmp_path: Path):
    results = list(convert_to_cif_files([sample_cif], tmp_path, copy_method="symlink"))

    assert len(results) == 1
    stats = results[0]
    assert stats.input_file == sample_cif
    assert stats.output_file.exists()
    assert stats.output_file.name == "3jrs_updated_B2A.cif"
    assert stats.output_file.stat().st_size > sample_cif.stat().st_size  # uncompressed file is larger
    assert not stats.output_file.is_symlink()
    assert stats.output_file.stat().st_nlink == 1  # not a hard link


def test_convert_bcif_to_cif(sample_cif: Path, tmp_path: Path):
    structure = read_structure(sample_cif)
    bcif_file = tmp_path / "bla.bcif"
    structure2bcif(structure, bcif_file)

    stats = convert_to_cif_file(bcif_file, tmp_path, copy_method="symlink")

    assert stats.output_file.exists()
    assert stats.output_file.name == "bla.cif"
    assert stats.output_file.stat().st_size > bcif_file.stat().st_size  # uncompressed file is larger
    assert not stats.output_file.is_symlink()
    assert stats.output_file.stat().st_nlink == 1  # not a hard link


def test_convert_cif_to_cif_in_other_dir(sample_cif: Path, tmp_path: Path):
    uncompressed_cif_file = tmp_path / "3JRS_B2A.cif"
    gunzip_file(sample_cif, uncompressed_cif_file)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    stats = convert_to_cif_file(uncompressed_cif_file, out_dir, copy_method="hardlink")

    assert stats.output_file.exists()
    assert out_dir in stats.output_file.parents
    assert stats.output_file.name == "3JRS_B2A.cif"
    assert stats.output_file.stat() == uncompressed_cif_file.stat()  # same file
    assert not stats.output_file.is_symlink()
    assert stats.output_file.stat().st_nlink == 2  # hard link created


def test_convert_cif_to_cif_in_same_dir(sample_cif: Path, tmp_path: Path):
    uncompressed_cif_file = tmp_path / "3JRS_B2A.cif"
    gunzip_file(sample_cif, uncompressed_cif_file)

    stats = convert_to_cif_file(uncompressed_cif_file, tmp_path, copy_method="hardlink")

    assert stats.output_file.exists()
    assert stats.output_file.name == "3JRS_B2A.cif"
    assert stats.output_file.stat() == uncompressed_cif_file.stat()  # same file
    assert not stats.output_file.is_symlink()
    assert stats.output_file.stat().st_nlink == 1  # no hard link created, same file


def test_convert_cif_to_cifgz(sample_cif: Path, tmp_path: Path):
    uncompressed_cif_file = tmp_path / "3JRS_B2A.cif"
    gunzip_file(sample_cif, uncompressed_cif_file)

    stats = convert_to_cif_file(uncompressed_cif_file, tmp_path, copy_method="symlink", output_format=".cif.gz")

    assert stats.output_file.exists()
    assert stats.output_file.name == "3JRS_B2A.cif.gz"
    assert stats.output_file.stat().st_size > 0
    assert (
        read_structure(stats.output_file).make_pdb_string() == read_structure(uncompressed_cif_file).make_pdb_string()
    )


def test_convert_many_cif_to_cifgz(sample_cif: Path, tmp_path: Path):
    uncompressed_cif_file = tmp_path / "3JRS_B2A.cif"
    gunzip_file(sample_cif, uncompressed_cif_file)

    results = list(
        convert_to_cif_files([uncompressed_cif_file], tmp_path, copy_method="symlink", output_format=".cif.gz")
    )

    assert len(results) == 1
    stats = results[0]
    assert stats.input_file == uncompressed_cif_file
    assert stats.output_file.exists()
    assert stats.output_file.name == "3JRS_B2A.cif.gz"
    assert (
        read_structure(stats.output_file).make_pdb_string() == read_structure(uncompressed_cif_file).make_pdb_string()
    )


def test_convert_pdb_to_cif(sample_cif: Path, tmp_path: Path):
    structure = read_structure(sample_cif)
    pdb_file = tmp_path / "bla.pdb"
    write_structure(structure, pdb_file)

    stats = convert_to_cif_file(pdb_file, tmp_path, copy_method="symlink")

    assert stats.output_file.exists()
    assert stats.output_file.name == "bla.cif"
    assert not stats.output_file.is_symlink()
    assert stats.output_file.stat().st_nlink == 1  # not a hard link


def test_convert_invalid_extension(tmp_path: Path):
    txt_file = tmp_path / "bla.txt"
    txt_file.write_text("some text")

    with pytest.raises(ValueError, match="Unsupported file extension"):
        convert_to_cif_file(txt_file, tmp_path, copy_method="symlink")


def test_gunzip_file_with_remove_original(sample_cif: Path, tmp_path: Path):
    sample_cif_copy = tmp_path / sample_cif.name
    sample_cif_copy.write_bytes(sample_cif.read_bytes())

    output_file = tmp_path / "3JRS_B2A.cif"
    result_file = gunzip_file(sample_cif_copy, output_file, keep_original=False)

    assert result_file == output_file
    assert output_file.exists()
    assert not sample_cif_copy.exists()
