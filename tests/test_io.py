from pathlib import Path

import pytest

from protein_quest.io import (
    convert_to_cif_file,
    convert_to_cif_files,
    glob_structure_files,
    gunzip_file,
    locate_structure_file,
    read_structure,
    split_name_and_extension,
    structure2bcif,
    write_structure,
)


@pytest.mark.parametrize("extension", [".pdb", ".pdb.gz", ".cif", ".cif.gz", ".bcif"])
def test_write_structure(sample2_cif: Path, tmp_path: Path, extension: str):
    structure = read_structure(sample2_cif)
    output_file = tmp_path / f"bla{extension}"

    write_structure(structure, output_file)

    found_files = list(glob_structure_files(tmp_path))
    assert len(found_files) == 1
    assert found_files[0].name == output_file.name


def test_write_structure_invalid_extension(sample2_cif: Path, tmp_path: Path):
    structure = read_structure(sample2_cif)
    output_file = tmp_path / "bla.txt"

    with pytest.raises(ValueError, match="Unsupported file extension"):
        write_structure(structure, output_file)


@pytest.mark.parametrize("extension", [".pdb", ".pdb.gz", ".cif", ".cif.gz", ".bcif"])
def test_read_structure(sample2_cif: Path, tmp_path: Path, extension: str):
    # We only have cif as fixture, so convert to other formats first
    structure_from_cif = read_structure(sample2_cif)
    thefile = tmp_path / f"foo{extension}"
    write_structure(structure_from_cif, thefile)

    structure_from_extension = read_structure(thefile)

    assert structure_from_extension.make_minimal_pdb() == structure_from_cif.make_minimal_pdb()


@pytest.mark.parametrize(
    "pdb_id, file_name",
    [
        # extensions
        ("2y29", "2y29.cif"),
        ("2y29", "2y29.cif.gz"),
        ("2y29", "2y29.pdb"),
        ("2y29", "2y29.pdb.gz"),
        ("2y29", "pdb2y29.ent"),
        ("2y29", "pdb2y29.ent.gz"),
        ("2y29", "2y29.bcif"),
        # cases
        ("1KVm", "1KVm.cif"),
        ("1KVm", "1kvm.cif"),
        ("1KVm", "1KVM.cif"),
    ],
)
def test_locate_structure_file(tmp_path: Path, pdb_id: str, file_name: str):
    test_input_file = tmp_path / file_name
    test_input_file.write_text("fake content")
    result = locate_structure_file(tmp_path, pdb_id)

    assert result == test_input_file


def test_locate_structure_file_notfound(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="No structure file found for nonexistent_id in"):
        locate_structure_file(tmp_path, "nonexistent_id")


def test_convert_cifgz_to_cif(sample_cif: Path, tmp_path: Path):
    output_file = convert_to_cif_file(sample_cif, tmp_path, copy_method="symlink")

    assert output_file.exists()
    assert output_file.name == "3JRS_B2A.cif"
    assert output_file.stat().st_size > sample_cif.stat().st_size  # uncompressed file is larger
    assert not output_file.is_symlink()
    assert output_file.stat().st_nlink == 1  # not a hard link


def test_convert_many_cifgz_to_cif(sample_cif: Path, tmp_path: Path):
    output_files = list(convert_to_cif_files([sample_cif], tmp_path, copy_method="symlink"))

    assert len(output_files) == 1
    input_file = output_files[0][0]
    assert input_file == sample_cif
    output_file = output_files[0][1]
    assert output_file.exists()
    assert output_file.name == "3JRS_B2A.cif"
    assert output_file.stat().st_size > sample_cif.stat().st_size  # uncompressed file is larger
    assert not output_file.is_symlink()
    assert output_file.stat().st_nlink == 1  # not a hard link


def test_convert_bcif_to_cif(sample_cif: Path, tmp_path: Path):
    structure = read_structure(sample_cif)
    bcif_file = tmp_path / "bla.bcif"
    structure2bcif(structure, bcif_file)

    output_file = convert_to_cif_file(bcif_file, tmp_path, copy_method="symlink")

    assert output_file.exists()
    assert output_file.name == "bla.cif"
    assert output_file.stat().st_size > bcif_file.stat().st_size  # uncompressed file is larger
    assert not output_file.is_symlink()
    assert output_file.stat().st_nlink == 1  # not a hard link


def test_convert_cif_to_cif_in_other_dir(sample_cif: Path, tmp_path: Path):
    uncompressed_cif_file = tmp_path / "3JRS_B2A.cif"
    gunzip_file(sample_cif, uncompressed_cif_file)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    output_file = convert_to_cif_file(uncompressed_cif_file, out_dir, copy_method="hardlink")

    assert output_file.exists()
    assert out_dir in output_file.parents
    assert output_file.name == "3JRS_B2A.cif"
    assert output_file.stat() == uncompressed_cif_file.stat()  # same file
    assert not output_file.is_symlink()
    assert output_file.stat().st_nlink == 2  # hard link created


def test_convert_cif_to_cif_in_same_dir(sample_cif: Path, tmp_path: Path):
    uncompressed_cif_file = tmp_path / "3JRS_B2A.cif"
    gunzip_file(sample_cif, uncompressed_cif_file)

    output_file = convert_to_cif_file(uncompressed_cif_file, tmp_path, copy_method="hardlink")

    assert output_file.exists()
    assert output_file.name == "3JRS_B2A.cif"
    assert output_file.stat() == uncompressed_cif_file.stat()  # same file
    assert not output_file.is_symlink()
    assert output_file.stat().st_nlink == 1  # no hard link created, same file


def test_convert_pdb_to_cif(sample_cif: Path, tmp_path: Path):
    structure = read_structure(sample_cif)
    pdb_file = tmp_path / "bla.pdb"
    write_structure(structure, pdb_file)

    output_file = convert_to_cif_file(pdb_file, tmp_path, copy_method="symlink")

    assert output_file.exists()
    assert output_file.name == "bla.cif"
    assert not output_file.is_symlink()
    assert output_file.stat().st_nlink == 1  # not a hard link


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


def test_gunzip_of_non_gz_file(tmp_path: Path):
    input_file = tmp_path / "bla.zip"
    input_file.write_text("some text")

    with pytest.raises(ValueError, match="must end with \\.gz"):
        gunzip_file(input_file)


def test_split_name_and_extension_without_extension(tmp_path: Path):
    file = tmp_path / "filename"
    file.write_text("some text")

    filename, extension = split_name_and_extension(file.name)

    assert filename == "filename"
    assert extension == ""
