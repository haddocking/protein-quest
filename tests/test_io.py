from pathlib import Path

import pytest

from protein_quest.io import glob_structure_files, locate_structure_file, read_structure, write_structure


@pytest.mark.parametrize("extension", [".pdb", ".pdb.gz", ".cif", ".cif.gz", ".bcif"])
def test_write_structure(sample2_cif: Path, tmp_path: Path, extension: str):
    structure = read_structure(sample2_cif)
    output_file = tmp_path / f"bla{extension}"

    write_structure(structure, output_file)

    found_files = list(glob_structure_files(tmp_path))
    assert len(found_files) == 1
    assert found_files[0].name == output_file.name


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
