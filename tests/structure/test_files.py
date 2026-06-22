from pathlib import Path

import pytest

from protein_quest.structure.files import glob_structure_files, locate_structure_file, split_name_and_extension
from protein_quest.structure.types import valid_structure_file_extensions


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
        ("2y29", "2y29.bcif.gz"),
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

def test_split_name_and_extension_without_extension(tmp_path: Path):
    file = tmp_path / "filename"
    file.write_text("some text")

    filename, extension = split_name_and_extension(file.name)

    assert filename == "filename"
    assert extension == ""


def test_glob_structure_files(tmp_path: Path):
    written_files = set()
    for ext in valid_structure_file_extensions:
        file = tmp_path / f"file{ext}"
        file.write_text("some text")
        written_files.add(file)

    non_structure_file = tmp_path / "file.txt"
    non_structure_file.write_text("some text")
    nested_structure_file = tmp_path / "nested" / "file.pdb"
    nested_structure_file.parent.mkdir()
    nested_structure_file.write_text("some text")

    found_files = set(glob_structure_files(tmp_path))

    assert found_files
    assert found_files == written_files
    assert non_structure_file not in found_files
    assert nested_structure_file not in found_files