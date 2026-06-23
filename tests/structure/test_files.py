from pathlib import Path

import pytest

from protein_quest.structure.files import (
    LocateStructureFilesByIdResult,
    glob_structure_files,
    locate_structure_file,
    locate_structure_files_by_id,
    split_name_and_extension,
)
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


class TestLocateStructureFilesByIdResult:
    found_cases = (
        ("2y29", "2y29.cif"),
        ("2Y29", "2y29.cif"),  # case-insensitive id
        ("2y29", "2Y29.cif"),  # case-insensitive file name
        ("2y29", "2y29.cif.gz"),  # different extension
        ("2y29", "2y29.pdb"),
        ("2y29", "2y29.pdb.gz"),
        ("2y29", "2y29.ent"),
        ("2y29", "2y29.ent.gz"),
        ("2y29", "2y29.bcif"),
        ("2y29", "2y29.bcif.gz"),
        ("2y29", "2y29_B2A.cif"),  # chain filtered aka stuff between id and extension
        ("2y29", "pdb2y29.cif"),  # prefix
        ("2y29", "2y29_updated.cif"),  # postfix
    )

    @pytest.mark.parametrize(
        "pdb_id, file_name",
        found_cases,
        ids=[f"id={pdb_id},fn={file_name}" for pdb_id, file_name in found_cases],
    )
    def test_found(self, tmp_path: Path, pdb_id: str, file_name: str):
        pdb_ids = {pdb_id}
        fn = tmp_path / file_name
        fn.write_text("some text")

        result = locate_structure_files_by_id(
            pdb_ids,
            input_dir=tmp_path,
        )

        expected = LocateStructureFilesByIdResult(
            found={(pdb_id, fn)},
        )
        assert result == expected

    def test_not_found(self, tmp_path: Path):
        pdb_ids = {"nonexistent_id"}

        result = locate_structure_files_by_id(
            pdb_ids,
            input_dir=tmp_path,
        )

        expected = LocateStructureFilesByIdResult(
            not_found={"nonexistent_id"},
        )
        assert result == expected

    def test_extras(self, tmp_path: Path):
        extra_file = tmp_path / "extra.cif"
        extra_file.write_text("some text")

        result = locate_structure_files_by_id(
            set(),
            input_dir=tmp_path,
        )

        expected = LocateStructureFilesByIdResult(
            extras={extra_file},
        )
        assert result == expected

    def test_real_found(self, tmp_path: Path, all_cifs: list[Path]):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        for fixture_name in all_cifs:
            (input_dir / fixture_name.name).symlink_to(fixture_name)

        pdb_ids = {"2y29", "1un5", "3jrs", "6o5i", "1amb", "8w77"}
        result = locate_structure_files_by_id(pdb_ids, input_dir)

        expected = LocateStructureFilesByIdResult(
            found={
                ("2y29", input_dir / "2Y29.cif.gz"),
                ("1un5", input_dir / "1un5.cif.gz"),
                ("3jrs", input_dir / "3JRS_B2A.cif.gz"),
                ("6o5i", input_dir / "6O5I.cif.gz"),
                ("1amb", input_dir / "1amb_updated.cif.gz"),
                ("8w77", input_dir / "8w77_updated.cif.gz"),
            },
            extras={input_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"},
        )
        assert result == expected

    def test_found_many(self, tmp_path: Path):
        nr_structures = 100
        # On my machine 10_000 took 2s and 100 took 0.02s
        pdb_ids = {f"{i:04d}" for i in range(nr_structures)}
        for i in range(nr_structures):
            file = tmp_path / f"{i:04d}.cif"
            file.write_text("some text")

        result = locate_structure_files_by_id(pdb_ids, tmp_path)

        expected = LocateStructureFilesByIdResult(
            found={(f"{i:04d}", tmp_path / f"{i:04d}.cif") for i in range(nr_structures)},
        )
        assert result == expected

    def test_files_with_same_id(self, tmp_path: Path):
        pdb_id = "2y29"
        file1 = tmp_path / "2y29_A2A.cif"
        file1.write_text("some text")
        file2 = tmp_path / "2y29_B2A.cif"
        file2.write_text("some text")

        result = locate_structure_files_by_id({pdb_id}, tmp_path)

        expected = LocateStructureFilesByIdResult(
            found={(pdb_id, file1), (pdb_id, file2)},
        )
        assert result == expected

    def test_file_with_multiple_ids(self, tmp_path: Path):
        file1 = tmp_path / "1aaa_2bbb.cif"
        file1.write_text("some text")

        result = locate_structure_files_by_id({"1aaa", "2bbb"}, tmp_path)

        expected = LocateStructureFilesByIdResult(found={("1aaa", file1), ("2bbb", file1)})
        assert result == expected
