from protein_quest.structure.types import valid_structure_file_extensions


def test_valid_structure_file_extensions():
    assert valid_structure_file_extensions == {
        ".cif",
        ".cif.gz",
        ".bcif",
        ".bcif.gz",
        ".pdb",
        ".pdb.gz",
        ".ent",
        ".ent.gz",
    }
