from pathlib import Path

import pytest

from protein_quest.utils import copyfile


def test_copyfile_copy(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    copyfile(src, dst, "copy")
    assert dst.read_text() == "Hello, World!"
    assert not dst.is_symlink()

def test_copyfile_symlink(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    copyfile(src, dst, "symlink")
    assert dst.read_text() == "Hello, World!"
    assert dst.is_symlink()
    assert dst.resolve() == src

def test_copyfile_invalid_method(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    with pytest.raises(ValueError, match="Unknown method"):
        copyfile(src, dst, "invalid") # pyright: ignore[reportArgumentType]
