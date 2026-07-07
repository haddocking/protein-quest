"""Path and filename helpers for structure files."""

from collections.abc import Generator
from pathlib import Path

from protein_quest.structure.types import valid_structure_file_extensions


def split_name_and_extension(name: str) -> tuple[str, str]:
    """Split a filename into its name and extension.

    `.gz` is considered part of the extension if present.

    Examples:
        Some example usages.

        >>> from protein_quest.structure.files import split_name_and_extension
        >>> split_name_and_extension("1234.pdb")
        ('1234', '.pdb')
        >>> split_name_and_extension("1234.pdb.gz")
        ('1234', '.pdb.gz')

    Args:
        name: The filename to split.

    Returns:
        A tuple containing the name and the extension."""
    ext = ""
    if name.endswith(".gz"):
        ext = ".gz"
        name = name.removesuffix(".gz")
    i = name.rfind(".")
    if 0 < i < len(name) - 1:
        ext = name[i:] + ext
        name = name[:i]
    return name, ext


def locate_structure_file(root: Path, pdb_id: str) -> Path:
    """Locate a structure file for a given PDB ID in the specified directory.

    Uses [StructureFileExtensions][protein_quest.structure.types.StructureFileExtensions] as potential extensions.
    Also tries different casing of the PDB ID.

    Args:
        root: The root directory to search in.
        pdb_id: The PDB ID to locate.

    Returns:
        The path to the located structure file.

    Raises:
        FileNotFoundError: If no structure file is found for the given PDB ID."""
    for ext in valid_structure_file_extensions:
        candidates = (
            root / f"{pdb_id}{ext}",
            root / f"{pdb_id.lower()}{ext}",
            root / f"{pdb_id.upper()}{ext}",
            root / f"pdb{pdb_id.lower()}{ext}",
            root / f"{pdb_id.lower()}_updated{ext}",
        )
        for candidate in candidates:
            if candidate.exists():
                return candidate
    msg = f"No structure file found for {pdb_id} in {root}"
    raise FileNotFoundError(msg)


def glob_structure_files(input_dir: Path) -> Generator[Path]:
    """Glob for structure files in a directory.

    Uses [StructureFileExtensions][protein_quest.structure.types.StructureFileExtensions] as valid extensions.
    Does not search recursively.

    Args:
        input_dir: The input directory to search for structure files.

    Yields:
        Paths to the found structure files."""
    for ext in valid_structure_file_extensions:
        yield from input_dir.glob(f"*{ext}")
