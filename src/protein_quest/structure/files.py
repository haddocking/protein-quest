"""Path and filename helpers for structure files."""

from collections.abc import Generator
from dataclasses import dataclass, field
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


@dataclass(frozen=True, slots=True)
class LocateStructureFilesByIdResult:
    """Result of locating structure files by PDB ID.

    Attributes:
        found: A set of PDB IDs and their located structure file path.
           A PDB ID may be associated with multiple structure files and
           a structure file may be associated with multiple PDB IDs.
        not_found: A set of PDB IDs that could not be located.
        extras: A set of structure files in the input directory that were not associated with any of the provided IDs.
    """

    found: set[tuple[str, Path]] = field(default_factory=set)
    not_found: set[str] = field(default_factory=set)
    extras: set[Path] = field(default_factory=set)


def locate_structure_files_by_id(ids: set[str], input_dir: Path) -> LocateStructureFilesByIdResult:
    """Locate structure files for a set of PDB IDs in the specified directory.

    Use presence of ID in filename to associate files with IDs.

    Args:
        ids: A set of PDB IDs to locate.
        input_dir: The directory to search for structure files.

    Returns:
        A LocateStructureFilesByIdResult containing found files, not found IDs, and extra files found.
    """
    all_files = list(glob_structure_files(input_dir))
    # nested loop is not very efficient, but works and is readable
    ids_lower = {pdb_id.lower(): pdb_id for pdb_id in ids}
    found = set()
    for file in all_files:
        file_name_lower = file.name.lower()
        for pdb_id_lower, pdb_id in ids_lower.items():
            if pdb_id_lower in file_name_lower:
                found.add((pdb_id, file))
                continue

    not_found = ids - {pdb_id for pdb_id, _ in found}
    extras = set(all_files) - {file for _, file in found}
    return LocateStructureFilesByIdResult(
        found=found,
        not_found=not_found,
        extras=extras,
    )
