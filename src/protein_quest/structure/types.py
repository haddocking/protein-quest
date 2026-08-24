"""Shared type declarations for structure modules."""

from typing import Literal, get_args

StructureFileExtensions = Literal[
    ".pdb",
    ".pdb.gz",
    ".ent",
    ".ent.gz",
    ".cif",
    ".cif.gz",
    ".bcif",
    ".bcif.gz",
]
"""Type of supported structure file extensions."""
ordered_structure_file_extensions: tuple[str, ...] = get_args(StructureFileExtensions)
"""Ordered structure file extensions in preferred lookup/write precedence."""
valid_structure_file_extensions: frozenset[str] = frozenset(ordered_structure_file_extensions)
"""Set of valid structure file extensions."""

CifOutputFormat = Literal[".cif", ".cif.gz"]
"""Supported output formats for CIF conversion."""
cif_output_formats: set[str] = set(get_args(CifOutputFormat))
"""Set of valid CIF conversion output formats."""

StructureMethod = Literal["EM", "NMR", "Predicted", "X-ray", "Other"]
"""Represents the method used to determine the structure."""
