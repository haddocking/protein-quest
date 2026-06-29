"""Shared type declarations for structure modules."""

from dataclasses import dataclass
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
valid_structure_file_extensions: set[str] = set(get_args(StructureFileExtensions))
"""Set of valid structure file extensions."""

CifOutputFormat = Literal[".cif", ".cif.gz"]
"""Supported output formats for CIF conversion."""
cif_output_formats: set[str] = set(get_args(CifOutputFormat))
"""Set of valid CIF conversion output formats."""

StructureMethod = Literal["EM", "NMR", "Predicted", "X-ray", "Other"]
"""Represents the method used to determine the structure."""

Pdb2UniprotMapping = dict[str, set[tuple[str, str]]]
"""Dictionary mapping PDB ID to set of tuples containing chain and UniProt accession.

The chain name is in the 'auth' [chain ID system][protein_quest.structure.chains.ChainIdSystem].
"""


@dataclass(frozen=True)
class StructRefSeq:
    """Collapsed `_struct_ref_seq` alignment information for one chain."""

    uniprot_accession: str
    uniprot_start: int
    uniprot_end: int
    chain_id: str
    sequence_identity: float
    aligned_residue_count: int
