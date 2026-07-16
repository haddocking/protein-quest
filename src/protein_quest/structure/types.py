"""Shared type declarations for structure modules."""

from dataclasses import dataclass
from typing import Literal, get_args

from protein_quest.csv_schema import ChainUniprotPair

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

Pdb2RawPairs = dict[str, set[ChainUniprotPair]]
"""Dictionary mapping PDB ID to set of raw (chain, UniProt accession) pairs.

These are simple pairs without residue-range information.
The chain name is in the 'auth' [chain ID system][protein_quest.structure.chains.ChainIdSystem].
"""


@dataclass(frozen=True)
class StructRefSeq:
    """Collapsed `_struct_ref_seq` alignment information for one chain.

    Attributes:

        uniprot_accession: The UniProt accession.
        uniprot_start: The start position of the alignment on the UniProt sequence.
        uniprot_end: The end position of the alignment on the UniProt sequence.
        chain_id: The chain ID in the 'auth' [chain ID system][protein_quest.structure.chains.ChainIdSystem].
        sequence_identity: The sequence identity of the alignment.
        aligned_residue_count: The number of aligned residues in the alignment.
    """

    uniprot_accession: str
    uniprot_start: int
    uniprot_end: int
    chain_id: str
    sequence_identity: float
    aligned_residue_count: int
