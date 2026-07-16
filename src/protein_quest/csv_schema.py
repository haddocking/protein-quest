"""Shared CSV schema definitions for protein-quest CLI commands."""

from collections import namedtuple

ChainUniprotPair = namedtuple("ChainUniprotPair", ["chain_id", "uniprot_accession"])

# Standard fieldnames for pair rows (chain_id + uniprot_accession).
# When a PDB ID is applicable, include pdb_id as well.
PAIR_FIELDNAMES = ("uniprot_accession", "chain_id", "pdb_id")
