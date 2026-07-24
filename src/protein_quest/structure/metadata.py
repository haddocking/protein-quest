"""Metadata extraction and ranking helpers for structures."""

import logging
from dataclasses import dataclass
from pathlib import Path

import gemmi
from gemmi import Structure

from protein_quest.structure.chains import (
    find_chain_in_structure,
    get_label2auth_chains,
    nr_of_residues_in_total,
)
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.types import StructureMethod
from protein_quest.structure.uniprot import (
    ChainUniprotPair,
    structure_to_uniprot,
)

logger = logging.getLogger(__name__)


def _metadata_without_uniprot(
    structure: gemmi.Structure,
    *,
    total_residue_count: int,
) -> "StructureMetadata":
    label2auth_chains = get_label2auth_chains(structure)
    if not label2auth_chains:
        msg = f"No chains found in structure {structure.name}"
        raise ValueError(msg)
    # Use first chain in alphabetical order if no UniProt mapping is found
    label_chain = min(label2auth_chains.keys())
    auth_chain = label2auth_chains[label_chain]
    chain = find_chain_in_structure(structure, auth_chain)
    if chain is None:
        raise ChainNotFoundError(auth_chain, None, set(label2auth_chains.values()))
    chain_length = len(chain)
    return StructureMetadata(
        id=structure.name,
        uniprot_accession=None,
        resolution=structure.resolution,
        total_residue_count=total_residue_count,
        is_alphafold=_source_is_alphafold(structure),
        uniprot_start=0,
        uniprot_end=0,
        sequence_identity=0.0,
        chain_length=chain_length,
        auth_chain=auth_chain,
        label_chain=label_chain,
        method=_structure_method(structure),
    )


@dataclass(frozen=True)
class StructureMetadata:
    """Metadata extracted from a structure file for ranking and grouping.

    Attributes:
        id: The structure ID.
        uniprot_accession: The UniProt accession if available and only one, otherwise None.
        resolution: The resolution of the structure in Angstroms.
        total_residue_count: The total number of residues in the structure.
        is_alphafold: True if the structure was predicted by AlphaFold, otherwise False.
        uniprot_start: The start position of the UniProt sequence in the structure.
        uniprot_end: The end position of the UniProt sequence in the structure.
        sequence_identity: The sequence identity between the structure and the UniProt sequence.
        chain_length: The length of the chain in the structure.
        auth_chain: The chain in 'auth' [id system][protein_quest.structure.chains.ChainIdSystem].
        label_chain: The chain in 'label' [id system][protein_quest.structure.chains.ChainIdSystem].
        method: The experimental method used to determine the structure.
    """

    id: str
    uniprot_accession: str | None
    resolution: float
    total_residue_count: int
    is_alphafold: bool
    uniprot_start: int
    uniprot_end: int
    sequence_identity: float
    chain_length: int
    auth_chain: str
    label_chain: str
    method: StructureMethod

    def __post_init__(self):
        object.__setattr__(self, "resolution", round(self.resolution, 3))
        object.__setattr__(self, "sequence_identity", round(self.sequence_identity, 3))


def _experimental_method(structure: gemmi.Structure) -> str | None:
    try:
        return structure.info["_exptl.method"].lower()
    except KeyError:
        return None


def _software_names(structure: gemmi.Structure) -> set[str]:
    return {s.name for s in structure.meta.software}


def _is_prediction_software(structure: gemmi.Structure) -> bool:
    software_names = _software_names(structure)
    prediction_software_names = {"AlphaFold", "alphafill"}
    return software_names.intersection(prediction_software_names) != set()


def _source_is_alphafold(structure: gemmi.Structure) -> bool:
    software_names = _software_names(structure)
    return "AlphaFold" in software_names and "alphafill" not in software_names


def _structure_method(structure: Structure) -> StructureMethod:
    has_resolution = structure.resolution > 0.0
    exp_method1 = _experimental_method(structure)
    prediction_software = _is_prediction_software(structure)

    if has_resolution:
        if exp_method1 and "x-ray" in exp_method1:
            return "X-ray"
        if exp_method1 and "electron microscopy" in exp_method1:
            return "EM"
    else:
        if exp_method1 and "nmr" in exp_method1:
            return "NMR"
        if prediction_software:
            return "Predicted"
    return "Other"


def structure_metadata(
    structure: gemmi.Structure,
    *,
    path: Path | None = None,
) -> "StructureMetadata":
    """Extract metadata from a Gemmi structure.

    If no UniProt accession is found, returns metadata with
    `uniprot_accession=None` and first chain in alphabetical order.

    If multiple accessions are found within one chain, chooses the accession
    with the highest aligned residue count from ``_struct_ref_seq``.
    Ties are resolved alphabetically by accession for deterministic behavior.

    If accessions map to multiple chains, logs a warning and returns metadata
    with `uniprot_accession=None` and first chain in alphabetical order.

    Args:
        structure: A Gemmi structure.
        path: Optional source path used only for error context.

    Returns:
        A ``StructureMetadata`` instance.

    Raises:
        ValueError: If no UniProt mapping is found in the structure.
        ChainNotFoundError: If the chain specified in the UniProt mapping is not found in the structure.
    """
    total_residue_count = nr_of_residues_in_total(structure)
    uniprot_mappings = structure_to_uniprot(structure)

    if not uniprot_mappings:
        msg = f"No UniProt mapping found in structure {structure.name}"
        logger.warning(msg)
        return _metadata_without_uniprot(
            structure,
            total_residue_count=total_residue_count,
        )

    if len(uniprot_mappings) > 1:
        uniprot_chain_pairs: list[ChainUniprotPair] = [
            ChainUniprotPair(m.chain_id, m.uniprot_accession) for m in uniprot_mappings
        ]
        logger.warning(
            "Multiple UniProt mappings found in structure %s: %s. Ignoring them all",
            structure.name,
            uniprot_chain_pairs,
        )
        return _metadata_without_uniprot(
            structure,
            total_residue_count=total_residue_count,
        )

    uniprot_mapping = next(iter(uniprot_mappings))
    auth_chain = uniprot_mapping.chain_id

    label2auth = get_label2auth_chains(structure)
    auth2label = {auth: label for label, auth in label2auth.items()}
    available_auth_chains = set(auth2label.keys())
    try:
        label_chain = auth2label[auth_chain]
    except KeyError as e:
        raise ChainNotFoundError(auth_chain, path, available_auth_chains) from e
    chain = find_chain_in_structure(structure, auth_chain)
    if chain is None:
        raise ChainNotFoundError(auth_chain, path, available_auth_chains)
    chain_length = len(chain)

    return StructureMetadata(
        id=structure.name,
        uniprot_accession=uniprot_mapping.uniprot_accession,
        resolution=structure.resolution,
        total_residue_count=total_residue_count,
        is_alphafold=_source_is_alphafold(structure),
        uniprot_start=uniprot_mapping.uniprot_start,
        uniprot_end=uniprot_mapping.uniprot_end,
        sequence_identity=uniprot_mapping.sequence_identity,
        chain_length=chain_length,
        auth_chain=auth_chain,
        label_chain=label_chain,
        method=_structure_method(structure),
    )
