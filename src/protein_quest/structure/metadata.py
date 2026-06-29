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
    selected_struct_ref_seqs_by_chain,
    structure2uniprot_accessions,
)

logger = logging.getLogger(__name__)


def _metadata_without_uniprot(
    structure: gemmi.Structure,
    *,
    total_residue_count: int,
    auth_chain: str,
    label_chain: str,
) -> "StructureMetadata":
    chain = find_chain_in_structure(structure, auth_chain)
    if chain is None:
        raise ChainNotFoundError(auth_chain, None, set(auth_chain))
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


def _build_multiple_accessions_message(structure: gemmi.Structure, accessions: set[str], path: Path | None) -> str:
    msg = (
        f"Multiple UniProt accessions found in structure {structure.name}: "
        f"{accessions}. Please resolve this ambiguity before using this "
        "structure. For example using `protein-quest filter chain` command."
    )
    if path is not None:
        msg = f"{msg} Source path: {path}."
    return msg


@dataclass(frozen=True)
class StructureMetadata:
    """Metadata extracted from a structure file for ranking and grouping.

    Attributes:
        id: The structure ID.
        uniprot_accession: The UniProt accession if available, otherwise None.
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
        ValueError: If UniProt accessions exist but no matching
            ``_struct_ref_seq`` row is found.
            or if no chains are found in the structure.
        ChainNotFoundError: If the mapped chain from ``_struct_ref_seq`` is
            missing in the structure."""
    accessions = structure2uniprot_accessions(structure)
    total_residue_count = nr_of_residues_in_total(structure)
    label2auth = get_label2auth_chains(structure)
    auth2label = {auth: label for label, auth in label2auth.items()}
    available_auth_chains = set(auth2label.keys())

    if not available_auth_chains:
        msg = f"No chains found in structure {structure.name}"
        raise ValueError(msg)

    if not accessions:
        first_auth_chain = sorted(available_auth_chains)[0]
        return _metadata_without_uniprot(
            structure,
            total_residue_count=total_residue_count,
            auth_chain=first_auth_chain,
            label_chain=auth2label[first_auth_chain],
        )

    struct_ref_seqs_by_chain = selected_struct_ref_seqs_by_chain(structure, accessions)
    if not struct_ref_seqs_by_chain:
        msg = f"No struct_ref_seq entry with any of {accessions} uniprot accessions found in {structure.name}"
        raise ValueError(msg)

    if len(struct_ref_seqs_by_chain) > 1:
        logger.warning(_build_multiple_accessions_message(structure, accessions, path))
        first_auth_chain = sorted(available_auth_chains)[0]
        return _metadata_without_uniprot(
            structure,
            total_residue_count=total_residue_count,
            auth_chain=first_auth_chain,
            label_chain=auth2label[first_auth_chain],
        )

    selected_struct_ref_seq = next(iter(struct_ref_seqs_by_chain.values()))

    uniprot_accession = selected_struct_ref_seq.uniprot_accession
    uniprot_start = selected_struct_ref_seq.uniprot_start
    uniprot_end = selected_struct_ref_seq.uniprot_end
    auth_chain = selected_struct_ref_seq.chain_id
    try:
        label_chain = auth2label[auth_chain]
    except KeyError as e:
        raise ChainNotFoundError(auth_chain, path, available_auth_chains) from e
    sequence_identity = selected_struct_ref_seq.sequence_identity

    chain = find_chain_in_structure(structure, auth_chain)
    if chain is None:
        raise ChainNotFoundError(auth_chain, path, available_auth_chains)
    chain_length = len(chain)

    return StructureMetadata(
        id=structure.name,
        uniprot_accession=uniprot_accession,
        resolution=structure.resolution,
        total_residue_count=total_residue_count,
        is_alphafold=_source_is_alphafold(structure),
        uniprot_start=uniprot_start,
        uniprot_end=uniprot_end,
        sequence_identity=sequence_identity,
        chain_length=chain_length,
        auth_chain=auth_chain,
        label_chain=label_chain,
        method=_structure_method(structure),
    )
