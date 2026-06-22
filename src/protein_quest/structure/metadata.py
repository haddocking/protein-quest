"""Metadata extraction and ranking helpers for structures."""

import logging
from dataclasses import dataclass
from pathlib import Path

import gemmi
from gemmi import Structure

from protein_quest.structure.chains import chains_in_structure, find_chain_in_structure, nr_of_residues_in_total
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.types import StructureMethod
from protein_quest.structure.uniprot import (
    _group_struct_ref_seqs_by_chain,
    _matching_struct_ref_seqs,
    _select_best_struct_ref_seq,
    structure2uniprot_accessions,
)

logger = logging.getLogger(__name__)


def _metadata_without_uniprot(
    structure: gemmi.Structure,
    *,
    total_residue_count: int,
) -> "StructureMetadata":
    return StructureMetadata(
        id=structure.name,
        uniprot_accession=None,
        resolution=structure.resolution,
        total_residue_count=total_residue_count,
        is_alphafold=_source_is_alphafold(structure),
        uniprot_start=0,
        uniprot_end=0,
        sequence_identity=0.0,
        chain_length=total_residue_count,
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
    """Metadata extracted from a structure file for ranking and grouping."""

    id: str
    uniprot_accession: str | None
    resolution: float
    total_residue_count: int
    is_alphafold: bool
    uniprot_start: int
    uniprot_end: int
    sequence_identity: float
    chain_length: int
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
    `uniprot_accession=None`.

    If multiple accessions are found within one chain, chooses the accession
    with the highest aligned residue count from ``_struct_ref_seq``.
    Ties are resolved alphabetically by accession for deterministic behavior.

    If accessions map to multiple chains, logs a warning and returns metadata
    with `uniprot_accession=None`.

    Args:
        structure: A Gemmi structure.
        path: Optional source path used only for error context.

    Returns:
        A ``StructureMetadata`` instance.

    Raises:
        ValueError: If UniProt accessions exist but no matching
            ``_struct_ref_seq`` row is found.
        ChainNotFoundError: If the mapped chain from ``_struct_ref_seq`` is
            missing in the structure."""
    accessions = structure2uniprot_accessions(structure)
    total_residue_count = nr_of_residues_in_total(structure)

    if not accessions:
        return _metadata_without_uniprot(
            structure,
            total_residue_count=total_residue_count,
        )

    matching_struct_ref_seqs = _matching_struct_ref_seqs(structure, accessions)
    if not matching_struct_ref_seqs:
        msg = f"No struct_ref_seq entry with any of {accessions} uniprot accessions found in {structure.name}"
        raise ValueError(msg)

    struct_ref_seqs_by_chain = _group_struct_ref_seqs_by_chain(matching_struct_ref_seqs)
    if len(struct_ref_seqs_by_chain) > 1:
        logger.warning(_build_multiple_accessions_message(structure, accessions, path))
        return _metadata_without_uniprot(
            structure,
            total_residue_count=total_residue_count,
        )

    selected_struct_ref_seq = _select_best_struct_ref_seq(next(iter(struct_ref_seqs_by_chain.values())))

    uniprot_accession = selected_struct_ref_seq.uniprot_accession
    uniprot_start = selected_struct_ref_seq.uniprot_start
    uniprot_end = selected_struct_ref_seq.uniprot_end
    chain_id = selected_struct_ref_seq.chain_id
    sequence_identity = selected_struct_ref_seq.sequence_identity

    chain = find_chain_in_structure(structure, chain_id)
    if chain is None:
        raise ChainNotFoundError(chain_id, path, {c.name for c in chains_in_structure(structure)})
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
        method=_structure_method(structure),
    )
