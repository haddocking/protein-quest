"""Module for filtering alphafold structures on confidence."""

import logging
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path

import gemmi

from protein_quest.pdbe.io import write_structure

"""
Methods to filter AlphaFoldDB structures on confidence scores.

In AlphaFold PDB files, the b-factor column has the
predicted local distance difference test (pLDDT).

See https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/
"""

logger = logging.getLogger(__name__)


def find_high_confidence_residues(structure: gemmi.Structure, confidence: float) -> Generator[int]:
    """Find residues in the structure with pLDDT confidence above the given threshold.

    Args:
        structure: The AlphaFoldDB structure to search.
        confidence: The confidence threshold (pLDDT) to use for filtering.

    Yields:
        The sequence numbers of residues with pLDDT above the confidence threshold.
    """
    for model in structure:
        for chain in model:
            for res in chain:
                res_confidence = res[0].b_iso
                if res_confidence > confidence:
                    seqid = res.seqid.num
                    if seqid is not None:
                        yield seqid


def filter_out_low_confidence_residues(structure: gemmi.Structure, allowed_residues: set[int]) -> gemmi.Structure:
    """Filter out residues from the structure that do not have high confidence.

    Args:
        structure: The AlphaFoldDB structure to filter.
        allowed_residues: The set of residue sequence numbers to keep.

    Returns:
        A new AlphaFoldDB structure with low confidence residues removed.
    """
    new_structure = structure.clone()
    for model in new_structure:
        new_chains = []
        for chain in model:
            new_chain = gemmi.Chain(chain.name)
            for res in chain:
                if res.seqid.num in allowed_residues:
                    new_chain.add_residue(res)
            new_chains.append(new_chain)
        for new_chain in new_chains:
            model.remove_chain(new_chain.name)
            model.add_chain(new_chain)
    return new_structure


@dataclass
class ConfidenceFilterQuery:
    """Query for filtering AlphaFoldDB structures based on confidence.

    Parameters:
        confidence: The confidence threshold for filtering residues.
            Residues with a pLDDT (b-factor) above this value are considered high confidence.
        min_threshold: The minimum number of high-confidence residues required to keep the structure.
        max_threshold: The maximum number of high-confidence residues required to keep the structure.
    """

    confidence: float
    min_threshold: int
    max_threshold: int


@dataclass
class ConfidenceFilterResult:
    """Result of filtering AlphaFoldDB structures based on confidence (pLDDT).

    Parameters:
        input_file: The name of the mmcif/PDB file that was processed.
        count: The number of residues with a pLDDT above the confidence threshold.
        filtered_file: The path to the filtered mmcif/PDB file, if passed filter.
    """

    input_file: str
    count: int
    filtered_file: Path | None = None


def filter_file_on_residues(file: Path, query: ConfidenceFilterQuery, filtered_dir: Path) -> ConfidenceFilterResult:
    """Filter a single AlphaFoldDB structure file (*.pdb[.gz], *.cif[.gz]) based on confidence.

    Args:
        file: The path to the PDB file to filter.
        query: The confidence filter query.
        filtered_dir: The directory to save the filtered PDB file.

    Returns:
        result with filtered_file property set to Path where filtered PDB file is saved.
            or None if structure was filtered out.
    """
    structure = gemmi.read_structure(str(file))
    residues = set(find_high_confidence_residues(structure, query.confidence))
    count = len(residues)
    if count < query.min_threshold or count > query.max_threshold:
        # Skip structure that is outside the min and max threshold
        # just return number of high confidence residues
        return ConfidenceFilterResult(
            input_file=file.name,
            count=count,
        )
    filtered_file = filtered_dir / file.name
    new_structure = filter_out_low_confidence_residues(
        structure,
        residues,
    )
    write_structure(new_structure, filtered_file)
    return ConfidenceFilterResult(
        input_file=file.name,
        count=count,
        filtered_file=filtered_file,
    )


def filter_files_on_confidence(
    alphafold_pdb_files: list[Path], query: ConfidenceFilterQuery, filtered_dir: Path
) -> Generator[ConfidenceFilterResult]:
    """Filter AlphaFoldDB structures based on confidence.

    Args:
        alphafold_pdb_files: List of mmcif/PDB files from AlphaFoldDB to filter.
        query: The confidence filter query containing the confidence thresholds.
        filtered_dir: Directory where the filtered mmcif/PDB files will be saved.

    Yields:
        For each mmcif/PDB files yields whether it was filtered or not,
            and number of residues with pLDDT above the confidence threshold.
    """
    # Note on why code looks duplicated:
    # In ../filter.py:filter_files_on_residues() we filter on number of residues on a file level
    # here we filter on file level and inside file remove low confidence residues
    for pdb_file in alphafold_pdb_files:
        yield filter_file_on_residues(pdb_file, query, filtered_dir)
