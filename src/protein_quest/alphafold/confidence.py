import logging
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path

"""
Methods to filter AlphaFoldDB structures on confidence scores.

In AlphaFold PDB files, the b-factor column has the
predicted local distance difference test (pLDDT).

See https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/
"""

logger = logging.getLogger(__name__)


def find_high_confidence_residues(pdb_file: Path, confidence: float) -> Generator[int]:
    prev_res_index = None
    with pdb_file.open("r") as f:
        for line in f:
            if line.startswith("ATOM"):
                # Extract the b-factor value from the PDB line
                b_factor = float(line[60:66].strip())
                if b_factor > confidence:
                    # Extract the residue index from the PDB line
                    res_index = int(line[22:26].strip())
                    # yield once per residue, not for every atom of the residue
                    if res_index != prev_res_index:
                        yield res_index
                        prev_res_index = res_index


def filter_out_low_confidence_residues(input_pdb_file: Path, allowed_residues: set[int], output_pdb_file: Path):
    # TODO if residue is removed from ATOM lines also remove it
    # elsewhere in the file (SEQRES, TER).
    # TODO do we need to take model/chain into account?
    # now assumes single model and single chain
    if output_pdb_file.exists():
        logger.info(f"Output file {output_pdb_file} already exists. Skipping filtering for {input_pdb_file}.")
        return
    with input_pdb_file.open("r") as input_file, output_pdb_file.open("w") as output_file:
        for line in input_file:
            if line.startswith("ATOM"):
                # Extract the residue index from the PDB line
                res_index = int(line[22:26].strip())
                if res_index in allowed_residues:
                    output_file.write(line)
            else:
                output_file.write(line)


@dataclass
class ConfidenceFilterQuery:
    """Query for filtering AlphaFoldDB structures based on confidence.

    Parameters:
        confidence: The confidence threshold for filtering residues.
            Residues with a pLDDT (b-factor) above this value are considered high confidence.
        min_threshold: The minimum number of high-confidence residues required to keep the structure.
        max_threshold: The maximum number of high-confidence residues allowed to keep the structure.
    """

    confidence: float
    min_threshold: int
    max_threshold: int


@dataclass
class ConfidenceFilterResult:
    """Result of filtering AlphaFoldDB structures based on confidence (pLDDT).

    Parameters:
        pdb_file: The name of the PDB file that was processed.
        count: The number of residues with a pLDDT above the confidence threshold.
        filtered_file: The path to the filtered PDB file, if passed filter.
    """

    pdb_file: str
    count: int
    filtered_file: Path | None = None


# TODO make work on cif file, not just pdb formatted files
def filter_files_on_confidence(
    alphafold_pdb_files: list[Path], query: ConfidenceFilterQuery, filtered_dir: Path
) -> Generator[ConfidenceFilterResult]:
    """Filter AlphaFoldDB structures based on confidence.

    Args:
        alphafold_pdb_files: List of PDB files from AlphaFoldDB to filter.
        query: The confidence filter query containing the confidence thresholds.
        filtered_dir: Directory where the filtered PDB files will be saved.

    Yields:
        For each PDB files yields whether it was filtered or not,
            and number of residues with pLDDT above the confidence threshold.
    """
    for pdb_file in alphafold_pdb_files:
        residues = set(find_high_confidence_residues(pdb_file, query.confidence))
        count = len(residues)
        if count < query.min_threshold or count > query.max_threshold:
            yield ConfidenceFilterResult(
                pdb_file=pdb_file.name,
                count=count,
            )
            # Skip structure that is outside the min and max threshold
            continue
        filtered_file = filtered_dir / pdb_file.name
        filter_out_low_confidence_residues(
            pdb_file,
            residues,
            filtered_file,
        )
        yield ConfidenceFilterResult(
            pdb_file=pdb_file.name,
            count=count,
            filtered_file=filtered_file,
        )
