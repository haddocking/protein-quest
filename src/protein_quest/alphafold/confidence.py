"""Module for filtering alphafold structures on confidence."""

import logging
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import gemmi
from dask.distributed import Client
from distributed.deploy.cluster import Cluster
from tqdm.auto import tqdm

from protein_quest.converter import Percentage, PositiveInt, converter
from protein_quest.io import read_structure, write_structure
from protein_quest.parallel import configure_dask_scheduler, dask_map_with_progress
from protein_quest.ss import nr_of_residues_in_total
from protein_quest.utils import CopyMethod, copyfile

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
        min_residues: The minimum number of high-confidence residues required to keep the structure.
        max_residues: The maximum number of high-confidence residues required to keep the structure.
    """

    confidence: Percentage
    min_residues: PositiveInt
    max_residues: PositiveInt


base_query_hook = converter.get_structure_hook(ConfidenceFilterQuery)


@converter.register_structure_hook
def confidence_filter_query_hook(val, _type) -> ConfidenceFilterQuery:
    result: ConfidenceFilterQuery = base_query_hook(val, _type)
    if result.min_residues > result.max_residues:
        msg = f"min_residues {result.min_residues} cannot be larger than max_residues {result.max_residues}"
        raise ValueError(msg)
    return result


@dataclass
class ConfidenceFilterResult:
    """Result of filtering AlphaFoldDB structures based on confidence (pLDDT).

    Parameters:
        input_file: The name of the mmcif/PDB file that was processed.
        count: The number of residues with a pLDDT above the confidence threshold.
        filtered_file: The path to the filtered mmcif/PDB file, if passed filter.
    """

    input_file: str
    count: PositiveInt
    filtered_file: Path | None = None


def filter_file_on_confidence(
    file: Path, query: ConfidenceFilterQuery, filtered_dir: Path, copy_method: CopyMethod = "copy"
) -> ConfidenceFilterResult:
    """Filter a single AlphaFoldDB structure file (*.pdb[.gz], *.cif[.gz]) based on confidence.

    Args:
        file: The path to the PDB file to filter.
        query: The confidence filter query.
        filtered_dir: The directory to save the filtered PDB file.
        copy_method: How to copy when no residues have to be removed.

    Returns:
        result with filtered_file property set to Path where filtered PDB file is saved.
            or None if structure was filtered out.
    """
    structure = read_structure(file)
    residues = set(find_high_confidence_residues(structure, query.confidence))
    count = len(residues)
    if count < query.min_residues or count > query.max_residues:
        # Skip structure that is outside the min and max threshold
        # just return number of high confidence residues
        return ConfidenceFilterResult(
            input_file=file.name,
            count=count,
        )
    total_residues = nr_of_residues_in_total(structure)
    filtered_file = filtered_dir / file.name
    if count == total_residues:
        # if no residues have to be removed then copy instead of slower gemmi writing
        copyfile(file, filtered_file, copy_method)
    else:
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


def _filter_files_on_confidence_sequentially(
    alphafold_pdb_files: list[Path],
    query: ConfidenceFilterQuery,
    filtered_dir: Path,
    copy_method: CopyMethod = "copy",
) -> list[ConfidenceFilterResult]:
    results = []
    for file in tqdm(
        alphafold_pdb_files,
        total=len(alphafold_pdb_files),
        desc="Filtering on confidence",
        unit="file",
    ):
        result = filter_file_on_confidence(file, query, filtered_dir, copy_method)
        results.append(result)
    return results


def filter_files_on_confidence(
    alphafold_pdb_files: list[Path],
    query: ConfidenceFilterQuery,
    filtered_dir: Path,
    copy_method: CopyMethod = "copy",
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
) -> list[ConfidenceFilterResult]:
    """Filter AlphaFoldDB structures based on confidence.

    Args:
        alphafold_pdb_files: List of mmcif/PDB files from AlphaFoldDB to filter.
        query: The confidence filter query containing the confidence thresholds.
        filtered_dir: Directory where the filtered mmcif/PDB files will be saved.
        copy_method: How to copy when a direct copy is possible.
        scheduler_address: The address of the Dask scheduler.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.

    Returns:
        For each mmcif/PDB files returns whether it was filtered or not,
            and number of residues with pLDDT above the confidence threshold.
    """
    filtered_dir.mkdir(parents=True, exist_ok=True)
    if scheduler_address == "sequential":
        return _filter_files_on_confidence_sequentially(
            alphafold_pdb_files,
            query,
            filtered_dir,
            copy_method=copy_method,
        )

    with configure_dask_scheduler(scheduler_address, name="filter-confidence") as cluster, Client(cluster) as client:
        client.forward_logging()
        return dask_map_with_progress(
            client,
            filter_file_on_confidence,
            alphafold_pdb_files,
            query=query,
            filtered_dir=filtered_dir,
            copy_method=copy_method,
        )
