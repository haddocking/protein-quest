"""Cluster PDB results by UniProt residue-range overlap.

Thin PDB-specific wrapper around [protein_quest.clustering][] that handles
[PdbResult][protein_quest.pdbe.result.PdbResult] entries with invalid chain
length metadata (these cannot participate in residue-range clustering and
are returned/appended separately).
"""

import logging

from protein_quest.clustering import (
    filter_structures_on_clustered_resolution,
    sort_structures,
)
from protein_quest.errors import ResolutionUnsetError
from protein_quest.pdbe.result import PdbChainLengthError, PdbResult

logging = logging.getLogger(__name__)


def _is_valid_pdb(pdb: PdbResult):
    _ = pdb.uniprot_start
    if pdb.resolution is None:
        raise ResolutionUnsetError(pdb.id)


def _separate_valid_invalid_pdbs(pdbs: list[PdbResult]) -> tuple[list[PdbResult], list[PdbResult]]:
    invalid_pdbs: list[PdbResult] = []
    valid_pdbs: list[PdbResult] = []
    for pdb in pdbs:
        try:
            _is_valid_pdb(pdb)
        except (PdbChainLengthError, ResolutionUnsetError) as e:
            logging.info(f"PDB {pdb.id} is invalid, placing last: {e}")
            invalid_pdbs.append(pdb)
        else:
            valid_pdbs.append(pdb)
    return sort_structures(invalid_pdbs), valid_pdbs


def filter_pdbs_on_clustered_resolution(pdbs: list[PdbResult], top: int) -> list[PdbResult]:
    """Filter PDB results by resolution within clusters.

    Clusters are formed by distances between PDB results based on
    Jaccard-like distance between two PDB residue ranges,
    see [structure_distance][protein_quest.clustering.structure_distance].

    The clustering is done by complete-linkage hierarchical clustering
    followed by flat clustering using a distance criterion.

    The cluster members are sorted by [sort_structures][protein_quest.clustering.sort_structures].

    The clusters are sorted by chain length, then by start and end residue, and then by PDB ID.

    The clusters are flattened to the returned list by taking the best of each cluster until clusters are exhausted.

    PDB results with invalid residue range (for example `A=-`) or unset resolution are
    always placed at the end of the returned list.

    Args:
        pdbs: PDB results to filter.
        top: Number of top results to retain per cluster.

    Returns:
        Filtered list of PDB results.
    """
    if top <= 0:
        msg = "Top must be a positive integer."
        raise ValueError(msg)

    invalid_pdbs, valid_pdbs = _separate_valid_invalid_pdbs(pdbs)
    filtered_pdbs: list[PdbResult] = []
    if valid_pdbs:
        filtered_pdbs.extend(filter_structures_on_clustered_resolution(valid_pdbs, top=top))
    if invalid_pdbs and len(filtered_pdbs) < top:
        filtered_pdbs.extend(invalid_pdbs[: top - len(filtered_pdbs)])
    return filtered_pdbs
