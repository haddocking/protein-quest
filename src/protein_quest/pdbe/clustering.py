"""Cluster PDB results by UniProt residue-range overlap.

Attributes:
    NO_OVERLAP_DISTANCE: Finite fallback distance used when two ranges do not overlap.
    MIN_OVERLAP_RESIDUES: Minimum required overlap size for assigning entries to the same cluster.
"""

from math import isinf

from scipy.cluster.hierarchy import fcluster, linkage

from protein_quest.pdbe.result import PdbChainLengthError, PdbResult

NO_OVERLAP_DISTANCE = 500_000
MIN_OVERLAP_RESIDUES = 1


def pdb_distance(a: PdbResult, b: PdbResult) -> float:
    """Compute inverse overlap distance between two PDB residue ranges.

    The distance is defined as the reciprocal of the number of overlapping
    UniProt residues. Non-overlapping ranges return positive infinity.

    Args:
        a: First PDB result.
        b: Second PDB result.

    Returns:
        Inverse overlap distance between the two residue ranges.
    """
    a_range = set(range(a.uniprot_start, a.uniprot_end + 1))
    b_range = set(range(b.uniprot_start, b.uniprot_end + 1))
    overlap = a_range.intersection(b_range)
    if not overlap:
        return float("+inf")
    return 1 / len(overlap)


def _cluster_sort_key(cluster: set[PdbResult]) -> tuple[int, int, str]:
    return min((member.uniprot_start, member.uniprot_end, member.id) for member in cluster)


def _cluster_distance(a: PdbResult, b: PdbResult) -> float:
    distance = pdb_distance(a, b)
    if isinf(distance):
        return NO_OVERLAP_DISTANCE
    return distance


def _cluster_many_valid_pdbs(valid_pdbs: list[PdbResult]) -> list[set[PdbResult]]:
    condensed_distances: list[float] = []
    for index, left in enumerate(valid_pdbs[:-1]):
        condensed_distances.extend(_cluster_distance(left, right) for right in valid_pdbs[index + 1 :])

    linkage_matrix = linkage(condensed_distances, method="complete")
    cluster_labels = fcluster(linkage_matrix, t=1 / MIN_OVERLAP_RESIDUES, criterion="distance")
    clusters_by_label: dict[int, set[PdbResult]] = {}
    for pdb, label in zip(valid_pdbs, cluster_labels, strict=True):
        clusters_by_label.setdefault(int(label), set()).add(pdb)

    return [
        clusters_by_label[label]
        for label in sorted(clusters_by_label, key=lambda label: _cluster_sort_key(clusters_by_label[label]))
    ]


def _cluster_valid_pdbs(valid_pdbs: list[PdbResult]) -> list[set[PdbResult]]:
    ordered_valid_pdbs = sorted(valid_pdbs, key=lambda pdb: (pdb.uniprot_start, pdb.uniprot_end, pdb.id))
    if not ordered_valid_pdbs:
        return []
    if len(ordered_valid_pdbs) == 1:
        return [{ordered_valid_pdbs[0]}]
    if len(ordered_valid_pdbs) == 2:
        a, b = ordered_valid_pdbs
        if _cluster_distance(a, b) <= 1 / MIN_OVERLAP_RESIDUES:
            return [{a, b}]
        return [{a}, {b}]
    return _cluster_many_valid_pdbs(ordered_valid_pdbs)


def cluster_pdbs(pdbs: list[PdbResult]) -> list[set[PdbResult]]:
    """Cluster PDB results by overlapping UniProt residue coverage.

    Results with valid chain length metadata are clustered by overlap. Results
    with invalid chain length metadata are grouped together into an additional
    cluster appended to the end of the result.

    Args:
        pdbs: PDB results to cluster.

    Returns:
        Ordered clusters of PDB results.
    """
    if not pdbs:
        return []

    invalid_pdbs: set[PdbResult] = set()
    valid_pdbs: list[PdbResult] = []
    for pdb in pdbs:
        try:
            _ = (pdb.uniprot_start, pdb.uniprot_end)
        except PdbChainLengthError:
            invalid_pdbs.add(pdb)
        else:
            valid_pdbs.append(pdb)

    clusters = _cluster_valid_pdbs(valid_pdbs)
    if invalid_pdbs:
        clusters.append(invalid_pdbs)
    return clusters
