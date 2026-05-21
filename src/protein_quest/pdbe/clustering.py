"""Cluster PDB results by UniProt residue-range overlap.

Attributes:
    CLUSTER_DISTANCE_THRESHOLD: Maximum pairwise distance considered connected in distance-based clustering.
    NO_OVERLAP_DISTANCE: Finite fallback distance used when two ranges do not overlap.
    MIN_OVERLAP_RESIDUES: Minimum required overlap size for assigning entries to the same cluster.
"""

from scipy.cluster.hierarchy import fcluster, linkage

from protein_quest.pdbe.result import PdbChainLengthError, PdbResult

MIN_OVERLAP_RESIDUES = 1
CLUSTER_DISTANCE_THRESHOLD = 1 / MIN_OVERLAP_RESIDUES
NO_OVERLAP_DISTANCE = CLUSTER_DISTANCE_THRESHOLD + 1e-6


def pdb_overlap(a: PdbResult, b: PdbResult) -> int:
    """Compute the number of overlapping UniProt residues between two PDB residue ranges.

    Args:
        a: First PDB result.
        b: Second PDB result.

    Returns:
        Number of overlapping UniProt residues between the two residue ranges.
    """
    a_range = set(range(a.uniprot_start, a.uniprot_end + 1))
    b_range = set(range(b.uniprot_start, b.uniprot_end + 1))
    return len(a_range.intersection(b_range))


def pdb_distance(a: PdbResult, b: PdbResult) -> float:
    """Compute inverse overlap distance between two PDB residue ranges.

    The distance is defined as the reciprocal of the number of overlapping
    UniProt residues.

    Args:
        a: First PDB result.
        b: Second PDB result.

    Returns:
        Inverse overlap distance between the two residue ranges.
        Non-overlapping ranges return NO_OVERLAP_DISTANCE sentinel value.
    """
    overlap = pdb_overlap(a, b)
    if overlap == 0:
        return NO_OVERLAP_DISTANCE

    return 1 / overlap


def _cluster_sort_key(cluster: set[PdbResult]) -> tuple[int, int, str]:
    return min((member.uniprot_start, member.uniprot_end, member.id) for member in cluster)


def _cluster_many_valid_pdbs(valid_pdbs: list[PdbResult]) -> list[set[PdbResult]]:
    condensed_distances: list[float] = []
    for index, left in enumerate(valid_pdbs[:-1]):
        condensed_distances.extend(pdb_distance(left, right) for right in valid_pdbs[index + 1 :])

    linkage_matrix = linkage(condensed_distances, method="complete")
    cluster_ids = fcluster(linkage_matrix, t=CLUSTER_DISTANCE_THRESHOLD, criterion="distance")
    clusters_by_id: dict[int, set[PdbResult]] = {}
    for pdb, cluster_id in zip(valid_pdbs, cluster_ids, strict=True):
        clusters_by_id.setdefault(int(cluster_id), set()).add(pdb)

    return [
        clusters_by_id[cluster_id]
        for cluster_id in sorted(clusters_by_id, key=lambda cluster_id: _cluster_sort_key(clusters_by_id[cluster_id]))
    ]


def _cluster_valid_pdbs(valid_pdbs: list[PdbResult]) -> list[set[PdbResult]]:
    ordered_valid_pdbs = sorted(valid_pdbs, key=lambda pdb: (pdb.uniprot_start, pdb.uniprot_end, pdb.id))
    if not ordered_valid_pdbs:
        return []
    if len(ordered_valid_pdbs) == 1:
        return [{ordered_valid_pdbs[0]}]
    if len(ordered_valid_pdbs) == 2:
        a, b = ordered_valid_pdbs
        if pdb_distance(a, b) <= CLUSTER_DISTANCE_THRESHOLD:
            return [{a, b}]
        return [{a}, {b}]
    return _cluster_many_valid_pdbs(ordered_valid_pdbs)


def _separate_valid_invalid_pdbs(pdbs: list[PdbResult]) -> tuple[set[PdbResult], list[PdbResult]]:
    invalid_pdbs: set[PdbResult] = set()
    valid_pdbs: list[PdbResult] = []
    for pdb in pdbs:
        try:
            _ = pdb.uniprot_start
        except PdbChainLengthError:
            invalid_pdbs.add(pdb)
        else:
            valid_pdbs.append(pdb)
    return invalid_pdbs, valid_pdbs


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

    invalid_pdbs, valid_pdbs = _separate_valid_invalid_pdbs(pdbs)
    clusters = _cluster_valid_pdbs(valid_pdbs)
    if invalid_pdbs:
        clusters.append(invalid_pdbs)
    return clusters
