"""Cluster PDB results by UniProt residue-range overlap.

Attributes:
    CLUSTER_DISTANCE_THRESHOLD: Maximum pairwise distance considered connected in distance-based clustering.
    NO_OVERLAP_DISTANCE: Finite fallback distance used when two ranges do not overlap.
    MIN_OVERLAP_RESIDUES: Minimum required overlap size for assigning entries to the same cluster.
"""

import logging

from scipy.cluster.hierarchy import fcluster, linkage

from protein_quest.pdbe.result import PdbChainLengthError, PdbResult

logger = logging.getLogger(__name__)

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


def pdb_union(a: PdbResult, b: PdbResult) -> int:
    """Compute the number of unique UniProt residues in the union of two PDB residue ranges.

    Args:
        a: First PDB result.
        b: Second PDB result.

    Returns:
        Number of unique UniProt residues in the union of the two residue ranges.
    """
    a_range = set(range(a.uniprot_start, a.uniprot_end + 1))
    b_range = set(range(b.uniprot_start, b.uniprot_end + 1))
    return len(a_range.union(b_range))


def pdb_distance(a: PdbResult, b: PdbResult) -> float:
    """Compute Jaccard-like distance between two PDB residue ranges.

    The distance is defined as 1 - (overlap / union), which is the Jaccard distance.

    Args:
        a: First PDB result.
        b: Second PDB result.

    Returns:
        Jaccard distance between the two residue ranges (0.0 for identical ranges).
        Non-overlapping ranges return NO_OVERLAP_DISTANCE sentinel value.
    """
    overlap = pdb_overlap(a, b)
    if overlap == 0:
        return NO_OVERLAP_DISTANCE

    union = pdb_union(a, b)
    return 1 - (overlap / union)


def _cluster_sort_key(cluster: set[PdbResult]) -> tuple[int, int, int, str]:
    # We know members have valid chain_length, so not checking for PdbChainLengthError here.
    max_chain_length = max(member.chain_length for member in cluster)
    start, end, pdb_id = min((member.uniprot_start, member.uniprot_end, member.id) for member in cluster)
    return (-max_chain_length, start, end, pdb_id)


def _cluster_member_sort_key(member: PdbResult) -> tuple[float, str | None, int, str]:
    """Return deterministic quality sort key for a cluster member.

    See [sort_pdbs][protein_quest.pdbe.clustering.sort_pdbs] for sort criteria.
    """
    try:
        chain_length = member.chain_length
    except PdbChainLengthError:
        chain_length = 0

    return (
        -member.sequence_identity,
        member.resolution,
        -chain_length,
        member.id,
    )


def sort_pdbs(pdbs: set[PdbResult]) -> list[PdbResult]:
    """Sort pdbs by quality criteria.

    1. Sequence identity descending (highest first)
    2. Resolution ascending (lowest first)
    3. Chain length descending (longest first)
    4. PDB ID ascending (deterministic tie-break)

    Args:
        pdbs: PDB results to sort.

    Returns:
        PDB results sorted by quality criteria.
    """
    return sorted(
        pdbs,
        key=_cluster_member_sort_key,
    )


def pdb_distances(pdbs: list[PdbResult]) -> list[float]:
    """Compute pairwise distances between PDB results.

    Args:
        pdbs: PDB results to compute distances for.

    Returns:
        Condensed distance matrix (upper triangle) of pairwise distances.
    """
    condensed_distances: list[float] = []
    for index, left in enumerate(pdbs[:-1]):
        condensed_distances.extend(pdb_distance(left, right) for right in pdbs[index + 1 :])
    return condensed_distances


def _cluster_many_valid_pdbs(valid_pdbs: list[PdbResult]) -> list[list[PdbResult]]:
    condensed_distances = pdb_distances(valid_pdbs)
    linkage_matrix = linkage(condensed_distances, method="complete")
    cluster_ids = fcluster(linkage_matrix, t=CLUSTER_DISTANCE_THRESHOLD, criterion="distance")

    logger.info("Formed %d clusters from %d PDBs with valid range(s)", len(set(cluster_ids)), len(valid_pdbs))

    clusters_by_id: dict[int, set[PdbResult]] = {}
    for pdb, cluster_id in zip(valid_pdbs, cluster_ids, strict=True):
        clusters_by_id.setdefault(int(cluster_id), set()).add(pdb)

    return [
        sort_pdbs(clusters_by_id[cluster_id])
        for cluster_id in sorted(clusters_by_id, key=lambda cluster_id: _cluster_sort_key(clusters_by_id[cluster_id]))
    ]


def _cluster_valid_pdbs(valid_pdbs: list[PdbResult]) -> list[list[PdbResult]]:
    ordered_valid_pdbs = sorted(valid_pdbs, key=lambda pdb: (pdb.uniprot_start, pdb.uniprot_end, pdb.id))
    if not ordered_valid_pdbs:
        return []
    if len(ordered_valid_pdbs) == 1:
        return [[ordered_valid_pdbs[0]]]
    if len(ordered_valid_pdbs) == 2:
        a, b = ordered_valid_pdbs
        if pdb_distance(a, b) <= CLUSTER_DISTANCE_THRESHOLD:
            return [[a, b]]
        return [[a], [b]]
    return _cluster_many_valid_pdbs(ordered_valid_pdbs)


def _separate_valid_invalid_pdbs(pdbs: list[PdbResult]) -> tuple[list[PdbResult], list[PdbResult]]:
    invalid_pdbs: set[PdbResult] = set()
    valid_pdbs: list[PdbResult] = []
    for pdb in pdbs:
        try:
            _ = pdb.uniprot_start
        except PdbChainLengthError:
            invalid_pdbs.add(pdb)
        else:
            valid_pdbs.append(pdb)
    return sort_pdbs(invalid_pdbs), valid_pdbs


def cluster_pdbs(pdbs: list[PdbResult]) -> tuple[list[list[PdbResult]], list[PdbResult]]:
    """Cluster PDB results by overlapping UniProt residue coverage.

    Results with valid chain length metadata are clustered by overlap. Results
    with invalid chain length metadata are returned separately.

    Args:
        pdbs: PDB results to cluster.

    Returns:
        Tuple of ordered valid clusters and invalid PDB results.
    """
    if not pdbs:
        return [], []

    invalid_pdbs, valid_pdbs = _separate_valid_invalid_pdbs(pdbs)
    valid_clusters = _cluster_valid_pdbs(valid_pdbs)
    return valid_clusters, invalid_pdbs


def filter_pdbs_on_clustered_resolution(pdbs: list[PdbResult], top: int) -> list[PdbResult]:
    """Filter PDB results by resolution within clusters.

    Clusters are formed by distances between PDB results based on
    Jaccard-like distance between two PDB residue ranges,
    see [pdb_distance][protein_quest.pdbe.clustering.pdb_distance].

    The clustering is done by complete-linkage hierarchical clustering
    followed by flat clustering using a distance criterion.

    The cluster members are sorted by [sort_pdbs][protein_quest.pdbe.clustering.sort_pdbs].

    The clusters are sorted by chain length, then by start and end residue, and then by PDB ID.

    The clusters are flattend to the returned list by taking the best of each cluster until clusters are exhausted.

    PDB results with invalid residue range (for example `A=-`) are always placed at the end of the returned list.

    Args:
        pdbs: PDB results to filter.
        top: Number of top results to retain per cluster.

    Returns:
        Filtered list of PDB results.
    """
    if top <= 0:
        msg = "Top must be a positive integer."
        raise ValueError(msg)

    valid_clusters, invalid_cluster = cluster_pdbs(pdbs)
    filtered_pdbs: list[PdbResult] = []
    cluster_id = 0
    while len(filtered_pdbs) < top and any(valid_clusters):
        members = valid_clusters[cluster_id]
        if members:
            filtered_pdbs.append(members.pop(0))
        cluster_id += 1
        if cluster_id >= len(valid_clusters):
            cluster_id = 0
    if invalid_cluster and len(filtered_pdbs) < top:
        filtered_pdbs.extend(invalid_cluster[: top - len(filtered_pdbs)])

    return filtered_pdbs
