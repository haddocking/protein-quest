"""Generic clustering of UniProt-mapped structures by residue-range overlap.

Source-agnostic clustering algorithm shared by the SPARQL PDB search
([protein_quest.pdbe.clustering][]) and the 3D Beacons structure search
([protein_quest.pdbe_3dbeacons.clustering][]).

Any object satisfying the [ClusterableStructure][protein_quest.clustering.ClusterableStructure]
protocol can be clustered.
"""

import logging
from collections.abc import Iterable, Iterator
from itertools import islice
from typing import Protocol

from scipy.cluster.hierarchy import fcluster, linkage

logger = logging.getLogger(__name__)

MIN_OVERLAP_RESIDUES = 1
"""Minimum required overlap size for assigning entries to the same cluster."""
CLUSTER_DISTANCE_THRESHOLD = 1 / MIN_OVERLAP_RESIDUES
"""Maximum pairwise distance considered connected in distance-based clustering."""
NO_OVERLAP_DISTANCE = CLUSTER_DISTANCE_THRESHOLD + 1e-6
"""Finite fallback distance used when two ranges do not overlap."""


class ClusterableStructure(Protocol):
    """Protocol describing the minimum interface required for clustering.

    Attributes:
        id: Identifier of the structure, used as a deterministic tie-breaker.
        uniprot_start: Lowest UniProt residue position covered by the structure.
        uniprot_end: Highest UniProt residue position covered by the structure.
        resolution_value: Resolution in Angstrom. ``0.0`` means
            missing/undesirable so entries with a real resolution rank first.
        sequence_identity: Sequence identity of the structure to the UniProt
            sequence in range ``[0, 1]``.
            For example gaps or mutations in structure versus UniProt sequence will lower this value.
        chain_length: Number of residues in the chain mapped to the UniProt sequence.
    """

    @property
    def id(self) -> str: ...
    @property
    def uniprot_start(self) -> int: ...
    @property
    def uniprot_end(self) -> int: ...
    @property
    def resolution_value(self) -> float: ...
    @property
    def sequence_identity(self) -> float: ...
    @property
    def chain_length(self) -> int: ...


def structure_overlap(a: ClusterableStructure, b: ClusterableStructure) -> int:
    """Number of overlapping UniProt residues between two structures.

    Args:
        a: First structure.
        b: Second structure.

    Returns:
        Number of overlapping residues.
    """
    a_range = set(range(a.uniprot_start, a.uniprot_end + 1))
    b_range = set(range(b.uniprot_start, b.uniprot_end + 1))
    return len(a_range.intersection(b_range))


def structure_union(a: ClusterableStructure, b: ClusterableStructure) -> int:
    """Number of unique UniProt residues in the union of two structures.

    Args:
        a: First structure.
        b: Second structure.

    Returns:
        Number of unique residues in the union of the two structures.
    """
    a_range = set(range(a.uniprot_start, a.uniprot_end + 1))
    b_range = set(range(b.uniprot_start, b.uniprot_end + 1))
    return len(a_range.union(b_range))


def structure_distance(a: ClusterableStructure, b: ClusterableStructure) -> float:
    """Jaccard-like distance between two structures' UniProt residue ranges.

    Non-overlapping ranges return [NO_OVERLAP_DISTANCE][protein_quest.clustering.NO_OVERLAP_DISTANCE].

    Args:
        a: First structure.
        b: Second structure.

    Returns:
        Distance between the two structures in range ``[0, NO_OVERLAP_DISTANCE]``
    """
    overlap = structure_overlap(a, b)
    if overlap == 0:
        return NO_OVERLAP_DISTANCE
    union = structure_union(a, b)
    return 1 - (overlap / union)


def structure_distances[T: ClusterableStructure](items: list[T]) -> list[float]:
    """Condensed (upper-triangle) pairwise distance matrix for a list of structures.

    Args:
        items: Structures to compute distances for.

    Returns:
        Condensed distance matrix as a flat list, suitable for input to
            [scipy.cluster.hierarchy.linkage][scipy.cluster.hierarchy.linkage].
    """
    condensed: list[float] = []
    for index, left in enumerate(items[:-1]):
        condensed.extend(structure_distance(left, right) for right in items[index + 1 :])
    return condensed


def _cluster_sort_key[T: ClusterableStructure](cluster: set[T]) -> tuple[int, int, int, str]:
    max_chain_length = max(member.chain_length for member in cluster)
    start, end, ident = min((member.uniprot_start, member.uniprot_end, member.id) for member in cluster)
    return (-max_chain_length, start, end, ident)


def _member_sort_key(member: ClusterableStructure) -> tuple[float, float, int, str]:
    """Deterministic quality sort key for a cluster member.

    A failing ``chain_length`` access (for example for PDB results with
    unparsable chain metadata) is treated as ``0`` so such entries can still
    be sorted alongside valid ones.

    See [sort_structures][protein_quest.clustering.sort_structures] for criteria.
    """
    try:
        chain_length = member.chain_length
    except Exception:  # noqa: BLE001 - sort-key fallback for adapters that derive chain_length lazily
        chain_length = 0
    return (
        -member.sequence_identity,
        member.resolution_value,
        -chain_length,
        member.id,
    )


def sort_structures[T: ClusterableStructure](items: set[T] | list[T]) -> list[T]:
    """Sort structures by quality criteria.

    1. Sequence identity descending (highest first)
    2. Resolution ascending (lowest first; ``0.0`` is treated as missing and ranks first)
    3. Chain length descending (longest first)
    4. Identifier ascending (deterministic tie-break)

    Args:
        items: Structures to sort.

    Returns:
        List of structures sorted by the criteria above.
    """
    return sorted(items, key=_member_sort_key)


def _cluster_many_valid[T: ClusterableStructure](items: list[T]) -> list[list[T]]:
    condensed_distances = structure_distances(items)
    linkage_matrix = linkage(condensed_distances, method="complete")
    cluster_ids = fcluster(linkage_matrix, t=CLUSTER_DISTANCE_THRESHOLD, criterion="distance")

    logger.info("Formed %d clusters from %d structures", len(set(cluster_ids)), len(items))

    clusters_by_id: dict[int, set[T]] = {}
    for item, cluster_id in zip(items, cluster_ids, strict=True):
        clusters_by_id.setdefault(int(cluster_id), set()).add(item)

    return [
        sort_structures(clusters_by_id[cluster_id])
        for cluster_id in sorted(clusters_by_id, key=lambda cid: _cluster_sort_key(clusters_by_id[cid]))
    ]


def cluster_structures[T: ClusterableStructure](items: list[T]) -> list[list[T]]:
    """Cluster structures by overlapping UniProt residue coverage.

    Clustering is done using [linkage(distances, method="complete")][scipy.cluster.hierarchy.linkage] and
    [fcluster(linkage_matrix, t=CLUSTER_DISTANCE_THRESHOLD, criterion="distance")][scipy.cluster.hierarchy.fcluster].

    Args:
        items: Structures to cluster. All items must have valid residue ranges;
            callers responsible for filtering out structures with missing
            range information beforehand.

    Returns:
        Ordered list of clusters; each cluster is a list of members sorted
            by [sort_structures][protein_quest.clustering.sort_structures].
    """
    ordered = sorted(items, key=lambda item: (item.uniprot_start, item.uniprot_end, item.id))
    if not ordered:
        return []
    if len(ordered) == 1:
        return [[ordered[0]]]
    if len(ordered) == 2:
        a, b = ordered
        if structure_distance(a, b) <= CLUSTER_DISTANCE_THRESHOLD:
            return [[a, b]]
        return [[a], [b]]
    return _cluster_many_valid(ordered)


def interleave_longest[T](*iterables: Iterable[T]) -> Iterator[T]:
    """Yield values round-robin from each iterable until all are exhausted.

    Examples:
        Scalar example.

        >>> list(interleave_longest([1, 2, 3], [4, 5], [6, 7, 8]))
        [1, 4, 6, 2, 5, 7, 3, 8]

    Args:
        iterables: Iterables to round-robin.

    Yields:
        Elements interleaved from the provided iterables.
    """
    active = [iter(iterable) for iterable in iterables]
    while active:
        next_round: list[Iterator[T]] = []
        for iterator in active:
            try:
                yield next(iterator)
            except StopIteration:
                continue
            next_round.append(iterator)
        active = next_round


def top_members_of_clusters[T](clusters: list[list[T]], top: int) -> list[T]:
    """Return up to ``top`` members by interleaving cluster members round-robin.

    Args:
        clusters: Ordered clusters whose members are also sorted.
            First cluster and first member of each cluster are considered best.
        top: Maximum number of members to return.

    Returns:
        Interleaved members truncated to ``top``.
    """
    return list(islice(interleave_longest(*clusters), top))


def filter_structures_on_clustered_resolution[T: ClusterableStructure](items: list[T], top: int) -> list[T]:
    """Filter structures by resolution within residue-range clusters.

    Interleaves cluster members round-robin (per
    [interleave_longest][protein_quest.clustering.interleave_longest]),
    taking the best remaining member of each cluster (per
    [sort_structures][protein_quest.clustering.sort_structures]) until
    ``top`` entries have been collected or all clusters are exhausted.

    Clusters themselves are ordered by chain length descending, then by start
    and end residue, then identifier.

    Args:
        items: Structures to filter. Must all have valid residue ranges.
        top: Number of top results to retain.

    Returns:
        Filtered list of up to ``top`` structures.
    """
    if top <= 0:
        msg = "Top must be a positive integer."
        raise ValueError(msg)

    clusters = cluster_structures(items)
    return top_members_of_clusters(clusters, top)
