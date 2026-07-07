"""Generic clustering of UniProt-mapped structures by residue-range overlap.

Any object satisfying the [ClusterableStructure][protein_quest.clustering.ClusterableStructure]
protocol can be clustered.
"""

import logging
from collections.abc import Hashable, Iterable, Iterator
from itertools import islice
from typing import Protocol

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

logger = logging.getLogger(__name__)

MIN_OVERLAP_RESIDUES = 1
"""Minimum required overlap size for assigning entries to the same cluster."""
CLUSTER_DISTANCE_THRESHOLD = 1 / MIN_OVERLAP_RESIDUES
"""Maximum pairwise distance considered connected in distance-based clustering."""
NO_OVERLAP_DISTANCE = CLUSTER_DISTANCE_THRESHOLD + 1e-6
"""Finite fallback distance used when two ranges do not overlap."""


class SortableStructure(Hashable, Protocol):
    """Protocol describing the minimum interface required for sorting.

    Attributes:
        id: Identifier of the structure, used as a deterministic tie-breaker.
        resolution_value: Resolution in Angstrom. ``0.0`` means
            missing/undesirable so entries with a real resolution rank first.
        sequence_identity: Sequence identity of the structure to the UniProt
            sequence in range ``[0, 1]``.
            For example gaps or mutations in structure versus UniProt sequence will lower this value.
        chain_length: Number of residues in the chain mapped to the UniProt sequence.
        geometry_quality: Geometry quality score (``0.0`` - ``100.0``) or
            ``None`` if unavailable. Higher is better.
    """

    @property
    def id(self) -> str: ...
    @property
    def resolution_value(self) -> float: ...
    @property
    def sequence_identity(self) -> float: ...
    @property
    def chain_length(self) -> int: ...
    @property
    def geometry_quality(self) -> float | None: ...


class ClusterableStructure(SortableStructure, Protocol):
    """Protocol describing the minimum interface required for clustering.

    Extends [SortableStructure][protein_quest.clustering.SortableStructure] with
    residue-range information needed for overlap-based clustering.

    Attributes:
        uniprot_start: Lowest UniProt residue position covered by the structure.
        uniprot_end: Highest UniProt residue position covered by the structure.
    """

    @property
    def uniprot_start(self) -> int: ...
    @property
    def uniprot_end(self) -> int: ...


def structure_overlap(a: ClusterableStructure, b: ClusterableStructure) -> int:
    """Number of overlapping UniProt residues between two structures.

    Both arguments must satisfy
    [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

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

    Both arguments must satisfy
    [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

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

    Both arguments must satisfy
    [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

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


def structure_distances[T: ClusterableStructure](structures: list[T]) -> list[float]:
    """Condensed (upper-triangle) pairwise distance matrix for a list of structures.

    Args:
        structures: Structures to compute distances for.
            Each structure must satisfy
            [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Returns:
        Condensed distance matrix as a flat list, suitable for input to
            [scipy.cluster.hierarchy.linkage][scipy.cluster.hierarchy.linkage].
    """
    condensed: list[float] = []
    for index, left in enumerate(structures[:-1]):
        condensed.extend(structure_distance(left, right) for right in structures[index + 1 :])
    return condensed


def _cluster_sort_key[T: ClusterableStructure](cluster: set[T]) -> tuple[int, int, int, str]:
    max_chain_length = max(member.chain_length for member in cluster)
    start, end, ident = min((member.uniprot_start, member.uniprot_end, member.id) for member in cluster)
    return (-max_chain_length, start, end, ident)


def structure_sort_key(member: SortableStructure) -> tuple[float, int, float, int, float, int, str]:
    """Deterministic quality sort key for a cluster member.

    1. Sequence identity descending (highest first)
    2. Resolution ascending (lowest first)
    3. Geometry quality descending (highest first; ``None`` sorts after valid values)
    4. Chain length descending (longest first)
    5. Identifier ascending (deterministic tie-break)

    A failing ``chain_length`` access (for example for PDB results with
    unparsable chain metadata) is treated as ``0`` so such entries can still
    be sorted alongside valid ones.

    Structures with lower resolution are preferred.
    If resolution is missing aka 0.0, those structures are undesirable.

    """
    try:
        chain_length = member.chain_length
    except Exception:  # noqa: BLE001 - sort-key fallback for adapters that derive chain_length lazily
        chain_length = 0
    resolution_kind: int = 0 if member.resolution_value != 0.0 else 1
    geometry_quality = member.geometry_quality
    geometry_quality_value = geometry_quality if geometry_quality is not None else 0.0
    geometry_quality_kind: int = 0 if geometry_quality is not None else 1
    return (
        -member.sequence_identity,
        resolution_kind,
        member.resolution_value,
        geometry_quality_kind,
        -geometry_quality_value,
        -chain_length,
        member.id,
    )


def sort_structures[T: SortableStructure](structures: Iterable[T]) -> list[T]:
    """Sort structures by quality criteria.

    See [structure_sort_key][protein_quest.clustering.structure_sort_key] for sort criteria.

    Args:
        structures: Structures to sort.
            Each structure must satisfy
            [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Returns:
        List of structures sorted by the criteria above.
    """
    return sorted(structures, key=structure_sort_key)


def hierarchical_clustering(condensed_distances: list[float]) -> np.ndarray:
    """Wrapper around [scipy.cluster.hierarchy.linkage][scipy.cluster.hierarchy.linkage] with complete method.

    Args:
        condensed_distances: Condensed distance matrix as a flat list

    Returns:
        Linkage matrix as returned by scipy's linkage function.
    """
    return linkage(condensed_distances, method="complete")


def flatten_hierarchical_clusters[T: ClusterableStructure](
    linkage_matrix: np.ndarray, structures: list[T]
) -> list[list[T]]:
    """Form flat clusters from a hierarchical clustering linkage matrix.

    Wrapper around [scipy.cluster.hierarchy.fcluster][scipy.cluster.hierarchy.fcluster] with
    distance criterion and [CLUSTER_DISTANCE_THRESHOLD][protein_quest.clustering.CLUSTER_DISTANCE_THRESHOLD].
    followed by
    mapping back to the original structures and sorting of clusters and their members.

    Clusters themselves are ordered by chain length descending, then by start
    and end residue, then identifier.

    Args:
        linkage_matrix: Linkage matrix as returned by scipy's linkage function.
        structures: Original list of structures corresponding to the distance matrix used to compute the linkage matrix.

    Returns:
        Sorted list of clusters with members also sorted.
    """
    cluster_ids = fcluster(linkage_matrix, t=CLUSTER_DISTANCE_THRESHOLD, criterion="distance")

    logger.info("Formed %d clusters from %d structures", len(set(cluster_ids)), len(structures))

    clusters_by_id: dict[int, set[T]] = {}
    for structure, cluster_id in zip(structures, cluster_ids, strict=True):
        clusters_by_id.setdefault(int(cluster_id), set()).add(structure)

    return [
        sort_structures(clusters_by_id[cluster_id])
        for cluster_id in sorted(clusters_by_id, key=lambda cid: _cluster_sort_key(clusters_by_id[cid]))
    ]


def cluster_structures_with_intermediates[T: ClusterableStructure](
    structures: list[T],
) -> tuple[list[list[T]], list[float], np.ndarray | None]:
    """Cluster structures and return reusable intermediate artifacts.

    Args:
        structures: Structures to cluster. All structures must have valid residue ranges;
            callers responsible for filtering out structures with missing
            range information beforehand.
            Each structure must satisfy
            [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Returns:
        Tuple of:
        1. Sorted clusters with sorted members.
        2. Condensed pairwise distances used for hierarchical clustering.
        3. Linkage matrix or ``None`` when fewer than two structures are provided.
    """
    if not structures:
        return [], [], None
    if len(structures) == 1:
        return [[structures[0]]], [], None

    condensed_distances = structure_distances(structures)
    linkage_matrix = hierarchical_clustering(condensed_distances)
    return flatten_hierarchical_clusters(linkage_matrix, structures), condensed_distances, linkage_matrix


def cluster_structures[T: ClusterableStructure](structures: list[T]) -> list[list[T]]:
    """Cluster structures by overlapping UniProt residue coverage.

    Args:
        structures: Structures to cluster. All structures must have valid residue ranges;
            callers responsible for filtering out structures with missing
            range information beforehand.
            Each structure must satisfy
            [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Returns:
        Sorted list of clusters with members also sorted.
    """
    clusters, _, _ = cluster_structures_with_intermediates(structures)
    return clusters


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


class ClusterCoverageError(ValueError):
    """Raised when not all clusters are represented in the top results."""

    def __init__(self, nr_clusters: int, top: int) -> None:
        msg = f"Not all {nr_clusters} clusters are represented in the top {top} results."
        super().__init__(msg)


def _is_each_cluster_represented_in_top[T](clusters: list[list[T]], top: int) -> None:
    nr_clusters = len(clusters)
    if nr_clusters > top:
        raise ClusterCoverageError(nr_clusters, top)


def top_members_of_clusters[T](clusters: list[list[T]], top: int) -> list[T]:
    """Return up to ``top`` members by interleaving cluster members round-robin.

    Args:
        clusters: Ordered clusters whose members are also sorted.
            First cluster and its first member is considered best.
        top: Maximum number of members to return.

    Returns:
        Interleaved members truncated to ``top``.

    Raises:
        ClusterCoverageError: If not all clusters are represented in the top results.
    """
    if top <= 0:
        msg = "Top must be a positive integer."
        raise ValueError(msg)
    _is_each_cluster_represented_in_top(clusters, top)
    return list(islice(interleave_longest(*clusters), top))


def filter_structures_on_clustered_resolution[T: ClusterableStructure](structures: list[T], top: int) -> list[T]:
    """Filter structures by resolution within residue-range clusters.

    Looks at how structures uniprot ranges overlap and clusters them by similarity of covered residue ranges.
    Then returns up to ``top`` structures by interleaving cluster members round-robin,

    Args:
        structures: Structures to filter. Must all have valid residue ranges.
            Each structure must satisfy
            [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.
        top: Number of top results to retain.

    Returns:
        Filtered list of up to ``top`` structures.

    Raises:
        ClusterCoverageError: If not all clusters are represented in the top results.
    """
    clusters = cluster_structures(structures)
    return top_members_of_clusters(clusters, top)
