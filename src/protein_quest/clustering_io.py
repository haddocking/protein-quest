"""Helpers for writing clustering outputs for CLI commands."""

import csv
import logging
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import ClusterNode, to_tree
from tqdm.auto import tqdm

from protein_quest.clustering import ClusterableStructure, cluster_structures_with_intermediates
from protein_quest.filters.resolution import ResolutionFilterStatistics

logger = logging.getLogger(__name__)


@dataclass
class AccessionClusters:
    """Clusters for a single UniProt accession."""

    uniprot_accession: str
    structures: list[ResolutionFilterStatistics]
    clusters: list[list[ResolutionFilterStatistics]]
    condensed_distances: list[float]
    linkage_matrix: np.ndarray | None


def _flatten_accession_clusters(result: AccessionClusters) -> list[ResolutionFilterStatistics]:
    """Flatten clustered structures in cluster/rank order.

    Args:
        result: Clustered structures for one accession.

    Returns:
        Structures in cluster order, preserving rank order within each cluster.
    """
    return [member for cluster in result.clusters for member in cluster]


def cluster_results_by_accession(stats: Iterable[ResolutionFilterStatistics]) -> list[AccessionClusters]:
    """Group and cluster structures by UniProt accession.

    Args:
        stats: Structure statistics to group and cluster.

    Returns:
        Cluster results grouped by accession, sorted by accession.
    """
    grouped: dict[str, list[ResolutionFilterStatistics]] = {}
    for stat in stats:
        if not stat.uniprot_accession:
            logger.warning("Skipping %s because it has no UniProt accession.", stat.input_file)
            continue
        grouped.setdefault(stat.uniprot_accession, []).append(stat)

    results: list[AccessionClusters] = []
    for accession, group_results in tqdm(sorted(grouped.items()), unit="acc"):
        clusters, condensed_distances, linkage_matrix = cluster_structures_with_intermediates(group_results)
        results.append(
            AccessionClusters(
                uniprot_accession=accession,
                structures=group_results,
                clusters=clusters,
                condensed_distances=condensed_distances,
                linkage_matrix=linkage_matrix,
            )
        )
    return results


def write_clusters_csv(results: list[AccessionClusters], output: Path) -> None:
    """Write structure-level cluster assignments.

    Args:
        results: Cluster results to serialize.
        output: Destination path for the CSV file.

    Returns:
        None.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "uniprot_accession",
                "cluster_id",
                "rank_in_cluster",
                "structure_id",
                "input_file",
                "resolution",
                "sequence_identity",
                "chain_length",
                "uniprot_start",
                "uniprot_end",
            ],
        )
        writer.writeheader()

        for result in results:
            for cluster_id, cluster in enumerate(result.clusters, start=1):
                for rank, structure in enumerate(cluster, start=1):
                    writer.writerow(
                        {
                            "uniprot_accession": result.uniprot_accession,
                            "cluster_id": cluster_id,
                            "rank_in_cluster": rank,
                            "structure_id": structure.id,
                            "input_file": structure.input_file,
                            "resolution": structure.resolution,
                            "sequence_identity": f"{structure.sequence_identity:.3f}",
                            "chain_length": structure.chain_length,
                            "uniprot_start": structure.uniprot_start,
                            "uniprot_end": structure.uniprot_end,
                        }
                    )


def write_stats_csv(results: list[AccessionClusters], output: Path) -> None:
    """Write accession-level clustering statistics.

    Args:
        results: Cluster results to summarize.
        output: Destination path for the CSV file.

    Returns:
        None.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "uniprot_accession",
                "num_structures",
                "num_clusters",
                "min_uniprot_start",
                "max_uniprot_end",
                "max_cluster_size",
                "min_resolution",
                "max_sequence_identity",
                "cluster_sizes",
            ],
        )
        writer.writeheader()

        for result in results:
            members = _flatten_accession_clusters(result)
            non_zero_resolutions = [member.resolution for member in members if member.resolution != 0.0]
            writer.writerow(
                {
                    "uniprot_accession": result.uniprot_accession,
                    "num_structures": len(members),
                    "num_clusters": len(result.clusters),
                    "min_uniprot_start": min(member.uniprot_start for member in members),
                    "max_uniprot_end": max(member.uniprot_end for member in members),
                    "max_cluster_size": max(len(cluster) for cluster in result.clusters),
                    "min_resolution": min(non_zero_resolutions) if non_zero_resolutions else 0.0,
                    "max_sequence_identity": max(member.sequence_identity for member in members),
                    "cluster_sizes": ";".join(str(len(cluster)) for cluster in result.clusters),
                }
            )


def write_condensed_distances_csv[T: ClusterableStructure](
    condensed_distances: list[float], structures: list[T], output: Path
) -> None:
    """Write pairwise condensed distances in long CSV form.

    Args:
        condensed_distances: Precomputed condensed distances matching ``structures`` order.
            See [protein_quest.clustering.structure_distances][] how to compute these distances.
        structures: Structures corresponding to ``condensed_distances`` order.
        output: Destination path for the CSV file.

    Returns:
        None.
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    # Sanity check
    expected_len_distances = (len(structures) * (len(structures) - 1)) // 2
    if len(condensed_distances) != expected_len_distances:
        msg = (
            "Condensed distance length does not match structure count: "
            f"expected {expected_len_distances}, got {len(condensed_distances)}."
        )
        raise ValueError(msg)

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["structure_i", "structure_j", "distance"])
        writer.writeheader()

        index = 0
        for i, left in enumerate(structures[:-1]):
            for right in structures[i + 1 :]:
                writer.writerow(
                    {
                        "structure_i": left.id,
                        "structure_j": right.id,
                        "distance": condensed_distances[index],
                    }
                )
                index += 1


def _node_label[T: ClusterableStructure](node_id: int, structures: list[T]) -> str:
    """Return a human-readable label for a linkage node id.

    Args:
        node_id: Node index from a scipy linkage matrix.
        structures: Input structures used to build the linkage matrix.

    Returns:
        Structure id for leaf nodes, otherwise an internal node label.
    """
    if node_id < len(structures):
        return structures[node_id].id
    return f"internal_{node_id}"


def write_linkage_matrix_csv[T: ClusterableStructure](
    linkage_matrix: np.ndarray, structures: list[T], output: Path
) -> None:
    """Write a linkage matrix with human-readable node labels.

    Args:
        linkage_matrix: Linkage matrix produced by hierarchical clustering.
        structures: Structures corresponding to leaf rows in the linkage matrix.
        output: Destination path for the CSV file.

    Returns:
        None.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["merge_id", "left", "right", "distance", "count", "left_label", "right_label"],
        )
        writer.writeheader()

        nr_leaves = len(structures)
        for row_index, row in enumerate(linkage_matrix):
            left = int(row[0])
            right = int(row[1])
            writer.writerow(
                {
                    "merge_id": nr_leaves + row_index,
                    "left": left,
                    "right": right,
                    "distance": float(row[2]),
                    "count": int(row[3]),
                    "left_label": _node_label(left, structures),
                    "right_label": _node_label(right, structures),
                }
            )


def _to_newick[T: ClusterableStructure](node: ClusterNode, parent_distance: float, structures: list[T]) -> str:
    """Convert a clustering node and descendants to a Newick subtree.

    Args:
        node: Current clustering node.
        parent_distance: Distance of the parent node.
        structures: Input structures used to build the clustering tree.

    Returns:
        Newick-formatted subtree rooted at ``node``.

    Raises:
        ValueError: If an internal node is missing either child.
    """
    branch_length = max(parent_distance - float(node.dist), 0.0)

    if node.is_leaf():
        return f"{structures[int(node.id)].id}:{branch_length:.6f}"

    left_node = node.get_left()
    right_node = node.get_right()
    if left_node is None or right_node is None:
        msg = f"Internal clustering node {int(node.id)} is missing a child."
        raise ValueError(msg)

    left = _to_newick(left_node, float(node.dist), structures)
    right = _to_newick(right_node, float(node.dist), structures)
    return f"({left},{right}):{branch_length:.6f}"


def linkage_to_newick[T: ClusterableStructure](linkage_matrix: np.ndarray, structures: list[T]) -> str:
    """Convert a scipy linkage matrix into Newick text.

    Args:
        linkage_matrix: Linkage matrix produced by hierarchical clustering.
        structures: Structures corresponding to leaf rows in the linkage matrix.

    Returns:
        Newick-formatted tree string.

    Raises:
        ValueError: If the root node is missing either child.
    """
    if not structures:
        return ";"
    if len(structures) == 1:
        return f"{structures[0].id};"

    root = to_tree(linkage_matrix, rd=False)
    left_node = root.get_left()
    right_node = root.get_right()
    if left_node is None or right_node is None:
        msg = "Root clustering node is missing a child."
        raise ValueError(msg)

    left = _to_newick(left_node, float(root.dist), structures)
    right = _to_newick(right_node, float(root.dist), structures)
    return f"({left},{right});"


def write_dendrogram_nwk[T: ClusterableStructure](
    linkage_matrix: np.ndarray, structures: list[T], output: Path
) -> None:
    """Write a dendrogram in Newick format.

    Args:
        linkage_matrix: Linkage matrix produced by hierarchical clustering.
        structures: Structures corresponding to leaf rows in the linkage matrix.
        output: Destination path for the Newick file.

    Returns:
        None.

    Raises:
        ValueError: If the linkage tree contains an internal node missing a child.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(linkage_to_newick(linkage_matrix, structures), encoding="utf-8")
