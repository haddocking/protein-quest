from dataclasses import dataclass

import numpy as np

from protein_quest.clustering import hierarchical_clustering, structure_distances
from protein_quest.clustering_io import linkage_to_newick


@dataclass(frozen=True)
class SimpleStructure:
    id: str
    uniprot_start: int
    uniprot_end: int
    resolution_value: float = 0.0
    sequence_identity: float = 1.0
    chain_length: int = 0


def make_structure(ident: str, start: int, end: int) -> SimpleStructure:
    return SimpleStructure(id=ident, uniprot_start=start, uniprot_end=end, chain_length=end - start + 1)


def test_linkage_to_newick_with_intermediate_nodes_manual_linkage():
    structures = [
        make_structure("A", 1, 10),
        make_structure("B", 1, 10),
        make_structure("C", 50, 60),
        make_structure("D", 50, 60),
    ]
    linkage_matrix = np.array(
        [
            [0, 1, 0.1, 2],
            [2, 3, 0.2, 2],
            [4, 5, 0.5, 4],
        ],
        dtype=float,
    )

    newick = linkage_to_newick(linkage_matrix, structures)

    assert newick == "((A:0.100000,B:0.100000):0.400000,(C:0.200000,D:0.200000):0.300000);"
    assert newick.count("(") == 3


def test_linkage_to_newick_from_distances_has_intermediate_nodes():
    structures = [
        make_structure("A", 1, 10),
        make_structure("B", 2, 11),
        make_structure("C", 100, 110),
        make_structure("D", 101, 111),
    ]

    condensed = structure_distances(structures)
    linkage_matrix = hierarchical_clustering(condensed)
    newick = linkage_to_newick(linkage_matrix, structures)

    assert newick.endswith(";")
    assert newick.count("(") == len(structures) - 1
    for structure in structures:
        assert f"{structure.id}:" in newick
