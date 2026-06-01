from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from protein_quest.clustering import hierarchical_clustering, structure_distances
from protein_quest.clustering_io import cluster_results_by_accession, linkage_to_newick, write_condensed_distances_csv
from protein_quest.filters.resolution import ResolutionFilterStatistics


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


def make_stat(ident: str, accession: str | None, start: int, end: int) -> ResolutionFilterStatistics:
    return ResolutionFilterStatistics(
        input_file=Path(f"{ident}.cif"),
        id=ident,
        uniprot_accession=accession,
        resolution=1.5,
        total_residue_count=end - start + 1,
        is_alphafold=False,
        uniprot_start=start,
        uniprot_end=end,
        sequence_identity=1.0,
        chain_length=end - start + 1,
        passed=False,
        output_file=None,
    )


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


def test_cluster_results_by_accession_populates_artifacts_once():
    stats = [
        make_stat("A", "P11111", 1, 50),
        make_stat("B", "P11111", 25, 75),
        make_stat("C", "P22222", 100, 120),
    ]

    results = cluster_results_by_accession(stats)

    assert [result.uniprot_accession for result in results] == ["P11111", "P22222"]

    p11111 = results[0]
    assert len(p11111.clusters) == 1
    assert len(p11111.structures) == 2
    assert p11111.linkage_matrix is not None
    assert p11111.condensed_distances == structure_distances(p11111.structures)

    p22222 = results[1]
    assert p22222.clusters == [[p22222.structures[0]]]
    assert p22222.condensed_distances == []
    assert p22222.linkage_matrix is None


def test_cluster_results_by_accession_skips_missing_accession(caplog: pytest.LogCaptureFixture):
    stats = [
        make_stat("A", None, 1, 10),
        make_stat("B", "P33333", 1, 10),
    ]

    results = cluster_results_by_accession(stats)

    assert len(results) == 1
    assert results[0].uniprot_accession == "P33333"
    assert "has no UniProt accession" in caplog.text


def test_write_condensed_distances_csv_validates_precomputed_length(tmp_path: Path):
    output = tmp_path / "distances.csv"
    structures = [
        make_structure("A", 1, 10),
        make_structure("B", 1, 10),
        make_structure("C", 20, 30),
    ]

    with pytest.raises(ValueError, match="Condensed distance length"):
        write_condensed_distances_csv([0.1], structures, output)
