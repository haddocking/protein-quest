from dataclasses import dataclass

import pytest

from protein_quest.clustering import (
    NO_OVERLAP_DISTANCE,
    ClusterableStructure,
    ClusterCoverageError,
    cluster_structures,
    filter_structures_on_clustered_resolution,
    sort_structures,
    structure_distance,
    structure_union,
    top_members_across_clusters,
    top_members_per_cluster,
)


@dataclass(frozen=True)
class SimpleStructure:
    """Minimal `ClusterableStructure` used by the generic clustering tests."""

    id: str
    uniprot_start: int
    uniprot_end: int
    resolution_value: float = 0.0
    sequence_identity: float = 1.0
    chain_length: int = 0


def make_structure(ident: str, start: int, end: int) -> SimpleStructure:
    return SimpleStructure(id=ident, uniprot_start=start, uniprot_end=end, chain_length=end - start + 1)


@dataclass(frozen=True)
class BrokenChainLengthStructure:
    id: str
    uniprot_start: int
    uniprot_end: int
    resolution_value: float = 0.0
    sequence_identity: float = 1.0

    @property
    def chain_length(self) -> int:
        msg = "chain length unavailable"
        raise RuntimeError(msg)


@pytest.mark.parametrize(
    "a,b,expected",
    [
        pytest.param(make_structure("1AAA", 1, 1), make_structure("2BBB", 1, 1), 0.0, id="one_residue_identical"),
        pytest.param(make_structure("1AAA", 1, 250), make_structure("2BBB", 1, 250), 0.0, id="identical_ranges"),
        pytest.param(
            make_structure("1AAA", 1, 250), make_structure("3CCC", 1, 54), 1 - 54 / 250, id="partial_overlap_subset"
        ),
        pytest.param(
            make_structure("1AAA", 1, 250), make_structure("4DDD", 300, 400), NO_OVERLAP_DISTANCE, id="no_overlap"
        ),
        pytest.param(
            make_structure("1AAA", 50, 100), make_structure("3CCC", 1, 100), 1 - 51 / 100, id="subset_overlap"
        ),
        pytest.param(
            make_structure("1AAA", 50, 100), make_structure("4DDD", 50, 120), 1 - 51 / 71, id="partial_overlap"
        ),
    ],
)
def test_structure_distance(a: ClusterableStructure, b: ClusterableStructure, expected: float):
    assert structure_distance(a, b) == expected


@pytest.mark.parametrize(
    "a,b,expected",
    [
        pytest.param(make_structure("1AAA", 1, 1), make_structure("2BBB", 1, 1), 1, id="union_one_residue"),
        pytest.param(make_structure("1AAA", 1, 250), make_structure("2BBB", 1, 250), 250, id="union_identical_ranges"),
        pytest.param(make_structure("1AAA", 1, 250), make_structure("3CCC", 1, 54), 250, id="union_subset"),
        pytest.param(make_structure("1AAA", 1, 100), make_structure("2BBB", 50, 150), 150, id="union_partial_overlap"),
        pytest.param(make_structure("1AAA", 1, 100), make_structure("2BBB", 300, 400), 201, id="union_no_overlap"),
    ],
)
def test_structure_union(a: ClusterableStructure, b: ClusterableStructure, expected: int):
    assert structure_union(a, b) == expected


class TestSortStructures:
    @pytest.mark.parametrize(
        "structures, expected_order",
        [
            pytest.param(
                {
                    SimpleStructure("3CCC", 1, 250, sequence_identity=0.8, chain_length=250),
                    SimpleStructure("1AAA", 1, 250, sequence_identity=1.0, chain_length=250),
                },
                ["1AAA", "3CCC"],
                id="sort_by_sequence_identity",
            ),
            pytest.param(
                {
                    SimpleStructure("2BBB", 1, 250, resolution_value=5.4, chain_length=250),
                    SimpleStructure("3CCC", 1, 250, resolution_value=2.1, chain_length=250),
                },
                ["3CCC", "2BBB"],
                id="sort_by_resolution",
            ),
            pytest.param(
                {
                    SimpleStructure("2BBB", 1, 100, chain_length=100),
                    SimpleStructure("3CCC", 1, 200, chain_length=200),
                },
                ["3CCC", "2BBB"],
                id="sort_by_chain_length",
            ),
            pytest.param(
                {
                    SimpleStructure("2BBB", 1, 250, chain_length=250),
                    SimpleStructure("1AAA", 1, 250, chain_length=250),
                },
                ["1AAA", "2BBB"],
                id="sort_by_id_tiebreak",
            ),
            pytest.param(
                {
                    SimpleStructure("2BBB", 1, 250, resolution_value=0.0, chain_length=250),
                    SimpleStructure("3CCC", 1, 250, resolution_value=2.1, chain_length=250),
                },
                ["3CCC", "2BBB"],
                id="sort_by_resolution_unset_last",
            ),
        ],
    )
    def test_sort_order(self, structures: list[ClusterableStructure], expected_order: list[str]):
        sorted_members = sort_structures(structures)
        sorted_ids = [member.id for member in sorted_members]
        assert sorted_ids == expected_order

    def test_chain_length_failure_falls_back_to_zero(self):
        structures = [
            BrokenChainLengthStructure("2BBB", 1, 250, resolution_value=1.0, sequence_identity=1.0),
            SimpleStructure("1AAA", 1, 250, resolution_value=1.0, sequence_identity=1.0, chain_length=250),
        ]

        sorted_ids = [member.id for member in sort_structures(structures)]
        assert sorted_ids == ["1AAA", "2BBB"]


@pytest.mark.parametrize(
    "structures,expected",
    [
        pytest.param([], [], id="empty_input_returns_empty_list"),
        pytest.param(
            [make_structure("1AAA", 1, 100)],
            [["1AAA"]],
            id="single_structure_returns_singleton_cluster",
        ),
        pytest.param(
            [
                make_structure("2BBB", 51, 150),
                make_structure("1AAA", 1, 100),
            ],
            [["1AAA", "2BBB"]],
            id="overlap_same_length_1cluster",
        ),
        pytest.param(
            [
                make_structure("2BBB", 75, 150),
                make_structure("1AAA", 1, 100),
            ],
            [["1AAA", "2BBB"]],
            id="overlap_diff_length_1cluster",
        ),
        pytest.param(
            [
                make_structure("2BBB", 301, 400),
                make_structure("1AAA", 1, 100),
            ],
            [["1AAA"], ["2BBB"]],
            id="no_overlap_same_length_2clusters",
        ),
        pytest.param(
            [
                make_structure("2BBB", 321, 400),
                make_structure("1AAA", 1, 100),
            ],
            [["1AAA"], ["2BBB"]],
            id="no_overlap_diff_length_2clusters",
        ),
        pytest.param(
            [
                SimpleStructure("2BBB", 1, 200, resolution_value=3.0, sequence_identity=1.0, chain_length=200),
                SimpleStructure("1AAA", 1, 200, resolution_value=2.0, sequence_identity=1.0, chain_length=200),
                SimpleStructure("4DDD", 300, 350, resolution_value=2.5, sequence_identity=0.8, chain_length=51),
                SimpleStructure("3CCC", 300, 350, resolution_value=1.5, sequence_identity=0.9, chain_length=51),
            ],
            [["1AAA", "2BBB"], ["3CCC", "4DDD"]],
            id="many_structures_orders_clusters_and_members",
        ),
        # First example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                SimpleStructure("2BBB", 1, 250, resolution_value=3.7, sequence_identity=1.0, chain_length=250),
                SimpleStructure("1AAA", 1, 250, resolution_value=3.6, sequence_identity=1.0, chain_length=250),
                SimpleStructure("3CCC", 1, 250, resolution_value=2.1, sequence_identity=1.0, chain_length=250),
                SimpleStructure("4DDD", 300, 400, resolution_value=8.1, sequence_identity=1.0, chain_length=101),
                SimpleStructure("5EEE", 300, 400, resolution_value=4.6, sequence_identity=1.0, chain_length=101),
                SimpleStructure("6FFF", 500, 1000, resolution_value=1.3, sequence_identity=1.0, chain_length=501),
                SimpleStructure("7GGG", 500, 1000, resolution_value=1.4, sequence_identity=1.0, chain_length=501),
                SimpleStructure("8HHH", 500, 1000, resolution_value=1.6, sequence_identity=1.0, chain_length=501),
            ],
            [
                ["6FFF", "7GGG", "8HHH"],
                [
                    "3CCC",
                    "1AAA",
                    "2BBB",
                ],
                [
                    "5EEE",
                    "4DDD",
                ],
            ],
            id="three_domains_split",
        ),
        # Second example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                SimpleStructure("1AAA", 1, 250, resolution_value=3.6, sequence_identity=1.0, chain_length=250),
                SimpleStructure("4DDD", 200, 400, resolution_value=8.1, sequence_identity=1.0, chain_length=201),
                SimpleStructure("6FFF", 500, 1000, resolution_value=1.3, sequence_identity=1.0, chain_length=501),
                SimpleStructure("9III", 1, 600, resolution_value=4.2, sequence_identity=1.0, chain_length=600),
                SimpleStructure("10JJJ", 1, 1000, resolution_value=1.4, sequence_identity=1.0, chain_length=1000),
            ],
            [["10JJJ", "1AAA", "9III", "4DDD"], ["6FFF"]],
            id="overlap_merges",
        ),
        pytest.param(
            [
                SimpleStructure("1AAA", 1, 100, resolution_value=1.0, sequence_identity=1.0, chain_length=100),
                SimpleStructure("2BBB", 100, 200, resolution_value=2.0, sequence_identity=1.0, chain_length=100),
            ],
            [["1AAA", "2BBB"]],
            id="boundary_overlap_merges",
        ),
        pytest.param(
            [
                SimpleStructure("1AAA", 1, 100, resolution_value=1.0, sequence_identity=1.0, chain_length=100),
                SimpleStructure("2BBB", 101, 200, resolution_value=2.0, sequence_identity=1.0, chain_length=100),
            ],
            [["1AAA"], ["2BBB"]],
            id="adjacent_no_overlap_separates",
        ),
        pytest.param(
            [
                SimpleStructure("1AAA", 1, 100, resolution_value=1.0, sequence_identity=1.0, chain_length=100),
                SimpleStructure("2BBB", 50, 150, resolution_value=2.0, sequence_identity=1.0, chain_length=101),
                SimpleStructure("3CCC", 120, 200, resolution_value=3.0, sequence_identity=1.0, chain_length=81),
            ],
            [["1AAA", "2BBB"], ["3CCC"]],
            id="complete_linkage_bridge_splits",
        ),
        pytest.param(
            [
                SimpleStructure("1AAA", 1, 250, resolution_value=3.6, sequence_identity=1.0, chain_length=250),
                SimpleStructure("2BBB", 2, 249, resolution_value=5.4, sequence_identity=1.0, chain_length=248),
                SimpleStructure("3CCC", 3, 248, resolution_value=2.1, sequence_identity=1.0, chain_length=246),
            ],
            [["3CCC", "1AAA", "2BBB"]],
            id="complete_linkage_bridge_splits",
        ),
    ],
)
def test_cluster_structures(structures: list[ClusterableStructure], expected: list[list[str]]):
    clusters = cluster_structures(structures)
    cluster_ids = [[member.id for member in cluster] for cluster in clusters]
    assert cluster_ids == expected


class TestTopMembersAcrossClusters:
    def test_balanced_clusters_round_robin(self):
        result = top_members_across_clusters([["A1", "A2"], ["B1", "B2"], ["C1", "C2"]], top=6)
        assert result == ["A1", "B1", "C1", "A2", "B2", "C2"]

    def test_uneven_clusters_skip_exhausted_inputs(self):
        result = top_members_across_clusters([["A1"], ["B1", "B2", "B3"]], top=4)
        assert result == ["A1", "B1", "B2", "B3"]

    def test_matches_documented_example(self):
        result = top_members_across_clusters([[1, 2, 3], [4, 5], [6, 7, 8]], top=8)
        assert result == [1, 4, 6, 2, 5, 7, 3, 8]

    def test_empty_clusters_return_empty(self):
        assert top_members_across_clusters([[], [], []], top=3) == []

    def test_top_truncates_result(self):
        result = top_members_across_clusters([["A1", "A2"], ["B1", "B2"]], top=3)
        assert result == ["A1", "B1", "A2"]

    def test_clusters_under_represented_in_top_raises(self):
        with pytest.raises(ClusterCoverageError, match="Not all 2 clusters are represented in the top 1 results"):
            top_members_across_clusters([["A2"], ["B1"]], top=1)


class TestTopMembersPerCluster:
    def test_balanced_clusters_round_robin(self):
        result = top_members_per_cluster([["A1", "A2"], ["B1", "B2"], ["C1", "C2"]], top=2)
        assert result == ["A1", "B1", "C1", "A2", "B2", "C2"]

    def test_top_limits_each_cluster_before_interleaving(self):
        result = top_members_per_cluster([["A1", "A2", "A3"], ["B1", "B2"]], top=2)
        assert result == ["A1", "B1", "A2", "B2"]

    def test_uneven_clusters_skip_exhausted_inputs(self):
        result = top_members_per_cluster([["A1"], ["B1", "B2", "B3"]], top=2)
        assert result == ["A1", "B1", "B2"]

    def test_empty_clusters_return_empty(self):
        assert top_members_per_cluster([[], [], []], top=3) == []

    @pytest.mark.parametrize("top", [0, -1])
    def test_non_positive_top_raises(self, top: int):
        with pytest.raises(ValueError, match="Top must be a positive integer"):
            top_members_per_cluster([["A1"]], top=top)


@pytest.fixture
def sample_structures() -> list[SimpleStructure]:
    return [
        SimpleStructure("A1", 1, 200, resolution_value=1.0, sequence_identity=1.0, chain_length=200),
        SimpleStructure("A2", 1, 200, resolution_value=2.0, sequence_identity=0.9, chain_length=200),
        SimpleStructure("B1", 300, 360, resolution_value=1.5, sequence_identity=1.0, chain_length=61),
        SimpleStructure("B2", 300, 360, resolution_value=2.5, sequence_identity=0.9, chain_length=61),
    ]


class TestFilterStructuresOnClusteredResolution:
    @pytest.mark.parametrize("top", [0, -1])
    def test_non_positive_top_raises(self, top: int):
        with pytest.raises(ValueError, match="Top must be a positive integer"):
            filter_structures_on_clustered_resolution([make_structure("1AAA", 1, 100)], top=top)

    def test_interleaved_cluster_order_smoke(self, sample_structures: list[SimpleStructure]):
        filtered = filter_structures_on_clustered_resolution(sample_structures, top=1)
        assert [member.id for member in filtered] == ["A1", "B1"]

    def test_across_clusters_selection_strategy(self, sample_structures: list[SimpleStructure]):
        filtered = filter_structures_on_clustered_resolution(
            sample_structures, top=3, selection_strategy="across_clusters_top"
        )
        assert [member.id for member in filtered] == ["A1", "B1", "A2"]

    def test_per_cluster_selection_strategy(self, sample_structures: list[SimpleStructure]):
        filtered = filter_structures_on_clustered_resolution(
            sample_structures, top=1, selection_strategy="per_cluster_top"
        )
        assert [member.id for member in filtered] == ["A1", "B1"]

    def test_top_above_available_returns_all_without_duplicates(self):
        structures = [
            make_structure("A1", 1, 100),
            make_structure("A2", 1, 100),
            make_structure("B1", 300, 350),
        ]

        filtered = filter_structures_on_clustered_resolution(structures, top=10)
        filtered_ids = [member.id for member in filtered]

        assert filtered_ids == ["A1", "B1", "A2"]
        assert len(filtered_ids) == len(set(filtered_ids))
