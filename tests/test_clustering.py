from dataclasses import dataclass

import pytest

from protein_quest.clustering import (
    NO_OVERLAP_DISTANCE,
    cluster_structures,
    filter_structures_on_clustered_resolution,
    sort_structures,
    structure_distance,
    structure_union,
    top_members_of_clusters,
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
def test_structure_distance(a, b, expected):
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
def test_structure_union(a, b, expected):
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
        ],
    )
    def test_sort_order(self, structures, expected_order):
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


class TestClusterStructures:
    def test_empty_input_returns_empty_list(self):
        assert cluster_structures([]) == []

    def test_single_item_returns_singleton_cluster(self):
        structure = make_structure("1AAA", 1, 100)
        assert cluster_structures([structure]) == [[structure]]

    @pytest.mark.parametrize(
        "items,expected",
        [
            pytest.param(
                [
                    make_structure("2BBB", 50, 150),
                    make_structure("1AAA", 1, 100),
                ],
                [["1AAA", "2BBB"]],
                id="two_items_overlap_merge",
            ),
            pytest.param(
                [
                    make_structure("2BBB", 300, 400),
                    make_structure("1AAA", 1, 100),
                ],
                [["1AAA"], ["2BBB"]],
                id="two_items_no_overlap_split",
            ),
        ],
    )
    def test_two_item_shortcut_paths(self, items, expected):
        clusters = cluster_structures(items)
        assert [[member.id for member in cluster] for cluster in clusters] == expected

    def test_many_items_orders_clusters_and_members(self):
        # Cluster 1: longer chains (1-200), sorted by sequence identity then resolution.
        # Cluster 2: shorter chains (300-350), sorted by quality.
        items = [
            SimpleStructure("2BBB", 1, 200, resolution_value=3.0, sequence_identity=1.0, chain_length=200),
            SimpleStructure("1AAA", 1, 200, resolution_value=2.0, sequence_identity=1.0, chain_length=200),
            SimpleStructure("4DDD", 300, 350, resolution_value=2.5, sequence_identity=0.8, chain_length=51),
            SimpleStructure("3CCC", 300, 350, resolution_value=1.5, sequence_identity=0.9, chain_length=51),
        ]

        clusters = cluster_structures(items)

        assert [[member.id for member in cluster] for cluster in clusters] == [
            ["1AAA", "2BBB"],
            ["3CCC", "4DDD"],
        ]


class TestTopMembersOfClusters:
    def test_balanced_clusters_round_robin(self):
        result = top_members_of_clusters([["A1", "A2"], ["B1", "B2"], ["C1", "C2"]], top=6)
        assert result == ["A1", "B1", "C1", "A2", "B2", "C2"]

    def test_uneven_clusters_skip_exhausted_inputs(self):
        result = top_members_of_clusters([["A1"], ["B1", "B2", "B3"]], top=4)
        assert result == ["A1", "B1", "B2", "B3"]

    def test_matches_documented_example(self):
        result = top_members_of_clusters([[1, 2, 3], [4, 5], [6, 7, 8]], top=8)
        assert result == [1, 4, 6, 2, 5, 7, 3, 8]

    def test_empty_clusters_return_empty(self):
        assert top_members_of_clusters([[], [], []], top=3) == []

    def test_top_truncates_result(self):
        result = top_members_of_clusters([["A1", "A2"], ["B1", "B2"]], top=3)
        assert result == ["A1", "B1", "A2"]


class TestFilterStructuresOnClusteredResolution:
    @pytest.mark.parametrize("top", [0, -1])
    def test_non_positive_top_raises(self, top):
        with pytest.raises(ValueError, match="Top must be a positive integer"):
            filter_structures_on_clustered_resolution([make_structure("1AAA", 1, 100)], top=top)

    def test_interleaved_cluster_order_smoke(self):
        items = [
            SimpleStructure("A1", 1, 200, resolution_value=1.0, sequence_identity=1.0, chain_length=200),
            SimpleStructure("A2", 1, 200, resolution_value=2.0, sequence_identity=0.9, chain_length=200),
            SimpleStructure("B1", 300, 360, resolution_value=1.5, sequence_identity=1.0, chain_length=61),
            SimpleStructure("B2", 300, 360, resolution_value=2.5, sequence_identity=0.9, chain_length=61),
        ]

        filtered = filter_structures_on_clustered_resolution(items, top=3)
        assert [member.id for member in filtered] == ["A1", "B1", "A2"]

    def test_top_above_available_returns_all_without_duplicates(self):
        items = [
            make_structure("A1", 1, 100),
            make_structure("A2", 1, 100),
            make_structure("B1", 300, 350),
        ]

        filtered = filter_structures_on_clustered_resolution(items, top=10)
        filtered_ids = [member.id for member in filtered]

        assert filtered_ids == ["A1", "B1", "A2"]
        assert len(filtered_ids) == len(set(filtered_ids))
