import pytest

from protein_quest.pdbe.clustering import (
    NO_OVERLAP_DISTANCE,
    cluster_pdbs,
    filter_pdbs_on_clustered_resolution,
    pdb_distance,
    pdb_union,
    sort_pdbs,
)
from protein_quest.pdbe.result import PdbResult


def make_pdb(pdb_id: str, uniprot_chains: str, resolution: float = 5.0) -> PdbResult:
    return PdbResult(
        id=pdb_id, method="X-Ray_Crystallography", resolution=str(resolution), uniprot_chains=uniprot_chains
    )


@pytest.mark.parametrize(
    "a,b,expected",
    [
        pytest.param(
            make_pdb("1AAA", "A=1-1", 3.6),
            make_pdb("2BBB", "A=1-1", 5.4),
            0.0,
            id="one_residue_identical",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("2BBB", "A=1-250", 5.4),
            0.0,
            id="identical_ranges",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("3CCC", "A=1-54", 2.1),
            1 - 54 / 250,
            id="partial_overlap_subset",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("4DDD", "A=300-400", 8.1),
            NO_OVERLAP_DISTANCE,
            id="no_overlap",
        ),
        pytest.param(
            make_pdb("1AAA", "A=50-100", 3.6),
            make_pdb("3CCC", "A=1-100", 2.1),
            1 - 51 / 100,
            id="subset_overlap",
        ),
        pytest.param(
            make_pdb("1AAA", "A=50-100", 3.6),
            make_pdb("4DDD", "A=50-120", 8.1),
            1 - 51 / 71,
            id="partial_overlap",
        ),
    ],
)
def test_pdb_distance(a, b, expected):
    assert pdb_distance(a, b) == expected


@pytest.mark.parametrize(
    "a,b,expected",
    [
        pytest.param(
            make_pdb("1AAA", "A=1-1", 3.6),
            make_pdb("2BBB", "A=1-1", 5.4),
            1,
            id="union_one_residue",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("2BBB", "A=1-250", 5.4),
            250,
            id="union_identical_ranges",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("3CCC", "A=1-54", 2.1),
            250,
            id="union_subset",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-100", 3.6),
            make_pdb("2BBB", "A=50-150", 5.4),
            150,
            id="union_partial_overlap",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-100", 3.6),
            make_pdb("2BBB", "A=300-400", 5.4),
            201,
            id="union_no_overlap",
        ),
    ],
)
def test_pdb_union(a, b, expected):
    assert pdb_union(a, b) == expected


@pytest.mark.parametrize(
    "pdbs, expected_clusters",
    [
        pytest.param([], ([], set()), id="empty_input"),
        pytest.param(
            [make_pdb("1AAA", "A=1-250", 3.6)],
            ([{"1AAA"}], set()),
            id="single_valid_shortcut",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-250", 2.1),
            ],
            ([{"1AAA", "2BBB", "3CCC"}], set()),
            id="same_range_merges",
        ),
        # First example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-250", 2.1),
                make_pdb("4DDD", "A=300-400", 8.1),
                make_pdb("5EEE", "A=300-400", 4.6),
                make_pdb("6FFF", "A=500-1000", 1.3),
                make_pdb("7GGG", "A=500-1000", 1.4),
                make_pdb("8HHH", "A=500-1000", 1.6),
            ],
            ([{"1AAA", "2BBB", "3CCC"}, {"4DDD", "5EEE"}, {"6FFF", "7GGG", "8HHH"}], set()),
            id="three_domains_split",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-54,A=90-250", 2.1),
            ],
            ([{"1AAA", "2BBB", "3CCC"}], set()),
            id="split_range_merges",
        ),
        pytest.param(
            [
                make_pdb("1BAD", "A=-", 1.0),
                make_pdb("2BAD", "A=-", 2.0),
            ],
            ([], {"1BAD", "2BAD"}),
            id="all_invalid_shortcut",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("1BAD", "A=-", 2.0),
            ],
            ([{"1AAA"}], {"1BAD"}),
            id="mixed_valid_invalid_separates",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("2BBB", "A=100-200", 2.0),
            ],
            ([{"1AAA", "2BBB"}], set()),
            id="boundary_overlap_merges",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("2BBB", "A=101-200", 2.0),
            ],
            ([{"1AAA"}, {"2BBB"}], set()),
            id="adjacent_no_overlap_separates",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("2BBB", "A=50-150", 2.0),
                make_pdb("3CCC", "A=120-200", 3.0),
            ],
            ([{"1AAA", "2BBB"}, {"3CCC"}], set()),
            id="complete_linkage_bridge_splits",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=2-249", 5.4),
                make_pdb("3CCC", "A=3-248", 2.1),
            ],
            ([{"1AAA", "2BBB", "3CCC"}], set()),
            id="nested_ranges_merge",
        ),
        # Second example from https://github.com/haddocking/protein-quest/issues/102
        # TODO are expected clusters what we want?
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("4DDD", "A=200-400", 8.1),
                make_pdb("6FFF", "A=500-1000", 1.3),
                make_pdb("9III", "A=1-600", 4.2),
                make_pdb("10JJJ", "A=1-1000", 1.4),
            ],
            ([{"1AAA", "4DDD", "9III", "10JJJ"}, {"6FFF"}], set()),
            id="overlap_merges",
        ),
    ],
)
def test_cluster_pdbs(pdbs, expected_clusters):
    expected_valid_member_ids, expected_invalid_member_ids = expected_clusters
    valid_clusters, invalid_cluster = cluster_pdbs(pdbs)

    cluster_member_ids = {frozenset(p.id for p in cluster) for cluster in valid_clusters}
    expected = {frozenset(ids) for ids in expected_valid_member_ids}
    assert cluster_member_ids == expected

    invalid_member_ids = {p.id for p in invalid_cluster}
    assert invalid_member_ids == expected_invalid_member_ids


@pytest.mark.parametrize(
    "pdbs, expected_order",
    [
        pytest.param(
            {
                make_pdb("3CCC", "A=1-50, 200-250", 2.1),
                make_pdb("1AAA", "A=1-250", 3.6),
            },
            ["1AAA", "3CCC"],
            id="sort_by_sequence_identity",
        ),
        pytest.param(
            {
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-250", 2.1),
            },
            ["3CCC", "2BBB"],
            id="sort_by_resolution",
        ),
        pytest.param(
            {
                make_pdb("2BBB", "A=1-100", 2.1),
                make_pdb("3CCC", "A=1-200", 2.1),
            },
            ["3CCC", "2BBB"],
            id="sort_by_chain_length",
        ),
        pytest.param(
            {
                make_pdb("2BBB", "A=1-250", 3.6),
                make_pdb("1AAA", "A=1-250", 3.6),
            },
            ["1AAA", "2BBB"],
            id="sort_by_id_tiebreak",
        ),
    ],
)
def test_sort_pdbs(pdbs, expected_order):
    sorted_member = sort_pdbs(pdbs)
    sorted_ids = [member.id for member in sorted_member]
    assert sorted_ids == expected_order


@pytest.mark.parametrize(
    "pdbs, expected_ids",
    [
        # First example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-250", 2.1),
                make_pdb("4DDD", "A=300-400", 8.1),
                make_pdb("5EEE", "A=300-400", 4.6),
                make_pdb("6FFF", "A=500-1000", 1.3),
                make_pdb("7GGG", "A=500-1000", 1.4),
                make_pdb("8HHH", "A=500-1000", 1.6),
            ],
            ["6FFF", "3CCC", "5EEE", "7GGG", "1AAA", "4DDD", "8HHH", "2BBB"],
            id="three_domains_split",
        ),
        # Second example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("4DDD", "A=200-400", 8.1),
                make_pdb("6FFF", "A=500-1000", 1.3),
                make_pdb("9III", "A=1-600", 4.2),
                make_pdb("10JJJ", "A=1-1000", 1.4),
            ],
            ["10JJJ", "6FFF", "1AAA", "9III", "4DDD"],
            id="overlap_merges",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=-", 1.4),
            ],
            ["1AAA", "2BBB"],
            id="invalid_chains_last",
        ),
    ],
)
def test_filter_pdbs_on_clustered_resolution(pdbs: list[PdbResult], expected_ids: set[str]):
    filtered_pdbs = filter_pdbs_on_clustered_resolution(pdbs, top=len(expected_ids))

    filtered_ids = [pdb.id for pdb in filtered_pdbs]
    assert filtered_ids == expected_ids
