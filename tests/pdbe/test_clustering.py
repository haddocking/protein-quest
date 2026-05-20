import pytest

from protein_quest.pdbe.clustering import cluster_pdbs, pdb_distance
from protein_quest.pdbe.result import PdbResult


def make_pdb(pdb_id: str, uniprot_chains: str, resolution: float = 5.0) -> PdbResult:
    return PdbResult(
        id=pdb_id, method="X-Ray_Crystallography", resolution=str(resolution), uniprot_chains=uniprot_chains
    )


def case_three_domains() -> list[PdbResult]:
    return [
        make_pdb("1AAA", "A=1-250", 3.6),
        make_pdb("2BBB", "A=1-250", 5.4),
        make_pdb("3CCC", "A=1-250", 2.1),
        make_pdb("4DDD", "A=300-400", 8.1),
        make_pdb("5EEE", "A=300-400", 4.6),
        make_pdb("6FFF", "A=500-1000", 1.3),
        make_pdb("7GGG", "A=500-1000", 1.4),
        make_pdb("8HHH", "A=500-1000", 1.6),
    ]


def case_nested_ranges() -> list[PdbResult]:
    return [
        make_pdb("1AAA", "A=1-250", 3.6),
        make_pdb("2BBB", "A=2-249", 5.4),
        make_pdb("3CCC", "A=3-248", 2.1),
    ]


def case_split_range_overlap() -> list[PdbResult]:
    return [
        make_pdb("1AAA", "A=1-250", 3.6),
        make_pdb("2BBB", "A=1-250", 5.4),
        make_pdb("3CCC", "A=1-54,A=90-250", 2.1),
    ]


def case_overlap_merges() -> list[PdbResult]:
    # # Second example from https://github.com/haddocking/protein-quest/issues/102
    return [
        make_pdb("1AAA", "A=1-250", 3.6),
        make_pdb("4DDD", "A=200-400", 8.1),
        make_pdb("6FFF", "A=500-1000", 1.3),
        make_pdb("9III", "A=1-600", 4.2),
        make_pdb("10JJJ", "A=1-1000", 1.4),
    ]


@pytest.mark.parametrize(
    "a,b,expected",
    [
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("2BBB", "A=1-250", 5.4),
            1 / 250,
            id="identical_ranges",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("3CCC", "A=1-54", 2.1),
            1 / 54,
            id="partial_overlap",
        ),
        pytest.param(
            make_pdb("1AAA", "A=1-250", 3.6),
            make_pdb("4DDD", "A=300-400", 8.1),
            float("+inf"),
            id="no_overlap",
        ),
    ],
)
def test_pdb_distance(a, b, expected):
    assert pdb_distance(a, b) == expected


@pytest.mark.parametrize(
    "pdbs, expected_member_ids",
    [
        pytest.param([], [], id="empty_input"),
        pytest.param(
            [make_pdb("1AAA", "A=1-250", 3.6)],
            [{"1AAA"}],
            id="single_valid_shortcut",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-250", 2.1),
            ],
            [{"1AAA", "2BBB", "3CCC"}],
            id="same_range_merges",
        ),
        pytest.param(
            case_three_domains(),
            [{"1AAA", "2BBB", "3CCC"}, {"4DDD", "5EEE"}, {"6FFF", "7GGG", "8HHH"}],
            id="three_domains_split",
        ),
        pytest.param(case_split_range_overlap(), [{"1AAA", "2BBB", "3CCC"}], id="split_range_merges"),
        pytest.param(
            [
                make_pdb("1BAD", "A=-", 1.0),
                make_pdb("2BAD", "A=-", 2.0),
            ],
            [{"1BAD", "2BAD"}],
            id="all_invalid_shortcut",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("1BAD", "A=-", 2.0),
            ],
            [{"1AAA"}, {"1BAD"}],
            id="mixed_valid_invalid_separates",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("2BBB", "A=100-200", 2.0),
            ],
            [{"1AAA", "2BBB"}],
            id="boundary_overlap_merges",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("2BBB", "A=101-200", 2.0),
            ],
            [{"1AAA"}, {"2BBB"}],
            id="adjacent_no_overlap_separates",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-100", 1.0),
                make_pdb("2BBB", "A=50-150", 2.0),
                make_pdb("3CCC", "A=120-200", 3.0),
            ],
            [{"1AAA", "2BBB"}, {"3CCC"}],
            id="complete_linkage_bridge_splits",
        ),
        pytest.param(case_nested_ranges(), [{"1AAA", "2BBB", "3CCC"}], id="nested_ranges_merge"),
        # TODO check clusters of case_overlap are correct
        pytest.param(case_overlap_merges(), [{"1AAA", "4DDD", "9III", "10JJJ"}, {"6FFF"}], id="overlap_merges"),
    ],
)
def test_cluster_pdbs(pdbs, expected_member_ids):
    clusters = cluster_pdbs(pdbs)

    cluster_member_ids = {frozenset(p.id for p in cluster) for cluster in clusters}
    expected = {frozenset(ids) for ids in expected_member_ids}
    assert cluster_member_ids == expected
