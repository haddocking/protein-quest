import pytest

from protein_quest.uniprot_chains import (
    UniprotChainRange,
    all_chain_ids,
    chain_length,
    covered_uniprot_range,
    first_chain_id,
    format_uniprot_chains,
    parse_chain_ids,
    parse_uniprot_chains,
    preferred_chain_id,
)


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        pytest.param(
            "A/B=2-459,A/B=520-610",
            (
                UniprotChainRange(chain_ids=("A", "B"), start=2, end=459),
                UniprotChainRange(chain_ids=("A", "B"), start=520, end=610),
            ),
            id="discontinuous-shared-chain-group",
        ),
        pytest.param(
            "A/B=268-443,C=268-373",
            (
                UniprotChainRange(chain_ids=("A", "B"), start=268, end=443),
                UniprotChainRange(chain_ids=("C",), start=268, end=373),
            ),
            id="split-over-multiple-chain-groups",
        ),
    ],
)
def test_parse_uniprot_chains(query: str, expected: tuple[UniprotChainRange, ...]):
    assert parse_uniprot_chains(query) == expected


@pytest.mark.parametrize(
    ("parsed", "expected"),
    [
        pytest.param(
            (
                UniprotChainRange(chain_ids=("A", "B"), start=2, end=459),
                UniprotChainRange(chain_ids=("C",), start=520, end=610),
            ),
            "A/B=2-459,C=520-610",
            id="multiple-groups",
        ),
    ],
)
def test_format_uniprot_chains(parsed: tuple[UniprotChainRange, ...], expected: str):
    assert format_uniprot_chains(parsed) == expected


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        pytest.param("1/B/G=1-124", "B", id="skip-leading-numeric"),
        pytest.param("A/B=268-443,C=268-373", "A", id="first-chain-of-first-group"),
    ],
)
def test_first_chain_id(query: str, expected: str):
    assert first_chain_id(parse_uniprot_chains(query)) == expected


@pytest.mark.parametrize(
    ("query", "expected_chain_ids", "expected_preferred_chain_id"),
    [
        pytest.param("A=-", ("A",), "A", id="single-chain-missing-range"),
        pytest.param("1/B=-", ("1", "B"), "B", id="numeric-leader-missing-range"),
    ],
)
def test_parse_chain_ids(query: str, expected_chain_ids: tuple[str, ...], expected_preferred_chain_id: str):
    chain_ids = parse_chain_ids(query)

    assert chain_ids == expected_chain_ids
    assert preferred_chain_id(chain_ids) == expected_preferred_chain_id


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        pytest.param("A/B=2-459,A/B=520-610,C=700-710", ("A", "B", "C"), id="deduplicate-in-order"),
        pytest.param("A/B=268-443,C=268-373", ("A", "B", "C"), id="multiple-groups"),
    ],
)
def test_all_chain_ids(query: str, expected: tuple[str, ...]):
    assert all_chain_ids(parse_uniprot_chains(query)) == expected


@pytest.mark.parametrize(
    ("query", "expected_length", "expected_range"),
    [
        pytest.param("A/B=2-459,A/B=520-610", 549, (2, 610), id="discontinuous-shared-chain-group"),
        pytest.param("A/B=268-443,C=268-373", 282, (268, 443), id="multiple-chain-groups"),
    ],
)
def test_chain_length_and_covered_range(query: str, expected_length: int, expected_range: tuple[int, int]):
    parsed = parse_uniprot_chains(query)

    assert chain_length(parsed) == expected_length
    assert covered_uniprot_range(parsed) == expected_range


@pytest.mark.parametrize("query", ["A=-", "", "A=10-9"])
def test_parse_uniprot_chains_invalid(query: str):
    with pytest.raises(ValueError):
        parse_uniprot_chains(query)
