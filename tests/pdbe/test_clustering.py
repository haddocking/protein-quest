from dataclasses import dataclass

import pytest

from protein_quest.pdbe.result import PdbResult


def case_1pdb() -> list[PdbResult]:
    m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
    return [m1]


def case_1domain_3pdbs_nooverlap_samerange() -> list[PdbResult]:
    m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
    m2 = PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="5.4", uniprot_chains="A=1-250")
    m3 = PdbResult(id="3CCC", method="X-Ray_Crystallography", resolution="2.1", uniprot_chains="A=1-250")
    return [m1, m2, m3]


def case_1domain_3pdbs_nooverlap_withnonperfectsequenceidentity() -> list[PdbResult]:
    m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
    m2 = PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="5.4", uniprot_chains="A=1-250")
    # m3 has non-perfect sequence identity due to gaps, the sequence identity must be ignored during clustering
    m3 = PdbResult(id="3CCC", method="X-Ray_Crystallography", resolution="2.1", uniprot_chains="A=1-54,A=90-250")
    return [m1, m2, m3]


def case_1domain_3pdbs_nested() -> list[PdbResult]:
    m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
    m2 = PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="5.4", uniprot_chains="A=2-249")
    m3 = PdbResult(id="3CCC", method="X-Ray_Crystallography", resolution="2.1", uniprot_chains="A=3-248")
    return [m1, m2, m3]


def case_3domains_3pdbs_nooverlap_samerange() -> list[PdbResult]:
    m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
    m2 = PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="5.4", uniprot_chains="A=1-250")
    m3 = PdbResult(id="3CCC", method="X-Ray_Crystallography", resolution="2.1", uniprot_chains="A=1-250")
    m4 = PdbResult(id="4DDD", method="X-Ray_Crystallography", resolution="8.1", uniprot_chains="A=300-400")
    m5 = PdbResult(id="5EEE", method="X-Ray_Crystallography", resolution="4.6", uniprot_chains="A=300-400")
    m6 = PdbResult(id="6FFF", method="X-Ray_Crystallography", resolution="1.3", uniprot_chains="A=500-1000")
    m7 = PdbResult(id="7GGG", method="X-Ray_Crystallography", resolution="1.4", uniprot_chains="A=500-1000")
    m8 = PdbResult(id="8HHH", method="X-Ray_Crystallography", resolution="1.6", uniprot_chains="A=500-1000")
    return [m1, m2, m3, m4, m5, m6, m7, m8]


def pdb_distance(a: PdbResult, b: PdbResult) -> float:
    """Counts the number of residues in command a and b that overlap

    and b are considered more similar if they have more overlapping residues, and less similar if they have fewer overlapping residues.
    The distance is defined as the inverse of the number of overlapping residues, so that a higher number of overlapping residues results in a smaller distance.
    """
    a_range = set(range(a.uniprot_start, a.uniprot_end + 1))
    b_range = set(range(b.uniprot_start, b.uniprot_end + 1))
    overlap = a_range.intersection(b_range)
    if not overlap:
        return float("+inf")  # No overlap, infinite distance
    return 1 / len(overlap)


@pytest.mark.parametrize(
    "a,b,expected",
    [
        (
            PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250"),
            PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="5.4", uniprot_chains="A=1-250"),
            1 / 250,
        ),
        (
            PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250"),
            PdbResult(id="3CCC", method="X-Ray_Crystallography", resolution="2.1", uniprot_chains="A=1-54"),
            1 / 54,
        ),
        (
            PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250"),
            PdbResult(id="4DDD", method="X-Ray_Crystallography", resolution="8.1", uniprot_chains="A=300-400"),
            float("+inf"),
        ),
    ],
)
def test_pdb_distance(a, b, expected):
    assert pdb_distance(a, b) == expected


@dataclass
class ClusterPdbOptions:
    pass


def cluster_pdbs(pdbs: list[PdbResult], options: ClusterPdbOptions | None = None) -> list[set[PdbResult]]:
    if options is None:
        options = ClusterPdbOptions()

    buckets: dict[tuple[int, int], set[PdbResult]] = {}
    for pdb in pdbs:
        interval = (pdb.uniprot_start, pdb.uniprot_end)
        buckets.setdefault(interval, set()).add(pdb)

    return [buckets[interval] for interval in sorted(buckets)]


@pytest.mark.parametrize(
    "pdbs, expected_member_ids",
    [
        (case_1pdb(), [{"1AAA"}]),
        (case_1domain_3pdbs_nooverlap_samerange(), [{"1AAA", "2BBB", "3CCC"}]),
        (
            case_3domains_3pdbs_nooverlap_samerange(),
            [{"1AAA", "2BBB", "3CCC"}, {"4DDD", "5EEE"}, {"6FFF", "7GGG", "8HHH"}],
        ),
        (case_1domain_3pdbs_nooverlap_withnonperfectsequenceidentity(), [{"1AAA", "2BBB", "3CCC"}]),
        (case_1domain_3pdbs_nested(), [{"1AAA", "2BBB", "3CCC"}]),
    ],
)
def test_cluster_pdbs(pdbs, expected_member_ids):
    clusters: list[set[PdbResult]] = cluster_pdbs(pdbs)

    cluster_member_ids = {frozenset(p.id for p in cluster) for cluster in clusters}
    expected = {frozenset(ids) for ids in expected_member_ids}
    assert cluster_member_ids == expected
