import pytest

from protein_quest.pdbe.result import (
    PdbChainLengthError,
    PdbResult,
    filter_pdb_results_on_chain_length,
    filter_pdb_results_on_resolution,
)

COMMON_CHAIN_CASES: list[tuple[str, str, int, int, int]] = [
    ("O=1-300", "O", 300, 1, 300),  # uniprot:A8MT69 pdb:7R5S
    ("B/D=1-81", "B", 81, 1, 81),  # uniprot:A8MT69 pdb:4E44
    (
        "B/D/H/L/M/N/U/V/W/X/Z/b/d/h/i/j/o/p/q/r=8-81",  # uniprot:A8MT69 pdb:4NE1
        "B",
        74,
        8,
        81,
    ),
    ("A/B=2-459,A/B=520-610", "A", 549, 2, 610),  # uniprot/O00255 pdb/3U84
    ("DD/Dd=1-1085", "DD", 1085, 1, 1085),  # uniprot/O00268 pdb/7ENA
    ("A=398-459,A=74-386,A=520-584,A=1-53", "A", 493, 1, 584),  # uniprot/O00255 pdb/7O9T
]
"""Format is
1. uniprot_chains string
2. expected chain (first chain in uniprot_chains)
3. expected chain length
4. expected uniprot start
5. expected uniprot end
"""


@pytest.mark.parametrize(
    "query,expected",
    [(query, chain) for query, chain, _length, _start, _end in COMMON_CHAIN_CASES]
    + [("A=-", "A")],  # uniprot/Q08499 pdb/1E9K
)
def test_pdbresult_chain(query, expected):
    pdb_result = PdbResult(id="DUMMY", method="DUMMY", uniprot_chains=query)
    result = pdb_result.chain

    assert result == expected


@pytest.mark.parametrize(
    "query,expected",
    [(query, length) for query, _chain, length, _start, _end in COMMON_CHAIN_CASES],
)
def test_pdb_result_chain_length(query, expected):
    pdb_result = PdbResult(id="DUMMY", method="DUMMY", uniprot_chains=query)
    result = pdb_result.chain_length

    assert result == expected


def test_pdb_result_chain_length_invalid():
    # uniprot/Q9NTW7 pdb/1X5W
    pdb_result = PdbResult(id="1X5W", method="NMR", uniprot_chains="A=-")

    with pytest.raises(PdbChainLengthError, match="Could not determine chain length of '1X5W' from 'A=-'"):
        _ = pdb_result.chain_length


@pytest.mark.parametrize(
    "query,expected_start,expected_end",
    [(query, start, end) for query, _chain, _length, start, end in COMMON_CHAIN_CASES],
)
def test_pdb_result_uniprot_range(query, expected_start, expected_end):
    pdb_result = PdbResult(id="DUMMY", method="DUMMY", uniprot_chains=query)

    assert pdb_result.uniprot_start == expected_start
    assert pdb_result.uniprot_end == expected_end


def test_pdb_result_uniprot_range_invalid_raises():
    pdb_result = PdbResult(id="1X5W", method="NMR", uniprot_chains="A=-")

    with pytest.raises(PdbChainLengthError, match="Could not determine chain length of '1X5W' from 'A=-'"):
        _ = pdb_result.uniprot_start

    with pytest.raises(PdbChainLengthError, match="Could not determine chain length of '1X5W' from 'A=-'"):
        _ = pdb_result.uniprot_end


class TestFilterPdbResultsOnChainLength:
    def test_unchanged(self):
        pdbs = {
            "P05067": {
                PdbResult(id="1AAP", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=287-344")
            },
        }
        result = filter_pdb_results_on_chain_length(pdbs, min_residues=None, max_residues=None)

        assert result is pdbs

    def test_badrange(self):
        with pytest.raises(
            ValueError, match="Maximum number of residues \\(13\\) must be > minimum number of residues \\(42\\)"
        ):
            filter_pdb_results_on_chain_length({}, min_residues=42, max_residues=13)

    def test_filtered(self):
        keeper = PdbResult(id="1AAP", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-100")
        pdbs = {
            "P05067": {
                keeper,
                PdbResult(id="2AAP", method="X-Ray_Crystallography", resolution="2.0", uniprot_chains="A=1-2000"),
            },
            "P12345": {
                PdbResult(id="3BBB", method="X-Ray_Crystallography", resolution="4.0", uniprot_chains="A=1-50"),
            },
        }
        result = filter_pdb_results_on_chain_length(pdbs, min_residues=75, max_residues=125)

        expected = {
            "P05067": {keeper},
        }
        assert result == expected


class TestFilterPdbResultsOnResolution:
    def test_keeps_top_per_accession(self):
        pdbs = {
            "P05067": {
                PdbResult(id="A", method="X-Ray_Crystallography", resolution="2.0", uniprot_chains="A=1-100"),
                PdbResult(id="B", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-80"),
                PdbResult(id="C", method="X-Ray_Crystallography", resolution="3.0", uniprot_chains="A=1-120"),
            },
            "Q8XYZ1": {
                PdbResult(id="D", method="X-Ray_Crystallography", resolution="4.0", uniprot_chains="A=1-100"),
                PdbResult(id="E", method="X-Ray_Crystallography", resolution="1.0", uniprot_chains="A=1-90"),
            },
        }

        result = filter_pdb_results_on_resolution(pdbs, top=1)

        expected = {
            "P05067": {PdbResult(id="B", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-80")},
            "Q8XYZ1": {PdbResult(id="E", method="X-Ray_Crystallography", resolution="1.0", uniprot_chains="A=1-90")},
        }
        assert result == expected

    def test_prefers_more_residues_on_tie(self):
        short = PdbResult(id="A", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-80")
        long = PdbResult(id="B", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-120")

        pdbs = {"P05067": {short, long}}

        result = filter_pdb_results_on_resolution(pdbs, top=1)

        assert result == {"P05067": {long}}

    def test_tie_breaks_on_pdb_id(self):
        a = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-120")
        b = PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-120")

        pdbs = {"P05067": {a, b}}

        result = filter_pdb_results_on_resolution(pdbs, top=1)

        assert result == {"P05067": {a}}

    def test_top_larger_than_group_keeps_all(self):
        a = PdbResult(id="A", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-120")
        b = PdbResult(id="B", method="X-Ray_Crystallography", resolution="2.5", uniprot_chains="A=1-80")

        pdbs = {"P05067": {a, b}}

        result = filter_pdb_results_on_resolution(pdbs, top=10)

        assert result == {"P05067": {a, b}}

    def test_deprioritizes_missing_resolution(self):
        missing_resolution = PdbResult(id="2K28", method="NMR_Spectroscopy", resolution=None, uniprot_chains="A=8-65")
        resolved = PdbResult(id="3I8Z", method="X-Ray_Crystallography", resolution="1.51", uniprot_chains="A=1-109")

        pdbs = {"O00257": {missing_resolution, resolved}}

        result = filter_pdb_results_on_resolution(pdbs, top=1)

        assert result == {"O00257": {resolved}}

    def test_deprioritizes_invalid_resolution_value(self):
        invalid_resolution = PdbResult(
            id="2BAD",
            method="X-Ray_Crystallography",
            resolution="not-a-number",
            uniprot_chains="A=1-120",
        )
        valid_resolution = PdbResult(
            id="1GOOD",
            method="X-Ray_Crystallography",
            resolution="1.5",
            uniprot_chains="A=1-80",
        )

        pdbs = {"P05067": {invalid_resolution, valid_resolution}}

        result = filter_pdb_results_on_resolution(pdbs, top=1)

        assert result == {"P05067": {valid_resolution}}

    def test_deprioritizes_invalid_chain_length_on_tie(self):
        invalid = PdbResult(id="1X5W", method="NMR_Spectroscopy", resolution="1.5", uniprot_chains="A=-")
        valid = PdbResult(id="2AAP", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=1-80")

        pdbs = {"Q9NTW7": {invalid, valid}}

        result = filter_pdb_results_on_resolution(pdbs, top=1)

        assert result == {"Q9NTW7": {valid}}

    @pytest.mark.parametrize("top", [0, -1])
    def test_nonpositive_top_raises(self, top):
        with pytest.raises(ValueError, match=f"Top must be a positive integer, got {top}"):
            filter_pdb_results_on_resolution({}, top=top)

    def test_on_sequence_domains(self):
        # First example from https://github.com/haddocking/protein-quest/issues/102
        # Given uniprot has length 1000
        m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
        m2 = PdbResult(id="2BBB", method="X-Ray_Crystallography", resolution="5.4", uniprot_chains="A=1-250")
        m3 = PdbResult(id="3CCC", method="X-Ray_Crystallography", resolution="2.1", uniprot_chains="A=1-250")
        m4 = PdbResult(id="4DDD", method="X-Ray_Crystallography", resolution="8.1", uniprot_chains="A=200-400")
        m5 = PdbResult(id="5EEE", method="X-Ray_Crystallography", resolution="4.6", uniprot_chains="A=200-400")
        m6 = PdbResult(id="6FFF", method="X-Ray_Crystallography", resolution="1.3", uniprot_chains="A=500-1000")
        m7 = PdbResult(id="7GGG", method="X-Ray_Crystallography", resolution="1.4", uniprot_chains="A=500-1000")
        m8 = PdbResult(id="8HHH", method="X-Ray_Crystallography", resolution="1.6", uniprot_chains="A=500-1000")
        pdbs = {
            "P12345": {
                m1,
                m2,
                m3,
                m4,
                m5,
                m6,
                m7,
                m8,
            }
        }
        results = filter_pdb_results_on_resolution(pdbs, top=3)

        expected = {"P12345": {m3, m5, m6}}
        assert results == expected

    def test_on_sequence_overlaps(self):
        # Second example from https://github.com/haddocking/protein-quest/issues/102
        m1 = PdbResult(id="1AAA", method="X-Ray_Crystallography", resolution="3.6", uniprot_chains="A=1-250")
        m4 = PdbResult(id="4DDD", method="X-Ray_Crystallography", resolution="8.1", uniprot_chains="A=200-400")
        m6 = PdbResult(id="6FFF", method="X-Ray_Crystallography", resolution="1.3", uniprot_chains="A=500-1000")
        m9 = PdbResult(id="9III", method="X-Ray_Crystallography", resolution="4.2", uniprot_chains="A=1-600")
        m10 = PdbResult(id="10JJJ", method="X-Ray_Crystallography", resolution="1.4", uniprot_chains="A=1-1000")

        pdbs = {
            "P12345": {
                m1,
                m4,
                m6,
                m9,
                m10,
            }
        }
        results = filter_pdb_results_on_resolution(pdbs, top=3)

        expected = {"P12345": {m9, m10, m6}}
        assert results == expected
