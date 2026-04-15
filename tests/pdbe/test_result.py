import pytest

from protein_quest.pdbe.result import (
    PdbChainLengthError,
    PdbResult,
    filter_pdb_results_on_chain_length,
    filter_pdb_results_on_resolution,
)


@pytest.mark.parametrize(
    "query,expected",
    [
        ("O=1-300", "O"),  #  uniprot:A8MT69 pdb:7R5S
        ("B/D=1-81", "B"),  # uniprot:A8MT69 pdb:4E44
        (
            "B/D/H/L/M/N/U/V/W/X/Z/b/d/h/i/j/o/p/q/r=8-81",  # uniprot:A8MT69 pdb:4NE1
            "B",
        ),
        ("A/B=2-459,A/B=520-610", "A"),  # uniprot/O00255 pdb/3U84
        ("DD/Dd=1-1085", "DD"),  # uniprot/O00268 pdb/7ENA
        ("A=398-459,A=74-386,A=520-584,A=1-53", "A"),  # uniprot/O00255 pdb/7O9T
        ("A=-", "A"),  # uniprot/Q08499 pdb/1E9K
    ],
)
def test_pdbresult_chain(query, expected):
    pdb_result = PdbResult(id="DUMMY", method="DUMMY", uniprot_chains=query)
    result = pdb_result.chain

    assert result == expected


@pytest.mark.parametrize(
    "query,expected",
    [
        ("O=1-300", 300),  #  uniprot:A8MT69 pdb:7R5S
        ("B/D=1-81", 81),  # uniprot:A8MT69 pdb:4E44
        (
            "B/D/H/L/M/N/U/V/W/X/Z/b/d/h/i/j/o/p/q/r=8-81",  # uniprot:A8MT69 pdb:4NE1
            74,
        ),
        ("A/B=2-459,A/B=520-610", 549),  # uniprot/O00255 pdb/3U84
        ("DD/Dd=1-1085", 1085),  # uniprot/O00268 pdb/7ENA
        ("A=398-459,A=74-386,A=520-584,A=1-53", 493),  # uniprot/O00255 pdb/7O9T
    ],
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
