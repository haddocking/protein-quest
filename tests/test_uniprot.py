from textwrap import dedent

import pytest

from protein_quest.uniprot import (
    ComplexPortalEntry,
    PdbChainLengthError,
    PdbResult,
    Query,
    UniprotDetails,
    _append_subcellular_location_filters,
    _build_sparql_query_pdb,
    _build_sparql_query_uniprot,
    filter_pdb_results_on_chain_length,
    map_uniprot_accessions2uniprot_details,
    search4af,
    search4emdb,
    search4interaction_partners,
    search4macromolecular_complexes,
    search4pdb,
    search4uniprot,
)


def assertQueryEqual(actual, expected):
    """
    Helper function to assert that two SPARQL queries are equal.
    Strips leading whitespace for comparison.
    """
    actual_lines = [line.lstrip() for line in actual.split("\n")]
    expected_lines = [line.strip() for line in expected.split("\n")]
    assert actual_lines == expected_lines, f"Expected:\n{expected}\n\nActual:\n{actual}"


def test_build_sparql_query_uniprot():
    # Test with a simple query
    query = Query(
        taxon_id="9606",
        reviewed=True,
        subcellular_location_uniprot="nucleus",
        subcellular_location_go=["GO:0005634"],  # Cellular component - Nucleus
        molecular_function_go=["GO:0003677"],  # Molecular function - DNA binding
    )
    result = _build_sparql_query_uniprot(query, limit=10)

    expected = dedent("""
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX GO:<http://purl.obolibrary.org/obo/GO_>

        SELECT DISTINCT ?protein
        WHERE {

            # --- Protein Selection ---
            ?protein a up:Protein .
            ?protein up:organism taxon:9606 .
            ?protein up:reviewed true .

            {

            ?protein up:annotation ?subcellAnnotation .
            ?subcellAnnotation up:locatedIn/up:cellularComponent ?cellcmpt .
            ?cellcmpt skos:prefLabel "nucleus" .

            } UNION {

            ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005634 .

            }


            ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0003677 .


            }

        LIMIT 10
    """)

    assertQueryEqual(result, expected)


def test_build_sparql_query_pdb():
    result = _build_sparql_query_pdb(["O15178", "O15294"], limit=42)
    expected = dedent("""
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX GO:<http://purl.obolibrary.org/obo/GO_>

        SELECT ?protein ?pdb_db ?pdb_method ?pdb_resolution
        (GROUP_CONCAT(DISTINCT ?pdb_chain; separator=",") AS ?pdb_chains)

        WHERE {

            # --- Protein Selection ---
            VALUES (?ac) { ("O15178") ("O15294")}
            BIND (IRI(CONCAT("http://purl.uniprot.org/uniprot/",?ac)) AS ?protein)
            ?protein a up:Protein .


            # --- PDB Info ---
            ?protein rdfs:seeAlso ?pdb_db .
            ?pdb_db up:database <http://purl.uniprot.org/database/PDB> .
            ?pdb_db up:method ?pdb_method .
            ?pdb_db up:chainSequenceMapping ?chainSequenceMapping .
            BIND(STRAFTER(STR(?chainSequenceMapping), "isoforms/") AS ?isoformPart)
            FILTER(STRSTARTS(?isoformPart, CONCAT(?ac, "-")))
            ?chainSequenceMapping up:chain ?pdb_chain .
            OPTIONAL { ?pdb_db up:resolution ?pdb_resolution . }


        }
        GROUP BY ?protein ?pdb_db ?pdb_method ?pdb_resolution
        LIMIT 42
    """)
    assertQueryEqual(result, expected)


@pytest.mark.parametrize(
    "subcellular_location_uniprot,subcellular_location_go,expected",
    [
        # Test case 1: Neither filter provided
        (
            None,
            None,
            "",
        ),
        # Test case 2: Only UniProt subcellular location provided
        (
            "nucleus",
            None,
            dedent("""
            ?protein up:annotation ?subcellAnnotation .
            ?subcellAnnotation up:locatedIn/up:cellularComponent ?cellcmpt .
            ?cellcmpt skos:prefLabel "nucleus" .
        """),
        ),
        # Test case 3: Only single GO term provided
        (
            None,
            "GO:0005634",
            dedent("""
            ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005634 .
        """),
        ),
        # Test case 4: Only multiple GO terms provided (list)
        (
            None,
            ["GO:0005634", "GO:0005737"],
            "{ ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005634 . } UNION { ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005737 . }",
        ),
        # Test case 5: Only multiple GO terms provided (set)
        (
            None,
            {"GO:0005634", "GO:0005737"},
            "{ ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005634 . } UNION { ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005737 . }",
        ),
        # Test case 6: Both UniProt and single GO term provided
        (
            "nucleus",
            "GO:0005634",
            dedent("""
            {

            ?protein up:annotation ?subcellAnnotation .
            ?subcellAnnotation up:locatedIn/up:cellularComponent ?cellcmpt .
            ?cellcmpt skos:prefLabel "nucleus" .

            } UNION {

            ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005634 .

            }
        """),
        ),
        # Test case 7: Both UniProt and multiple GO terms provided
        (
            "cytoplasm",
            ["GO:0005634", "GO:0005737"],
            dedent("""
            {

            ?protein up:annotation ?subcellAnnotation .
            ?subcellAnnotation up:locatedIn/up:cellularComponent ?cellcmpt .
            ?cellcmpt skos:prefLabel "cytoplasm" .

            } UNION {
            { ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005634 . } UNION { ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) GO:0005737 . }
            }
        """),
        ),
    ],
)
def test_append_subcellular_location_filters(subcellular_location_uniprot, subcellular_location_go, expected):
    """Test _append_subcellular_location_filters with various input combinations."""
    query = Query(
        taxon_id=None,
        reviewed=None,
        subcellular_location_uniprot=subcellular_location_uniprot,
        subcellular_location_go=subcellular_location_go,
        molecular_function_go=None,
    )

    result = _append_subcellular_location_filters(query)

    # For sets, we need to handle the unordered nature
    if isinstance(subcellular_location_go, set) and len(subcellular_location_go) > 1:
        # Check that result contains both GO terms (in any order)
        assert "GO:0005634" in result
        assert "GO:0005737" in result
        assert "UNION" in result
        assert result.count("?protein up:classifiedWith") == 2
    else:
        # Normalize whitespace for comparison
        result_normalized = " ".join(result.split())
        expected_normalized = " ".join(expected.split())
        assert result_normalized == expected_normalized


def test_append_subcellular_location_filters_invalid_go_term():
    """Test that invalid GO terms raise ValueError."""
    query = Query(
        taxon_id=None,
        reviewed=None,
        subcellular_location_uniprot=None,
        subcellular_location_go=["INVALID:0005634"],  # Invalid GO term
        molecular_function_go=None,
    )

    with pytest.raises(ValueError, match="Subcellular location GO term must start with 'GO:'"):
        _append_subcellular_location_filters(query)


def test_append_subcellular_location_filters_invalid_go_term_in_list():
    """Test that invalid GO terms in a list raise ValueError."""
    query = Query(
        taxon_id=None,
        reviewed=None,
        subcellular_location_uniprot=None,
        subcellular_location_go=["GO:0005634", "INVALID:0005737"],  # One invalid GO term
        molecular_function_go=None,
    )

    with pytest.raises(ValueError, match="Subcellular location GO term must start with 'GO:'"):
        _append_subcellular_location_filters(query)


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
    pdb_result = PdbResult(id="DUMMY", method="DUMMY", uniprot_chains="A=-")

    with pytest.raises(PdbChainLengthError, match="Could not determine chain length of 'DUMMY' from 'A=-'"):
        _ = pdb_result.chain_length


def test_filter_pdb_results_on_chain_length_unchanged():
    pdbs = {
        "P05067": {PdbResult(id="1AAP", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A=287-344")},
    }
    result = filter_pdb_results_on_chain_length(pdbs, min_residues=None, max_residues=None)

    assert result is pdbs


def test_filter_pdb_results_on_chain_length_badrange():
    with pytest.raises(
        ValueError, match="Maximum number of residues \\(13\\) must be > minimum number of residues \\(42\\)"
    ):
        filter_pdb_results_on_chain_length({}, min_residues=42, max_residues=13)


def test_filter_pdb_results_on_chain_length_filtered():
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


@pytest.mark.vcr
def test_search4uniprot():
    query = Query(
        taxon_id="9606",
        reviewed=True,
        subcellular_location_uniprot="nucleus",
        subcellular_location_go=["GO:0005634"],  # Cellular component - Nucleus
        molecular_function_go=["GO:0003677"],  # Molecular function - DNA binding
    )

    results = search4uniprot(query, limit=1)

    expected = {"A0A087WUV0"}
    assert results == expected


@pytest.mark.vcr
def test_search4pdb():
    uniprot_accession = "P05067"

    results = search4pdb({uniprot_accession}, limit=1)

    expected = {
        uniprot_accession: {
            PdbResult(id="1AAP", method="X-Ray_Crystallography", resolution="1.5", uniprot_chains="A/B=287-344")
        }
    }
    assert results == expected
    assert next(iter(results[uniprot_accession])).chain == "A"


@pytest.mark.vcr
def test_search4af():
    uniprot_accession = "P05067"

    results = search4af({uniprot_accession}, limit=1)

    expected = {uniprot_accession: {uniprot_accession}}
    assert results == expected


# P05067 has a sequence length of 770 residues
@pytest.mark.vcr
def test_search4af_ok_sequence_length():
    uniprot_accession = "P05067"

    results = search4af({uniprot_accession}, limit=1, min_sequence_length=600, max_sequence_length=800)

    expected = {uniprot_accession: {uniprot_accession}}
    assert results == expected


@pytest.mark.vcr
def test_search4af_too_small_sequence_length():
    results = search4af({"P05067"}, limit=1, min_sequence_length=800)

    expected = {}
    assert results == expected


@pytest.mark.vcr
def test_search4af_too_big_sequence_length():
    results = search4af({"P05067"}, limit=1, max_sequence_length=600)

    expected = {}
    assert results == expected


def test_search4af_invalid_sequence_length():
    with pytest.raises(
        ValueError,
        match="Maximum sequence length \\(500\\) must be greater than minimum sequence length \\(600\\)",
    ):
        search4af({"P05067"}, limit=1, min_sequence_length=600, max_sequence_length=500)


class TestSearch4AfExternalIsoforms:
    # P42284 has P42284-2 as canonical isoform with 549 length
    # and several other isoforms based on external entries,
    # which should not be used for length filtering
    # one of them is Q7KQZ4-2 with length 787

    @pytest.mark.vcr
    def test_match_canonical_isoform(self):
        results = search4af({"P42284"}, min_sequence_length=540, max_sequence_length=560, limit=10)

        expected = {"P42284": {"P42284"}}
        assert results == expected

    @pytest.mark.vcr
    def test_do_not_match_external_isoform(self):
        # so setting min_sequence_length to 600 should exclude P42284
        # as only non-canonical isoforms are longer than that
        results = search4af({"P42284"}, min_sequence_length=600, limit=10)

        expected = {}
        assert results == expected


@pytest.mark.vcr
def test_search4emdb():
    uniprot_accession = "P05067"
    results = search4emdb({uniprot_accession}, limit=1)

    expected = {uniprot_accession: {"EMD-0405"}}
    assert results == expected


@pytest.mark.vcr
def test_search4macromolecular_complexes():
    uniprot_accession = "P60709"

    results = search4macromolecular_complexes({uniprot_accession}, limit=100)

    assert len(results) == 40
    first_expected = ComplexPortalEntry(
        complex_id="CPX-1203",
        complex_title="Brain-specific SWI/SNF ATP-dependent chromatin remodeling complex, ARID1A-SMARCA2 variant",
        complex_url="https://www.ebi.ac.uk/complexportal/complex/CPX-1203",
        members={
            "O94805",
            "P60709",
            "Q969G3",
            "P51531",
            "Q12824",
            "Q8TAQ2",
            "Q92925",
            "O14497",
        },
        query_protein="P60709",
    )
    first_result = results[0]
    assert first_result == first_expected


@pytest.mark.vcr
def test_search4interaction_partners():
    uniprot_accession = "P60709"
    excludes = {"Q92925", "O14497", "Q92922", "Q8TAQ2"}
    results = search4interaction_partners(uniprot_accession, excludes=excludes, limit=100)

    assert len(results) == 40
    expected_key = "O94805"
    first_expected = {
        "CPX-1203",
        "CPX-1210",
        "CPX-1220",
        "CPX-1218",
        "CPX-1221",
        "CPX-1217",
        "CPX-1196",
        "CPX-1209",
        "CPX-1207",
        "CPX-1211",
        "CPX-1228",
        "CPX-4224",
        "CPX-1227",
        "CPX-4223",
        "CPX-1216",
        "CPX-1202",
        "CPX-1219",
        "CPX-4225",
        "CPX-4226",
        "CPX-1226",
        "CPX-1225",
    }
    assert results[expected_key] == first_expected
    assert not results.keys() & excludes


@pytest.mark.vcr
def test_map_uniprot_accessions2uniprot_details():
    uniprot_accessions = ["P05067", "A6NGD5", "O14627", "P00697", "P42284", "A0A0B5AC95", "A0A0S2Z4R0"]
    results = set(map_uniprot_accessions2uniprot_details(uniprot_accessions))

    expected = {
        UniprotDetails(
            uniprot_accession="A0A0B5AC95",
            uniprot_id="INS1A_CONGE",
            sequence_length=115,
            reviewed=True,
            protein_name="Con-Ins G1a",
            taxon_id=6491,
            taxon_name="Conus geographus",
        ),
        UniprotDetails(
            uniprot_accession="A0A0S2Z4R0",
            uniprot_id="A0A0S2Z4R0_HUMAN",
            sequence_length=862,
            reviewed=False,
            protein_name="Axin-1",
            taxon_id=9606,
            taxon_name="Homo sapiens",
        ),
        UniprotDetails(
            uniprot_accession="A6NGD5",
            uniprot_id="ZSA5C_HUMAN",
            sequence_length=496,
            reviewed=True,
            protein_name="Zinc finger and SCAN domain-containing protein 5C",
            taxon_id=9606,
            taxon_name="Homo sapiens",
        ),
        UniprotDetails(
            uniprot_accession="O14627",
            uniprot_id="CDX4_HUMAN",
            sequence_length=284,
            reviewed=True,
            protein_name="Homeobox protein CDX-4",
            taxon_id=9606,
            taxon_name="Homo sapiens",
        ),
        UniprotDetails(
            uniprot_accession="P00697",
            uniprot_id="LYSC1_RAT",
            sequence_length=148,
            reviewed=True,
            protein_name="Lysozyme C-1",
            taxon_id=10116,
            taxon_name="Rattus norvegicus",
        ),
        UniprotDetails(
            uniprot_accession="P05067",
            uniprot_id="A4_HUMAN",
            sequence_length=770,
            reviewed=True,
            protein_name="Amyloid-beta precursor protein",
            taxon_id=9606,
            taxon_name="Homo sapiens",
        ),
        UniprotDetails(
            uniprot_accession="P42284",
            uniprot_id="LOLA2_DROME",
            sequence_length=549,
            reviewed=True,
            protein_name="Longitudinals lacking protein, isoforms H/M/V",
            taxon_id=7227,
            taxon_name="Drosophila melanogaster",
        ),
    }
    assert results == expected
