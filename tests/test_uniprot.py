from textwrap import dedent

import pytest

from protein_quest.uniprot import (
    PdbResult,
    Query,
    _append_subcellular_location_filters,
    _build_sparql_query_pdb,
    _build_sparql_query_uniprot,
    _first_chain_from_uniprot_chains,
    search4af,
    search4emdb,
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
        subcellular_location_go="GO:0005634",  # Cellular component - Nucleus
        molecular_function_go="GO:0003677",  # Molecular function - DNA binding
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
        subcellular_location_go="INVALID:0005634",  # Invalid GO term
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
    ],
)
def test_first_chain_from_uniprot_chains(query, expected):
    result = _first_chain_from_uniprot_chains(query)

    assert result == expected


@pytest.mark.vcr
def test_search4uniprot():
    query = Query(
        taxon_id="9606",
        reviewed=True,
        subcellular_location_uniprot="nucleus",
        subcellular_location_go="GO:0005634",  # Cellular component - Nucleus
        molecular_function_go="GO:0003677",  # Molecular function - DNA binding
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


@pytest.mark.vcr
def test_search4emdb():
    uniprot_accession = "P05067"
    results = search4emdb({uniprot_accession}, limit=1)

    expected = {uniprot_accession: {"EMD-0405"}}
    assert results == expected
