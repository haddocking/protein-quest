"""Helpers to construct and execute SPARQL queries against the UniProtKB SPARQL endpoint."""

import logging
from collections.abc import Iterable
from textwrap import dedent

from SPARQLWrapper import JSON, SPARQLWrapper

logger = logging.getLogger(__name__)


def build_sparql_generic_query(select_clause: str, where_clause: str, limit: int = 10_000, groupby_clause="") -> str:
    """
    Builds a generic SPARQL query with the given select and where clauses.

    Args:
        select_clause: The SELECT clause body of the SPARQL query, without the "SELECT" keyword.
        where_clause: The WHERE clause body of the SPARQL query, without the "WHERE" keyword.
        limit: The maximum number of results to return. Default is 10,000.
        groupby_clause: The GROUP BY clause of the SPARQL query, without the "GROUP BY" keyword.
    """
    groupby = f" GROUP BY {groupby_clause}" if groupby_clause else ""
    return dedent(f"""
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX GO:<http://purl.obolibrary.org/obo/GO_>

        SELECT {select_clause}
        WHERE {{
            {where_clause}
        }}
        {groupby}
        LIMIT {limit}
    """)


def build_sparql_generic_by_uniprot_accessions_query(
    uniprot_accs: Iterable[str], select_clause: str, where_clause: str, limit: int = 10_000, groupby_clause=""
) -> str:
    """Builds a generic SPARQL query and filters by the given UniProt accessions.

    Args:
        uniprot_accs: An iterable of UniProt accessions to filter by.
        select_clause: The SELECT clause of the SPARQL query, without the "SELECT" keyword.
        where_clause: The WHERE clause of the SPARQL query, without the "WHERE" keyword.
        limit: The maximum number of results to return. Default is 10,000.
        groupby_clause: The GROUP BY clause of the SPARQL query, without the "GROUP BY" keyword.

    Returns:
        A string containing the complete SPARQL query.
    """
    values = " ".join(f'("{ac}")' for ac in uniprot_accs)
    where_clause2 = dedent(f"""
        # --- Protein Selection ---
        VALUES (?ac) {{ {values}}}
        BIND (IRI(CONCAT("http://purl.uniprot.org/uniprot/",?ac)) AS ?protein)
        ?protein a up:Protein .

        {where_clause}
    """)
    return build_sparql_generic_query(
        select_clause=select_clause,
        where_clause=where_clause2,
        limit=limit,
        groupby_clause=groupby_clause,
    )


def execute_sparql_search(
    sparql_query: str,
    timeout: int,
) -> list:
    """
    Execute a SPARQL query.

    Args:
        sparql_query: The SPARQL query to execute.
        timeout: Timeout for the SPARQL query in seconds. Must be less than 2700 seconds (45 minutes).

    Returns:
        A list of bindings returned by the SPARQL query. Each binding is a dict mapping
    """
    if timeout > 2_700:
        msg = "Uniprot SPARQL timeout is limited to 2700 seconds (45 minutes)."
        raise ValueError(msg)

    # Execute the query
    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql")
    sparql.setReturnFormat(JSON)
    sparql.setTimeout(timeout)

    # Default is GET method which can be cached by the server so is preferred.
    # Too prevent URITooLong errors, we use POST method for large queries.
    too_long_for_get = 5_000
    if len(sparql_query) > too_long_for_get:
        sparql.setMethod("POST")

    sparql.setQuery(sparql_query)
    rawresults = sparql.queryAndConvert()
    if not isinstance(rawresults, dict):
        msg = f"Expected rawresults to be a dict, but got {type(rawresults)}"
        raise TypeError(msg)

    bindings = rawresults.get("results", {}).get("bindings")
    if not isinstance(bindings, list):
        logger.warning("SPARQL query did not return 'bindings' list as expected.")
        return []

    logger.debug(bindings)
    return bindings
