"""Module for searching UniProtKB using SPARQL."""

import logging
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from itertools import batched
from textwrap import dedent

from SPARQLWrapper import JSON, SPARQLWrapper
from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


@dataclass
class Query:
    """Search query for UniProtKB.

    Parameters:
        taxon_id: NCBI Taxon ID to filter results by organism (e.g., "9606" for human).
        reviewed: Whether to filter results by reviewed status (True for reviewed, False for unreviewed).
        subcellular_location_uniprot: Subcellular location in UniProt format (e.g., "nucleus").
        subcellular_location_go: Subcellular location in GO format. Can be a single GO term
            (e.g., ["GO:0005634"]) or a collection of GO terms (e.g., ["GO:0005634", "GO:0005737"]).
        molecular_function_go: Molecular function in GO format. Can be a single GO term
            (e.g., ["GO:0003674"]) or a collection of GO terms (e.g., ["GO:0003674", "GO:0008150"]).
    """

    # TODO make taxon_id an int
    taxon_id: str | None
    reviewed: bool | None = None
    subcellular_location_uniprot: str | None = None
    subcellular_location_go: list[str] | None = None
    molecular_function_go: list[str] | None = None


def _first_chain_from_uniprot_chains(uniprot_chains: str) -> str:
    """Extracts the first chain identifier from a UniProt chains string.

    The UniProt chains string is formatted (with EBNF notation) as follows:

        chain_group(=range)?(,chain_group(=range)?)*

    where:
        chain_group := chain_id(/chain_id)*
        chain_id    := [A-Za-z]+
        range       := start-end
        start, end  := integer

    Args:
        uniprot_chains: A string representing UniProt chains, For example "B/D=1-81".
    Returns:
        The first chain identifier from the UniProt chain string. For example "B".
    """
    chains = uniprot_chains.split("=")
    parts = chains[0].split("/")
    chain = parts[0]
    try:
        # Workaround for Q9Y2Q5 │ 5YK3 │ 1/B/G=1-124, 1 does not exist but B does
        int(chain)
        if len(parts) > 1:
            return parts[1]
    except ValueError:
        # A letter
        pass
    return chain


@dataclass(frozen=True)
class PdbResult:
    """Result of a PDB search in UniProtKB.

    Parameters:
        id: PDB ID (e.g., "1H3O").
        method: Method used for the PDB entry (e.g., "X-ray diffraction").
        uniprot_chains: Chains in UniProt format (e.g., "A/B=1-42,A/B=50-99").
        resolution: Resolution of the PDB entry (e.g., "2.0" for 2.0 Å). Optional.
    """

    id: str
    method: str
    uniprot_chains: str
    resolution: str | None = None

    @property
    def chain(self) -> str:
        """The first chain from the UniProt chains aka self.uniprot_chains."""
        return _first_chain_from_uniprot_chains(self.uniprot_chains)


def _query2dynamic_sparql_triples(query: Query):
    parts: list[str] = []
    if query.taxon_id:
        parts.append(f"?protein up:organism taxon:{query.taxon_id} .")

    if query.reviewed:
        parts.append("?protein up:reviewed true .")
    elif query.reviewed is False:
        parts.append("?protein up:reviewed false .")

    parts.append(_append_subcellular_location_filters(query))

    if query.molecular_function_go:
        # Handle both single GO term (string) and multiple GO terms (list)
        if isinstance(query.molecular_function_go, str):
            go_terms = [query.molecular_function_go]
        else:
            go_terms = query.molecular_function_go

        molecular_function_filter = _create_go_filter(go_terms, "Molecular function")
        parts.append(molecular_function_filter)

    return "\n".join(parts)


def _create_go_filter(go_terms: Collection[str], term_type: str) -> str:
    """Create SPARQL filter for GO terms.

    Args:
        go_terms: Collection of GO terms to filter by.
        term_type: Type of GO terms for error messages (e.g., "Molecular function", "Subcellular location").

    Returns:
        SPARQL filter string.
    """
    # Validate all GO terms start with "GO:"
    for term in go_terms:
        if not term.startswith("GO:"):
            msg = f"{term_type} GO term must start with 'GO:', got: {term}"
            raise ValueError(msg)

    if len(go_terms) == 1:
        # Single GO term - get the first (and only) term
        term = next(iter(go_terms))
        return dedent(f"""
            ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) {term} .
        """)

    # Multiple GO terms - use UNION for OR logic
    union_parts = [
        dedent(f"""
            {{ ?protein up:classifiedWith|(up:classifiedWith/rdfs:subClassOf+) {term} . }}
        """).strip()
        for term in go_terms
    ]
    return " UNION ".join(union_parts)


def _append_subcellular_location_filters(query: Query) -> str:
    subcellular_location_uniprot_part = ""
    subcellular_location_go_part = ""

    if query.subcellular_location_uniprot:
        subcellular_location_uniprot_part = dedent(f"""
            ?protein up:annotation ?subcellAnnotation .
            ?subcellAnnotation up:locatedIn/up:cellularComponent ?cellcmpt .
            ?cellcmpt skos:prefLabel "{query.subcellular_location_uniprot}" .
        """)

    if query.subcellular_location_go:
        # Handle both single GO term (string) and multiple GO terms (list)
        if isinstance(query.subcellular_location_go, str):
            go_terms = [query.subcellular_location_go]
        else:
            go_terms = query.subcellular_location_go

        subcellular_location_go_part = _create_go_filter(go_terms, "Subcellular location")

    if subcellular_location_uniprot_part and subcellular_location_go_part:
        # If both are provided include results for both with logical OR
        return dedent(f"""
            {{
                {subcellular_location_uniprot_part}
            }} UNION {{
                {subcellular_location_go_part}
            }}
        """)

    return subcellular_location_uniprot_part or subcellular_location_go_part


def _build_sparql_generic_query(select_clause: str, where_clause: str, limit: int = 10_000, groupby_clause="") -> str:
    """
    Builds a generic SPARQL query with the given select and where clauses.
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


def _build_sparql_generic_by_uniprot_accessions_query(
    uniprot_accs: Iterable[str], select_clause: str, where_clause: str, limit: int = 10_000, groupby_clause=""
) -> str:
    values = " ".join(f'("{ac}")' for ac in uniprot_accs)
    where_clause2 = dedent(f"""
        # --- Protein Selection ---
        VALUES (?ac) {{ {values}}}
        BIND (IRI(CONCAT("http://purl.uniprot.org/uniprot/",?ac)) AS ?protein)
        ?protein a up:Protein .

        {where_clause}
    """)
    return _build_sparql_generic_query(
        select_clause=select_clause,
        where_clause=where_clause2,
        limit=limit,
        groupby_clause=groupby_clause,
    )


def _build_sparql_query_uniprot(query: Query, limit=10_000) -> str:
    dynamic_triples = _query2dynamic_sparql_triples(query)
    # TODO add usefull columns that have 1:1 mapping to protein
    # like uniprot_id with `?protein up:mnemonic ?mnemonic .`
    # and sequence, take care to take first isoform
    # ?protein up:sequence ?isoform .
    # ?isoform rdf:value ?sequence .
    select_clause = "DISTINCT ?protein"
    where_clause = dedent(f"""
        # --- Protein Selection ---
        ?protein a up:Protein .
        {dynamic_triples}
    """)
    return _build_sparql_generic_query(select_clause, dedent(where_clause), limit)


def _build_sparql_query_pdb(uniprot_accs: Iterable[str], limit=10_000) -> str:
    # For http://purl.uniprot.org/uniprot/O00268 + http://rdf.wwpdb.org/pdb/1H3O
    # the chainSequenceMapping are
    # http://purl.uniprot.org/isoforms/O00268-1#PDB_1H3O_tt872tt945
    # http://purl.uniprot.org/isoforms/Q16514-1#PDB_1H3O_tt57tt128
    # For http://purl.uniprot.org/uniprot/O00255 + http://rdf.wwpdb.org/pdb/3U84
    # the chainSequenceMapping are
    # http://purl.uniprot.org/isoforms/O00255-2#PDB_3U84_tt520tt610
    # http://purl.uniprot.org/isoforms/O00255-2#PDB_3U84_tt2tt459
    # To get the the chain belonging to the uniprot/pdb pair we need to
    # do some string filtering.
    # Also there can be multiple cnhins for the same uniprot/pdb pair, so we need to
    # do a group by and concat

    select_clause = dedent("""\
        ?protein ?pdb_db ?pdb_method ?pdb_resolution
         (GROUP_CONCAT(DISTINCT ?pdb_chain; separator=",") AS ?pdb_chains)
    """)

    where_clause = dedent("""
        # --- PDB Info ---
        ?protein rdfs:seeAlso ?pdb_db .
        ?pdb_db up:database <http://purl.uniprot.org/database/PDB> .
        ?pdb_db up:method ?pdb_method .
        ?pdb_db up:chainSequenceMapping ?chainSequenceMapping .
        BIND(STRAFTER(STR(?chainSequenceMapping), "isoforms/") AS ?isoformPart)
        FILTER(STRSTARTS(?isoformPart, CONCAT(?ac, "-")))
        ?chainSequenceMapping up:chain ?pdb_chain .
        OPTIONAL { ?pdb_db up:resolution ?pdb_resolution . }
    """)

    groupby_clause = "?protein ?pdb_db ?pdb_method ?pdb_resolution"
    return _build_sparql_generic_by_uniprot_accessions_query(
        uniprot_accs, select_clause, where_clause, limit, groupby_clause
    )


def _build_sparql_query_af(uniprot_accs: Iterable[str], limit=10_000) -> str:
    select_clause = "?protein ?af_db"
    where_clause = dedent("""
        # --- Protein Selection ---
        ?protein a up:Protein .

        # --- AlphaFoldDB Info ---
        ?protein rdfs:seeAlso ?af_db .
        ?af_db up:database <http://purl.uniprot.org/database/AlphaFoldDB> .
    """)
    return _build_sparql_generic_by_uniprot_accessions_query(uniprot_accs, select_clause, dedent(where_clause), limit)


def _build_sparql_query_emdb(uniprot_accs: Iterable[str], limit=10_000) -> str:
    select_clause = "?protein ?emdb_db"
    where_clause = dedent("""
        # --- Protein Selection ---
        ?protein a up:Protein .

        # --- EMDB Info ---
        ?protein rdfs:seeAlso ?emdb_db .
        ?emdb_db up:database <http://purl.uniprot.org/database/EMDB> .
    """)
    return _build_sparql_generic_by_uniprot_accessions_query(uniprot_accs, select_clause, dedent(where_clause), limit)


def _execute_sparql_search(
    sparql_query: str,
    timeout: int,
) -> list:
    """
    Execute a SPARQL query.
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


def _flatten_results_pdb(rawresults: Iterable) -> dict[str, set[PdbResult]]:
    pdb_entries: dict[str, set[PdbResult]] = {}
    for result in rawresults:
        protein = result["protein"]["value"].split("/")[-1]
        if "pdb_db" not in result:  # Should not happen with build_sparql_query_pdb
            continue
        pdb_id = result["pdb_db"]["value"].split("/")[-1]
        method = result["pdb_method"]["value"].split("/")[-1]
        uniprot_chains = result["pdb_chains"]["value"]
        pdb = PdbResult(id=pdb_id, method=method, uniprot_chains=uniprot_chains)
        if "pdb_resolution" in result:
            pdb = PdbResult(
                id=pdb_id,
                method=method,
                uniprot_chains=uniprot_chains,
                resolution=result["pdb_resolution"]["value"],
            )
        if protein not in pdb_entries:
            pdb_entries[protein] = set()
        pdb_entries[protein].add(pdb)

    return pdb_entries


def _flatten_results_af(rawresults: Iterable) -> dict[str, set[str]]:
    alphafold_entries: dict[str, set[str]] = {}
    for result in rawresults:
        protein = result["protein"]["value"].split("/")[-1]
        if "af_db" in result:
            af_id = result["af_db"]["value"].split("/")[-1]
            if protein not in alphafold_entries:
                alphafold_entries[protein] = set()
            alphafold_entries[protein].add(af_id)
    return alphafold_entries


def _flatten_results_emdb(rawresults: Iterable) -> dict[str, set[str]]:
    emdb_entries: dict[str, set[str]] = {}
    for result in rawresults:
        protein = result["protein"]["value"].split("/")[-1]
        if "emdb_db" in result:
            emdb_id = result["emdb_db"]["value"].split("/")[-1]
            if protein not in emdb_entries:
                emdb_entries[protein] = set()
            emdb_entries[protein].add(emdb_id)
    return emdb_entries


def limit_check(what: str, limit: int, len_raw_results: int):
    if len_raw_results >= limit:
        logger.warning(
            "%s returned %d results. "
            "There may be more results available, "
            "but they are not returned due to the limit of %d. "
            "Consider increasing the limit to get more results.",
            what,
            len_raw_results,
            limit,
        )


def search4uniprot(query: Query, limit: int = 10_000, timeout: int = 1_800) -> set[str]:
    """
    Search for UniProtKB entries based on the given query.

    Args:
        query: Query object containing search parameters.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.

    Returns:
        Set of uniprot accessions.
    """
    sparql_query = _build_sparql_query_uniprot(query, limit)
    logger.info("Executing SPARQL query for UniProt: %s", sparql_query)

    # Type assertion is needed because _execute_sparql_search returns a Union
    raw_results = _execute_sparql_search(
        sparql_query=sparql_query,
        timeout=timeout,
    )
    limit_check("Search for uniprot accessions", limit, len(raw_results))
    return {result["protein"]["value"].split("/")[-1] for result in raw_results}


def search4pdb(
    uniprot_accs: Collection[str], limit: int = 10_000, timeout: int = 1_800, batch_size: int = 10_000
) -> dict[str, set[PdbResult]]:
    """
    Search for PDB entries in UniProtKB accessions.

    Args:
        uniprot_accs: UniProt accessions.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.
        batch_size: Size of batches to process the UniProt accessions.

    Returns:
        Dictionary with protein IDs as keys and sets of PDB results as values.
    """
    all_raw_results = []
    total = len(uniprot_accs)
    with tqdm(total=total, desc="Searching for PDBs of uniprots", disable=total < batch_size, unit="acc") as pbar:
        for batch in batched(uniprot_accs, batch_size, strict=False):
            sparql_query = _build_sparql_query_pdb(batch, limit)
            logger.info("Executing SPARQL query for PDB: %s", sparql_query)

            raw_results = _execute_sparql_search(
                sparql_query=sparql_query,
                timeout=timeout,
            )
            all_raw_results.extend(raw_results)
            pbar.update(len(batch))

    limit_check("Search for pdbs on uniprot", limit, len(all_raw_results))
    return _flatten_results_pdb(all_raw_results)


def search4af(
    uniprot_accs: Collection[str], limit: int = 10_000, timeout: int = 1_800, batch_size: int = 10_000
) -> dict[str, set[str]]:
    """
    Search for AlphaFold entries in UniProtKB accessions.

    Args:
        uniprot_accs: UniProt accessions.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.
        batch_size: Size of batches to process the UniProt accessions.

    Returns:
        Dictionary with protein IDs as keys and sets of AlphaFold IDs as values.
    """
    all_raw_results = []
    total = len(uniprot_accs)
    with tqdm(total=total, desc="Searching for AlphaFolds of uniprots", disable=total < batch_size, unit="acc") as pbar:
        for batch in batched(uniprot_accs, batch_size, strict=False):
            sparql_query = _build_sparql_query_af(batch, limit)
            logger.info("Executing SPARQL query for AlphaFold: %s", sparql_query)

            raw_results = _execute_sparql_search(
                sparql_query=sparql_query,
                timeout=timeout,
            )
            all_raw_results.extend(raw_results)
            pbar.update(len(batch))

    limit_check("Search for alphafold entries on uniprot", limit, len(all_raw_results))
    return _flatten_results_af(all_raw_results)


def search4emdb(uniprot_accs: Iterable[str], limit: int = 10_000, timeout: int = 1_800) -> dict[str, set[str]]:
    """
    Search for EMDB entries in UniProtKB accessions.

    Args:
        uniprot_accs: UniProt accessions.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.

    Returns:
        Dictionary with protein IDs as keys and sets of EMDB IDs as values.
    """
    sparql_query = _build_sparql_query_emdb(uniprot_accs, limit)
    logger.info("Executing SPARQL query for EMDB: %s", sparql_query)

    raw_results = _execute_sparql_search(
        sparql_query=sparql_query,
        timeout=timeout,
    )
    limit_check("Search for EMDB entries on uniprot", limit, len(raw_results))
    return _flatten_results_emdb(raw_results)


def _build_complex_sparql_query(uniprot_accs: Iterable[str], limit: int) -> str:
    """Builds a SPARQL query to retrieve ComplexPortal information for given UniProt accessions.

    Example:

    ```sparql
    PREFIX up:   <http://purl.uniprot.org/core/>
    PREFIX rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT
    ?protein
    ?cp_db
    ?cp_comment
    (GROUP_CONCAT(DISTINCT ?member; separator=",") AS ?complex_members)
    (COUNT(DISTINCT ?member) AS ?member_count)
    WHERE {
    # Input UniProt accessions
    VALUES (?ac) { ("P05067") ("P60709") ("Q05471")}
    BIND (IRI(CONCAT("http://purl.uniprot.org/uniprot/", ?ac)) AS ?protein)

    # ComplexPortal cross-reference for each input protein
    ?protein a up:Protein ;
            rdfs:seeAlso ?cp_db .
    ?cp_db up:database <http://purl.uniprot.org/database/ComplexPortal> .
    OPTIONAL { ?cp_db rdfs:comment ?cp_comment . }

    # All member proteins of the same ComplexPortal complex
    ?member a up:Protein ;
            rdfs:seeAlso ?cp_db .
    }
    GROUP BY ?protein ?cp_db ?cp_comment
    ORDER BY ?protein ?cp_db
    LIMIT 500
    ```

    """
    select_clause = dedent("""\
        ?protein ?cp_db ?cp_comment
        (GROUP_CONCAT(DISTINCT ?member; separator=",") AS ?complex_members)
    """)
    where_clause = dedent("""
        # --- Complex Info ---
        ?protein a up:Protein ;
                rdfs:seeAlso ?cp_db .
        ?cp_db up:database <http://purl.uniprot.org/database/ComplexPortal> .
        OPTIONAL { ?cp_db rdfs:comment ?cp_comment . }
        # All member proteins of the same ComplexPortal complex
        ?member a up:Protein ;
        rdfs:seeAlso ?cp_db .
    """)
    group_by = dedent("""
       ?protein ?cp_db ?cp_comment
    """)
    return _build_sparql_generic_by_uniprot_accessions_query(
        uniprot_accs, select_clause, where_clause, limit, groupby_clause=group_by
    )


@dataclass(frozen=True)
class ComplexPortalEntry:
    """A ComplexPortal entry.

    Parameters:
        query_protein: The UniProt accession used to find entry.
        complex_id: The ComplexPortal identifier (for example "CPX-1234").
        complex_url: The URL to the ComplexPortal entry.
        complex_title: The title of the complex.
        members: UniProt accessions which are members of the complex.
    """

    query_protein: str
    complex_id: str
    complex_url: str
    complex_title: str
    members: set[str]


def _flatten_results_complex(raw_results) -> list[ComplexPortalEntry]:
    results = []
    for raw_result in raw_results:
        query_protein = raw_result["protein"]["value"].split("/")[-1]
        complex_id = raw_result["cp_db"]["value"].split("/")[-1]
        complex_url = f"https://www.ebi.ac.uk/complexportal/complex/{complex_id}"
        complex_title = raw_result.get("cp_comment", {}).get("value", "")
        members = {m.split("/")[-1] for m in raw_result["complex_members"]["value"].split(",")}
        results.append(
            ComplexPortalEntry(
                query_protein=query_protein,
                complex_id=complex_id,
                complex_url=complex_url,
                complex_title=complex_title,
                members=members,
            )
        )
    return results


def search4macromolecular_complexes(
    uniprot_accs: Iterable[str], limit: int = 10_000, timeout: int = 1_800
) -> list[ComplexPortalEntry]:
    """Search for macromolecular complexes by UniProtKB accessions.

    Queries for references to/from https://www.ebi.ac.uk/complexportal/ database in the Uniprot SPARQL endpoint.

    Args:
        uniprot_accs: UniProt accessions.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.

    Returns:
        List of ComplexPortalEntry objects.
    """
    sparql_query = _build_complex_sparql_query(uniprot_accs, limit)
    logger.info("Executing SPARQL query for macromolecular complexes: %s", sparql_query)
    raw_results = _execute_sparql_search(
        sparql_query=sparql_query,
        timeout=timeout,
    )
    limit_check("Search for complexes", limit, len(raw_results))
    return _flatten_results_complex(raw_results)


def search4interaction_partners(
    uniprot_acc: str, excludes: set[str] | None = None, limit: int = 10_000, timeout: int = 1_800
) -> dict[str, set[str]]:
    """Search for interaction partners of a given UniProt accession using ComplexPortal database references.

    Args:
        uniprot_acc: UniProt accession to search interaction partners for.
        excludes: Set of UniProt accessions to exclude from the results.
            For example already known interaction partners.
            If None then no complex members are excluded.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.

    Returns:
        Dictionary with UniProt accessions of interaction partners as keys and sets of ComplexPortal entry IDs
        in which the interaction occurs as values.
    """
    ucomplexes = search4macromolecular_complexes([uniprot_acc], limit=limit, timeout=timeout)
    hits: dict[str, set[str]] = {}
    if excludes is None:
        excludes = set()
    for ucomplex in ucomplexes:
        for member in ucomplex.members:
            if member != uniprot_acc and member not in excludes:
                if member not in hits:
                    hits[member] = set()
                hits[member].add(ucomplex.complex_id)
    return hits
