"""Module for searching UniProtKB using SPARQL."""

import logging
from collections.abc import Collection, Generator, Iterable
from dataclasses import dataclass
from functools import cached_property
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
        min_sequence_length: Minimum length of the canonical sequence.
        max_sequence_length: Maximum length of the canonical sequence.
    """

    # TODO make taxon_id an int
    taxon_id: str | None
    reviewed: bool | None = None
    subcellular_location_uniprot: str | None = None
    subcellular_location_go: list[str] | None = None
    molecular_function_go: list[str] | None = None
    min_sequence_length: int | None = None
    max_sequence_length: int | None = None


def _first_chain_from_uniprot_chains(uniprot_chains: str) -> str:
    """Extracts the first chain identifier from a UniProt chains string.

    The UniProt chains string is formatted (with EBNF notation) as follows:

        chain_group=range(,chain_group=range)*

    where:
        chain_group := chain_id(/chain_id)*
        chain_id    := [A-Za-z0-9]+
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


def _chain_length_from_uniprot_chains(uniprot_chains: str) -> int:
    """Calculates the total length of chain from a UniProt chains string.

    See `_first_chain_from_uniprot_chains` for the format of the UniProt chains string.

    Args:
        uniprot_chains: A string representing UniProt chains, For example "B/D=1-81".

    Returns:
        The length of the chain in the UniProt chain string. For example 81 for "B/D=1-81".
    """
    total_length = 0
    chains = uniprot_chains.split(",")
    for chain in chains:
        _, rangestr = chain.split("=")
        start, stop = rangestr.split("-")
        # Residue positions are 1-based so + 1
        total_length += int(stop) - int(start) + 1
    return total_length


class PdbChainLengthError(ValueError):
    """Raised when a UniProt chain description does not yield a chain length."""

    def __init__(self, pdb_id: str, uniprot_chains: str):
        msg = f"Could not determine chain length of '{pdb_id}' from '{uniprot_chains}'"
        super().__init__(msg)


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

    @cached_property
    def chain(self) -> str:
        """The first chain from the UniProt chains aka self.uniprot_chains."""
        return _first_chain_from_uniprot_chains(self.uniprot_chains)

    @cached_property
    def chain_length(self) -> int:
        """The length of the chain from the UniProt chains aka self.uniprot_chains."""
        try:
            return _chain_length_from_uniprot_chains(self.uniprot_chains)
        except ValueError as e:
            raise PdbChainLengthError(self.id, self.uniprot_chains) from e


type PdbResults = dict[str, set[PdbResult]]
"""Dictionary with uniprot accessions as keys and sets of PDB results as values."""


def filter_pdb_results_on_chain_length(
    pdb_results: PdbResults,
    min_residues: int | None,
    max_residues: int | None,
    keep_invalid: bool = False,
) -> PdbResults:
    """Filter PDB results based on chain length.

    Args:
        pdb_results: Dictionary with protein IDs as keys and sets of PDB results as values.
        min_residues: Minimum number of residues required in the chain mapped to the UniProt accession.
            If None, no minimum is applied.
        max_residues: Maximum number of residues allowed in chain mapped to the UniProt accession.
            If None, no maximum is applied.
        keep_invalid: If True, PDB results with invalid chain length (could not be determined) are kept.
            If False, PDB results with invalid chain length are filtered out.
            Warnings are logged when length can not be determined.

    Returns:
        Filtered dictionary with protein IDs as keys and sets of PDB results as values.
    """
    if min_residues is None and max_residues is None:
        # No filtering needed
        return pdb_results
    if min_residues is not None and max_residues is not None and max_residues <= min_residues:
        msg = f"Maximum number of residues ({max_residues}) must be > minimum number of residues ({min_residues})"
        raise ValueError(msg)
    results: PdbResults = {}
    for uniprot_accession, pdb_entries in pdb_results.items():
        filtered_pdb_entries = set()
        for pdb_entry in pdb_entries:
            try:
                if (min_residues is None or pdb_entry.chain_length >= min_residues) and (
                    max_residues is None or pdb_entry.chain_length <= max_residues
                ):
                    filtered_pdb_entries.add(pdb_entry)
            except PdbChainLengthError:
                if keep_invalid:
                    logger.warning(
                        f"Could not determine chain length of '{pdb_entry.id}' from '{pdb_entry.uniprot_chains}' "
                        f"belonging to uniprot accession '{uniprot_accession}', "
                        "for completeness not filtering it out"
                    )
                    filtered_pdb_entries.add(pdb_entry)
                else:
                    logger.warning(
                        f"Filtering out PDB entry '{pdb_entry.id}' belonging to uniprot accession "
                        f"'{uniprot_accession}' due to invalid chain length from '{pdb_entry.uniprot_chains}'"
                    )
        if filtered_pdb_entries:
            # Only include uniprot_accession if there are any pdb entries left after filtering
            results[uniprot_accession] = filtered_pdb_entries
    return results


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

    if query.min_sequence_length is not None or query.max_sequence_length is not None:
        length_filter = _build_sparql_query_sequence_length_filter(
            min_length=query.min_sequence_length,
            max_length=query.max_sequence_length,
        )
        parts.append(length_filter)

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


def _build_sparql_query_sequence_length_filter(min_length: int | None = None, max_length: int | None = None) -> str:
    """Builds a SPARQL filter for sequence length.

    See 107_uniprot_sequences_and_mark_which_is_cannonical_for_human
    on https://sparql.uniprot.org/.well-known/sparql-examples/ for similar query.

    Args:
        min_length: Minimum sequence length. If None, no minimum is applied.
        max_length: Maximum sequence length. If None, no maximum is applied.
    """
    if min_length is None and max_length is None:
        return ""
    # An uniprot entry can have multiple isoforms,
    # we want to check the length of the canonical isoform
    # We do this by selecting the isoform that is not based on another isoform
    # and excluding isoforms from other uniprot entries.
    # For example for http://purl.uniprot.org/uniprot/P42284:
    # - http://purl.uniprot.org/isoforms/P42284-2 is ok
    # - http://purl.uniprot.org/isoforms/P42284-1 is not ok, because it is based on P42284-2
    # - http://purl.uniprot.org/isoforms/Q7KQZ4-1 is not ok, because it is from another uniprot entry
    header = dedent("""\
        ?protein up:sequence ?isoform .
        ?isoform a up:Simple_Sequence .
        BIND (IRI(STRBEFORE(REPLACE(
            STR(?isoform), "http://purl.uniprot.org/isoforms/", "http://purl.uniprot.org/uniprot/"
        ), "-")) AS ?ac_of_isoform)
        FILTER (?protein = ?ac_of_isoform)
        ?isoform rdf:value ?sequence .
        BIND (STRLEN(?sequence) AS ?seq_length)
    """)
    if min_length is not None and max_length is not None:
        if max_length <= min_length:
            msg = f"Maximum sequence length ({max_length}) must be greater than minimum sequence length ({min_length})"
            raise ValueError(msg)
        return dedent(f"""\
            {header}
            FILTER (?seq_length >= {min_length} && ?seq_length <= {max_length})
        """)
    if min_length is not None:
        return dedent(f"""\
            {header}
            FILTER (?seq_length >= {min_length})
        """)
    if max_length is not None:
        return dedent(f"""\
            {header}
            FILTER (?seq_length <= {max_length})
        """)
    return ""


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
    # Also there can be multiple chains for the same uniprot/pdb pair, so we need to
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


def _build_sparql_query_af(
    uniprot_accs: Iterable[str],
    min_sequence_length: int | None = None,
    max_sequence_length: int | None = None,
    limit=10_000,
) -> str:
    select_clause = "?protein ?af_db"
    where_clause = dedent("""
        # --- Protein Selection ---
        ?protein a up:Protein .

        # --- AlphaFoldDB Info ---
        ?protein rdfs:seeAlso ?af_db .
        ?af_db up:database <http://purl.uniprot.org/database/AlphaFoldDB> .
    """)
    if min_sequence_length is not None or max_sequence_length is not None:
        length_filter = _build_sparql_query_sequence_length_filter(
            min_length=min_sequence_length,
            max_length=max_sequence_length,
        )
        where_clause += "\n" + length_filter
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


def _flatten_results_pdb(rawresults: Iterable) -> PdbResults:
    pdb_entries: PdbResults = {}
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
) -> PdbResults:
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
    uniprot_accs: Collection[str],
    min_sequence_length: int | None = None,
    max_sequence_length: int | None = None,
    limit: int = 10_000,
    timeout: int = 1_800,
    batch_size: int = 10_000,
) -> dict[str, set[str]]:
    """
    Search for AlphaFold entries in UniProtKB accessions.

    Args:
        uniprot_accs: UniProt accessions.
        min_sequence_length: Minimum length of the canonical sequence.
        max_sequence_length: Maximum length of the canonical sequence.
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
            sparql_query = _build_sparql_query_af(batch, min_sequence_length, max_sequence_length, limit)
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
    (GROUP_CONCAT(
        DISTINCT STRAFTER(STR(?member), "http://purl.uniprot.org/uniprot/"); separator=","
    ) AS ?complex_members)
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
        (GROUP_CONCAT(
            DISTINCT STRAFTER(STR(?member), "http://purl.uniprot.org/uniprot/"); separator=","
        ) AS ?complex_members)
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
        members = set(raw_result["complex_members"]["value"].split(","))
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
    uniprot_accession: str, excludes: set[str] | None = None, limit: int = 10_000, timeout: int = 1_800
) -> dict[str, set[str]]:
    """Search for interaction partners of a given UniProt accession using ComplexPortal database references.

    Args:
        uniprot_accession: UniProt accession to search interaction partners for.
        excludes: Set of UniProt accessions to exclude from the results.
            For example already known interaction partners.
            If None then no complex members are excluded.
        limit: Maximum number of results to return.
        timeout: Timeout for the SPARQL query in seconds.

    Returns:
        Dictionary with UniProt accessions of interaction partners as keys and sets of ComplexPortal entry IDs
        in which the interaction occurs as values.
    """
    ucomplexes = search4macromolecular_complexes([uniprot_accession], limit=limit, timeout=timeout)
    hits: dict[str, set[str]] = {}
    if excludes is None:
        excludes = set()
    for ucomplex in ucomplexes:
        for member in ucomplex.members:
            if member != uniprot_accession and member not in excludes:
                if member not in hits:
                    hits[member] = set()
                hits[member].add(ucomplex.complex_id)
    return hits


@dataclass(frozen=True)
class UniprotDetails:
    """Details of an UniProt entry.

    Parameters:
        uniprot_accession: UniProt accession.
        uniprot_id: UniProt ID (mnemonic).
        sequence_length: Length of the canonical sequence.
        reviewed: Whether the entry is reviewed (Swiss-Prot) or unreviewed (TrEMBL).
        protein_name: Recommended protein name.
        taxon_id: NCBI Taxonomy ID of the organism.
        taxon_name: Scientific name of the organism.
    """

    uniprot_accession: str
    uniprot_id: str
    sequence_length: int
    reviewed: bool
    protein_name: str
    taxon_id: int
    taxon_name: str


def map_uniprot_accessions2uniprot_details(
    uniprot_accessions: Collection[str], timeout: int = 1_800, batch_size: int = 1000
) -> Generator[UniprotDetails]:
    """Map UniProt accessions to UniProt details by querying the UniProt SPARQL endpoint.

    Example:

    SPARQL query to get details for 7 UniProt entries, run on [https://sparql.uniprot.org/sparql](https://sparql.uniprot.org/sparql).

    ```sparql
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
    PREFIX up:   <http://purl.uniprot.org/core/>
    PREFIX rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT
    (?ac AS ?uniprot_accession)
    ?uniprot_id
    (STRAFTER(STR(?organism), "taxonomy/") AS ?taxon_id)
    ?taxon_name
    ?reviewed
    ?protein_name
    (STRLEN(?sequence) AS ?seq_length)
    WHERE {
    # Input UniProt accessions
    VALUES (?ac) { ("P05067") ("A6NGD5") ("O14627") ("P00697") ("P42284") ("A0A0B5AC95") ("A0A0S2Z4R0")}
    BIND (IRI(CONCAT("http://purl.uniprot.org/uniprot/", ?ac)) AS ?protein)
    ?protein a up:Protein .
    ?protein up:mnemonic ?uniprot_id .
    ?protein up:organism ?organism .
    ?organism up:scientificName ?taxon_name .
    ?protein up:reviewed ?reviewed .
    ?protein up:recommendedName/up:fullName ?protein_name .
    ?protein up:sequence ?isoform .
    ?isoform a up:Simple_Sequence .
    ?isoform rdf:value ?sequence .
    BIND (IRI(STRBEFORE(REPLACE(
        STR(?isoform), "http://purl.uniprot.org/isoforms/", "http://purl.uniprot.org/uniprot/"
    ), "-")) AS ?ac_of_isoform)
    FILTER(?ac_of_isoform = ?protein)
    }
    ```

    Args:
        uniprot_accessions: Iterable of UniProt accessions.
        timeout: Timeout for the SPARQL query in seconds.
        batch_size: Size of batches to process the UniProt accessions.

    Yields:
        UniprotDetails objects in random order.
    """
    select_clause = dedent("""\
        (?ac AS ?uniprot_accession)
        ?uniprot_id
        (STRAFTER(STR(?organism), "taxonomy/") AS ?taxon_id)
        ?taxon_name
        ?reviewed
        ?protein_name
        (STRLEN(?sequence) AS ?seq_length)
    """)
    where_clause = dedent("""
        ?protein up:mnemonic ?uniprot_id .
        ?protein up:organism ?organism .
        ?organism up:scientificName ?taxon_name .
        ?protein up:reviewed ?reviewed .
        OPTIONAL {
        ?protein up:recommendedName/up:fullName ?protein_name .
        }
        ?protein up:sequence ?isoform .
        ?isoform a up:Simple_Sequence .
        ?isoform rdf:value ?sequence .
        BIND (IRI(STRBEFORE(REPLACE(
            STR(?isoform), "http://purl.uniprot.org/isoforms/", "http://purl.uniprot.org/uniprot/"
        ), "-")) AS ?ac_of_isoform)
        FILTER(?ac_of_isoform = ?protein)
    """)
    total = len(uniprot_accessions)
    with tqdm(
        total=total,
        desc="Retrieving UniProt details",
        disable=total < batch_size,
        unit="acc",
    ) as pbar:
        for batch in batched(uniprot_accessions, batch_size, strict=False):
            sparql_query = _build_sparql_generic_by_uniprot_accessions_query(
                batch, select_clause, where_clause, limit=batch_size
            )
            logger.info("Executing SPARQL query for UniProt details: %s", sparql_query)
            raw_results = _execute_sparql_search(
                sparql_query=sparql_query,
                timeout=timeout,
            )
            for raw_result in raw_results:
                protein_name = raw_result.get("protein_name", {}).get("value", "")
                result = UniprotDetails(
                    uniprot_accession=raw_result["uniprot_accession"]["value"],
                    uniprot_id=raw_result["uniprot_id"]["value"],
                    sequence_length=int(raw_result["seq_length"]["value"]),
                    reviewed=raw_result["reviewed"]["value"] == "true",
                    protein_name=protein_name,
                    taxon_id=int(raw_result["taxon_id"]["value"]),
                    taxon_name=raw_result["taxon_name"]["value"],
                )
                yield result
            pbar.update(len(batch))
