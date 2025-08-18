"""Module for Gene Ontology (GO) functions."""

import logging
from collections.abc import Generator
from dataclasses import dataclass
from typing import Literal

from cattrs.gen import make_dict_structure_fn, override
from cattrs.preconf.orjson import make_converter

from protein_quest.utils import friendly_session

logger = logging.getLogger(__name__)

Aspect = Literal["cellular_component", "biological_process", "molecular_function"]
"""The aspect of the GO term."""


@dataclass(frozen=True, slots=True)
class GoTerm:
    """A Gene Ontology (GO) term.

    Parameters:
        id: The unique identifier for the GO term, e.g., 'GO:0043293'.
        is_obsolete: Whether the GO term is obsolete.
        name: The name of the GO term.
        definition: The definition of the GO term.
        aspect: The aspect of the GO term.
    """

    id: str
    is_obsolete: bool
    name: str
    definition: str
    aspect: Aspect


@dataclass(frozen=True, slots=True)
class PageInfo:
    current: int
    total: int


@dataclass(frozen=True, slots=True)
class SearchResponse:
    results: list[GoTerm]
    number_of_hits: int
    page_info: PageInfo


converter = make_converter()


def flatten_definition(definition, _context) -> str:
    return definition["text"]


# Use hook to convert incoming camelCase to snake_case
# and to flatten definition {text} to text
# see https://catt.rs/en/stable/customizing.html#rename
converter.register_structure_hook(
    GoTerm,
    make_dict_structure_fn(
        GoTerm,
        converter,
        is_obsolete=override(rename="isObsolete"),
        definition=override(struct_hook=flatten_definition),
    ),
)
converter.register_structure_hook(
    SearchResponse,
    make_dict_structure_fn(
        SearchResponse, converter, number_of_hits=override(rename="numberOfHits"), page_info=override(rename="pageInfo")
    ),
)


async def search_go_term(
    term: str, aspect: Aspect | None = None, include_obsolete: bool = False, limit: int = 100
) -> list[GoTerm]:
    """Search for a Gene Ontology (GO) term by its name or ID.

    Calls the EBI QuickGO API at https://www.ebi.ac.uk/QuickGO/api/index.html .

    Examples:
        To search for `apoptosome` terms do.

        >>> from protein_quest.go import search_go_term
        >>> r = await search_go_term('apoptosome')
        >>> len(r)
        5
        >>> r[0]
        GoTerm(id='GO:0043293', is_obsolete=False, name='apoptosome', definition='A multisubunit protein ...')

    Args:
        term: The GO term to search for. For example `nucleus` or `GO:0006816`.
        aspect: The aspect to filter by. If not given, all aspects are included.
        include_obsolete: Whether to include obsolete terms. By default, obsolete terms are excluded.
        limit: The maximum number of results to return.

    Returns:
        List of GO terms
    """
    url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search"
    page_limit = 100
    params = {"query": term, "limit": str(page_limit), "page": "1"}
    logger.debug("Fetching GO terms from %s with params %s", url, params)
    async with friendly_session() as session:
        # Fetch first page to learn how many pages there are
        async with session.get(url, params=params) as response:
            response.raise_for_status()
            raw_data = await response.read()
            data = converter.loads(raw_data, SearchResponse)

        terms = list(_filter_go_terms(data.results, aspect, include_obsolete))
        if len(terms) >= limit:
            # Do not fetch additional pages if we have enough results
            return terms[:limit]
        total_pages = data.page_info.total
        logger.debug("GO search returned %s pages (current=%s)", total_pages, data.page_info.current)

        # Retrieve remaining pages (if any) and extend results
        if total_pages > 1:
            for page in range(2, total_pages + 1):
                params["page"] = str(page)
                logger.debug("Fetching additional GO terms page %s/%s with params %s", page, total_pages, params)
                async with session.get(url, params=params) as response:
                    response.raise_for_status()
                    raw_data = await response.read()
                    data = converter.loads(raw_data, SearchResponse)
                terms.extend(_filter_go_terms(data.results, aspect, include_obsolete))
                if len(terms) >= limit:
                    # Do not fetch additional pages if we have enough results
                    break

        return terms[:limit]


def _filter_go_terms(terms: list[GoTerm], aspect: Aspect | None, include_obsolete: bool) -> Generator[GoTerm]:
    for oboterm in terms:
        if not include_obsolete and oboterm.is_obsolete:
            continue
        if aspect and oboterm.aspect != aspect:
            continue
        yield oboterm
