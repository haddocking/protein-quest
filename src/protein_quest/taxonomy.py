"""Module for searching taxon information from UniProt."""

import csv
import gzip
import logging
from dataclasses import dataclass
from typing import Literal, get_args

from aiohttp.client import ClientResponse
from aiohttp_retry import RetryClient
from cattrs.gen import make_dict_structure_fn, override
from yarl import URL

from protein_quest.converter import converter
from protein_quest.go import TextIOWrapper
from protein_quest.utils import friendly_session

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class Taxon:
    """Dataclass representing a taxon.

    Arguments:
        taxon_id: The unique identifier for the taxon.
        scientific_name: The scientific name of the taxon.
        rank: The taxonomic rank of the taxon (e.g., species, genus).
        common_name: The common name of the taxon (if available).
        other_names: A set of other names for the taxon (if available).
    """

    taxon_id: str
    scientific_name: str
    rank: str
    common_name: str | None = None
    other_names: set[str] | None = None


@dataclass(frozen=True, slots=True)
class SearchTaxonResponse:
    results: list[Taxon]


converter.register_structure_hook(
    Taxon,
    make_dict_structure_fn(
        Taxon,
        converter,
        taxon_id=override(rename="taxonId"),
        scientific_name=override(rename="scientificName"),
        common_name=override(rename="commonName"),
        other_names=override(rename="otherNames"),
    ),
)

SearchField = Literal["tax_id", "scientific", "common", "parent"]
"""Type of search field"""
search_fields: set[SearchField | None] = set(get_args(SearchField)) | {None}
"""Set of valid search fields"""


def _get_next_page(response: ClientResponse) -> URL | str | None:
    next_page = response.links.getone("next", None)
    if next_page is None:
        return None
    return next_page.getone("url", None)


async def _fetch_page(url: URL | str, session: RetryClient) -> tuple[list[Taxon], URL | str | None]:
    async with session.get(url) as response:
        response.raise_for_status()
        gzipped_raw_data = await response.read()
    next_page = _get_next_page(response)
    raw_data = gzip.decompress(gzipped_raw_data)
    taxons = converter.loads(raw_data, SearchTaxonResponse).results
    return taxons, next_page


async def search_taxon(query: str, field: SearchField | None = None, limit: int = 100) -> list[Taxon]:
    """Search for taxon information in UniProt.

    Uses <https://www.uniprot.org/taxonomy?query=*>.

    Args:
        query: Search query for the taxon.
        field: Field to search in.
            If None, searches in all fields.
            If "tax_id" then searches by taxon ID.
            If "parent" then given a parent taxon ID returns all its children.
            For example, if the parent taxon ID is 9606 (Human), it will return Neanderthal and Denisovan.
        limit: Maximum number of results to return.

    Returns:
        List of Taxon objects matching the search query.

    Raises:
        ValueError: If the search field is invalid.
    """
    # https://rest.uniprot.org/taxonomy/search?compressed=true&format=json&query=%28Mouse%29&size=500
    page_limit = 100
    if field not in search_fields:
        msg = f"Invalid search field: {field}. Must be one of {search_fields}."
        raise ValueError(msg)
    if field is not None:
        # ((common:"house+mouse"))
        query = f'(({field}:"{query}"))'
    params = {"query": query, "limit": str(page_limit), "compressed": "true", "format": "json"}
    url = URL("https://rest.uniprot.org/taxonomy/search").with_query(params)
    logger.debug("Fetching uniprot taxonomy from %s with params %s", url, params)
    async with friendly_session() as session:
        # Fetch first page
        taxons, next_page = await _fetch_page(url, session)
        if len(taxons) >= limit:
            return taxons[:limit]
        if next_page is None:
            return taxons

        # Fetch next pages
        while next_page:
            logger.debug("Fetching next page of uniprot taxonomy from %s", next_page)
            next_data, next_page = await _fetch_page(next_page, session)
            taxons.extend(next_data)
            if len(taxons) >= limit:
                return taxons[:limit]
            if next_page is None:
                return taxons
    return taxons


def _write_taxonomy_csv(taxons: list[Taxon], output_csv: TextIOWrapper) -> None:
    """Write taxon information to a CSV file.

    Args:
        taxons: List of Taxon objects to write to the CSV file.
        output_csv: File object for the output CSV file.
    """
    writer = csv.writer(output_csv)
    writer.writerow(["taxon_id", "scientific_name", "common_name", "rank", "other_names"])
    for taxon in taxons:
        writer.writerow(
            [
                taxon.taxon_id,
                taxon.scientific_name,
                taxon.common_name,
                taxon.rank,
                ";".join(taxon.other_names) if taxon.other_names else "",
            ]
        )
