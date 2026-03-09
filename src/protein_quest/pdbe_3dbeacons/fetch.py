"""Module for searching structures from [3D beacons network](https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/).

3D beacons HUB API information:

* docs at <https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/>
* openapi spec at <https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/openapi.json>
* code at <https://github.com/3D-Beacons/3d-beacons-hub-api/>
"""

import logging
from asyncio import sleep
from collections.abc import Generator, Iterable
from itertools import batched
from typing import cast, get_args

from aiohttp_retry import RetryClient
from attrs import define, field, validators
from tqdm.rich import tqdm

from protein_quest.converter import PositiveInt, converter
from protein_quest.pdbe_3dbeacons.model import (
    AccessionListRequest,
    Overview,
    Provider,
    UniprotSummary,
)
from protein_quest.utils import friendly_session

logger = logging.getLogger(__name__)

search_structure_provider_choices = set(get_args(Provider))
"""All providers"""

provider_request2response = {
    "pdbe": "PDBe",
    "ped": "PED",
    "swissmodel": "SWISS-MODEL",
    "alphafold": "AlphaFold DB",
    "sasbdb": "SASBDB",
    "alphafill": "AlphaFill",
    "hegelab": "HEGELAB",
    "modelarchive": "ModelArchive",
    "isoformio": "isoform.io",
    "levylab": "levylab",
}
# When Provider in model module changes, update this mapping accordingly.
"""The response uses slightly different [provider names][protein_quest.pdbe_3dbeacons.model.Provider] than the request.
Use this mapping to convert between them."""

# From https://github.com/3D-Beacons/3d-beacons-hub-api/blob/ffc0648a80701e756cbe88742a96974f1022ac89/app/config/__init__.py#L11
MAX_POST_LIMIT = 10
"""Maximum number of uniprot accessions to query in each batch.
The 3D beacons HUB API has a maximum limit of 10 accessions per POST request."""


@define
class PruneOptions:
    """Options for pruning the summaries.

    Attributes:
        providers: Set of providers to keep in the summaries.
            Default is `pdbe` and `alphafold`.
            Use
            [search_structure_provider_choices][protein_quest.pdbe_3dbeacons.fetch.search_structure_provider_choices]
            for all providers.
        limit: Maximum number of structures per uniprot summary per provider to return.
        min_residues: Minimum number of residues a structure must have to be included.
        max_residues: Maximum number of residues a structure can have to be included.
    """

    providers: set[Provider] = field(
        validator=[
            validators.min_len(1),  # non-empty set required
        ],
        default=cast("set[Provider]", {"pdbe", "alphafold"}),
    )
    limit: PositiveInt = 10_000
    min_residues: PositiveInt | None = None
    max_residues: PositiveInt | None = None


def _prune_on_providers(structures: list[Overview], response_providers: set[str]) -> list[Overview]:
    """Keep structures from allowed response provider names."""
    return [structure for structure in structures if structure.summary.provider in response_providers]


def _prune_on_residues(
    structures: list[Overview],
    min_residues: int | None,
    max_residues: int | None,
) -> list[Overview]:
    """Keep structures with residue counts inside configured bounds."""
    if min_residues is None and max_residues is None:
        return structures

    filtered_on_residues: list[Overview] = []
    for structure in structures:
        residue_count = structure.summary.uniprot_end - structure.summary.uniprot_start + 1
        if (min_residues is None or residue_count >= min_residues) and (
            max_residues is None or residue_count <= max_residues
        ):
            filtered_on_residues.append(structure)
        else:
            logger.debug(
                "Pruning structure %s because residue count %d outside bounds [%s, %s]",
                structure.summary.model_identifier,
                residue_count,
                min_residues if min_residues is not None else "-inf",
                max_residues if max_residues is not None else "inf",
            )
    return filtered_on_residues


def prune_on_limit(
    structures: list[Overview],
    response_providers: set[str],
    limit: int,
    summary: UniprotSummary,
) -> list[Overview]:
    """Keep at most ``limit`` structures per provider."""
    filtered_on_limit: list[Overview] = []
    provider_counts: dict[str, int] = dict.fromkeys(response_providers, 0)
    for structure in structures:
        provider = structure.summary.provider
        if provider_counts[provider] < limit:
            filtered_on_limit.append(structure)
            provider_counts[provider] += 1
        else:
            logger.debug(
                "Pruning structure %s from provider %s for uniprot %s because limit of %d reached",
                structure.summary.model_identifier,
                provider,
                summary.uniprot_entry.ac if summary.uniprot_entry else None,
                limit,
            )
    return filtered_on_limit


def _prune_summary(
    summary: UniprotSummary,
    options: PruneOptions,
    response_providers: set[str],
) -> UniprotSummary | None:
    """Prune a single summary by provider, residue bounds, and per-provider limit."""
    if summary.structures is None:
        return None

    filtered_on_provider = _prune_on_providers(summary.structures, response_providers)
    filtered_on_residues = _prune_on_residues(
        filtered_on_provider,
        min_residues=options.min_residues,
        max_residues=options.max_residues,
    )
    filtered_on_limit = prune_on_limit(
        filtered_on_residues,
        response_providers=response_providers,
        limit=options.limit,
        summary=summary,
    )

    return UniprotSummary(
        uniprot_entry=summary.uniprot_entry,
        structures=filtered_on_limit,
    )


def _prune_summaries(summaries: list[UniprotSummary], options: PruneOptions) -> Generator[UniprotSummary]:
    """Remove structures from unwanted providers from the summaries.

    Webservice can only include or exclude one provider, so summaries includes one or all or all but one provider.

    Args:
        summaries: List of summaries to prune.
        options: Options that control provider filtering and per-provider limits.
    """
    response_providers = {provider_request2response[provider] for provider in options.providers}
    for summary in summaries:
        pruned_summary = _prune_summary(summary, options=options, response_providers=response_providers)
        if pruned_summary is not None:
            yield pruned_summary


async def fetch_summary_batch(
    batch: Iterable[str],
    session: RetryClient,
    prune_options: PruneOptions,
    url="https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/uniprot/summary",
) -> list[UniprotSummary]:
    providers = prune_options.providers
    request = AccessionListRequest(accessions=list(batch))
    provider = None
    if len(providers) == 1:
        provider = next(iter(providers))
    exclude_provider = None
    if len(providers) == len(search_structure_provider_choices) - 1:
        exclude_provider = next(iter(search_structure_provider_choices - providers))
    pdbe_provider: Provider = "pdbe"
    if pdbe_provider not in providers:
        # If pdbe can be excluded then do that, because pdbe has the most entries
        exclude_provider = pdbe_provider
    if provider is not None:
        request.provider = provider
    if exclude_provider is not None:
        request.exclude_provider = exclude_provider

    request = converter.dumps(request)
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    async with session.post(url, data=request, headers=headers) as response:
        response.raise_for_status()
        raw_response = await response.read()
        summaries = converter.loads(raw_response, list[UniprotSummary])
        return list(_prune_summaries(summaries, prune_options))


async def uniprots2structures(
    uniprot_accessions: set[str],
    prune_options: PruneOptions,
    timeout: int = 1800,
    batch_size: int = MAX_POST_LIMIT,
    sleep_between_batches: float = 0.1,
) -> list[UniprotSummary]:
    """For list of uniprot accessions, find structures using 3D beacons HUB API.

    Args:
        uniprot_accessions: Set of uniprot accessions to find structures for.
        prune_options: Options controlling provider filtering and per-provider limits.
        timeout: Maximum seconds to wait for each batch query to complete.
        batch_size: Number of uniprot accessions to query in each batch.
        sleep_between_batches: Seconds to wait between each batch query.

    Raises:
        ValueError: If batch_size is not between 1 and MAX_POST_LIMIT.
            Or if providers is not a subset of search_structure_provider_choices.
        aiohttp.ClientResponseError: If the 3D beacons HUB API returns an error response.

    Returns:
        List of summaries of structures for each uniprot accession.
    """
    if batch_size <= 0 or batch_size > MAX_POST_LIMIT:
        msg = f"Batch size {batch_size} must be between 1 and {MAX_POST_LIMIT}."
        raise ValueError(msg)

    total = len(uniprot_accessions)
    structures = []
    async with friendly_session(total_timeout=timeout) as session:
        with tqdm(
            total=total,
            desc="Searching for structures of uniprots",
            disable=total < batch_size,
            unit="acc",
        ) as pbar:
            for batch in batched(uniprot_accessions, batch_size, strict=False):
                batch_structures = await fetch_summary_batch(batch, session, prune_options=prune_options)
                structures.extend(batch_structures)
                pbar.update(len(batch))
                await sleep(sleep_between_batches)

    return structures


def flatten_structure_summaries(summaries: list[UniprotSummary]) -> Generator[dict[str, str | int]]:
    """Flatten the summaries to a list of dicts with uniprot accession and structure information.

    Args:
        summaries: List of summaries to flatten.

    Yields:
        Dict with uniprot accession, provider, model identifier, model url, model format,
        chain (first chain of first entity or 'A'), and number of residues.
    """
    provider_response2request = {v: k for k, v in provider_request2response.items()}
    for summary in summaries:
        if summary.uniprot_entry is None or summary.structures is None:
            continue
        uniprot_accession = summary.uniprot_entry.ac
        for structure_summary in summary.structures:
            s = structure_summary.summary
            provider = provider_response2request[s.provider]
            # Most providers have the POLYMER entity as first entity with single chain.
            # A PED structure has zero entities, but its cif file has chain A
            chain = s.entities[0].chain_ids[0] if len(s.entities) >= 1 and len(s.entities[0].chain_ids) else "A"
            residue_count = s.uniprot_end - s.uniprot_start + 1
            row = {
                "uniprot_accession": uniprot_accession,
                "provider": provider,
                "model_identifier": s.model_identifier,
                "model_url": s.model_url,
                "model_format": s.model_format,
                "chain": chain,
                "residue_count": residue_count,
            }
            yield row
