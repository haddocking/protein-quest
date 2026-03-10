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
from typing import TypedDict, get_args

from aiohttp_retry import RetryClient
from attrs import define, field, validators
from tqdm.rich import tqdm

from protein_quest.converter import PositiveInt, converter
from protein_quest.pdbe_3dbeacons.model import (
    AccessionListRequest,
    AppUniprotSchemaSummaryItems,
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

DEFAULT_PROVIDERS: set[Provider] = {"pdbe", "alphafold"}


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
        default=DEFAULT_PROVIDERS,
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


def _find_chain_for_uniprot(uniprot_accession: str, summary: AppUniprotSchemaSummaryItems) -> str:
    if entities_of_uniprot := [e for e in summary.entities if e.identifier == uniprot_accession]:
        # Surest way to get the correct chain,
        # only supported by PDBe, SWISS-model and Alphafold DB
        return entities_of_uniprot[0].chain_ids[0]
    if summary.entities:
        # Most providers have the POLYMER entity as first entity with single chain.
        return summary.entities[0].chain_ids[0]
    # A PED structure has zero entities, but its cif file has chain A
    return "A"


class FlattenedUniprotSummary(TypedDict):
    """Typed representation of a flattened structure summary.

    Attributes:
        uniprot_accession: Uniprot accession.
        provider: [Provider][protein_quest.pdbe_3dbeacons.model.Provider] of the structure.
        model_identifier: Model identifier of the structure.
        model_url: URL to download the structure.
        model_format: [Format][protein_quest.pdbe_3dbeacons.model.AppUniprotSchemaModelFormat]
          of the structure file
        chain: Chain identifier of the structure (first chain of first entity or "A" if no entities or chains).
        residue_count: Number of residues in the structure
    """

    uniprot_accession: str
    provider: str
    model_identifier: str
    model_url: str
    model_format: str
    chain: str
    residue_count: int


def _sum_residue_counts_for_same_models(summaries: list[FlattenedUniprotSummary]) -> list[FlattenedUniprotSummary]:
    """Sum the residue counts for summaries with the same provider and model identifier.

    Needed for uniprot acc P38634 as chain y in 8k0g pdb entry which has residues: 2 — 5, 6 — 48
    API returns two summaries with `2 — 5` and `6 — 48` as uniprot start - end.
    Here we sum them 4 + 43 to get the total residue count of 47 for the model.
    """
    model_chain_to_row: dict[tuple[str, str], FlattenedUniprotSummary] = {}
    for summary in summaries:
        provider_id = (summary["provider"], summary["model_identifier"])
        if provider_id in model_chain_to_row:
            model_chain_to_row[provider_id]["residue_count"] += summary["residue_count"]
        else:
            model_chain_to_row[provider_id] = summary
    return list(model_chain_to_row.values())


def flatten_structure_summaries(summaries: list[UniprotSummary]) -> list[FlattenedUniprotSummary]:
    """Flatten the summaries to a list of dicts with uniprot accession and structure information.

    Args:
        summaries: List of summaries to flatten.

    Returns:
        List of dicts.
    """
    provider_response2request = {v: k for k, v in provider_request2response.items()}
    rows: list[FlattenedUniprotSummary] = []
    for summary in summaries:
        if summary.uniprot_entry is None or summary.structures is None:
            continue
        uniprot_accession = summary.uniprot_entry.ac
        for structure_summary in summary.structures:
            s = structure_summary.summary
            provider = provider_response2request[s.provider]
            chain = _find_chain_for_uniprot(uniprot_accession, s)
            residue_count = s.uniprot_end - s.uniprot_start + 1
            row: FlattenedUniprotSummary = {
                "uniprot_accession": uniprot_accession,
                "provider": provider,
                "model_identifier": s.model_identifier,
                "model_url": s.model_url,
                "model_format": s.model_format,
                "chain": chain,
                "residue_count": residue_count,
            }
            rows.append(row)
    return _sum_residue_counts_for_same_models(rows)
