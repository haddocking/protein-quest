"""Wrappers around [PDBe REST API](https://www.ebi.ac.uk/pdbe/api/v2/doc/)."""

from asyncio import sleep
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import batched

from tqdm.rich import tqdm

from protein_quest.converter import converter
from protein_quest.utils import RetryClient, friendly_session


@dataclass(frozen=True, slots=True)
class Scores:
    """PDBe validation summary quality scores

    These scores are harmonic means of percentile-based quality metrics for macromolecular structures:

    * geometry_quality: from geometric percentiles (Ramachandran, clashscore, sidechains).
    * data_quality: from diffraction-data percentiles (R-free, RSRZ).
    * overall_quality: combines geometry and data percentiles.
      If any contributing percentile is 0, the harmonic mean is 0.
      If all percentiles are unavailable, the score is null.
    * experiment_data_available: indicates whether experimental diffraction data were used.

    Attributes:
        geometry_quality: Harmonic mean of all absolute validation percentile metrics related to model geometry.
        data_quality: Harmonic mean of absolute percentiles related to model geometry.
        overall_quality: Harmonic mean of absolute percentiles related to experimental data and its fit to the model.
        experiment_data_available: Sometimes data quality is absent due to unavailability of experimental data itself.
            This flag tells whether the data are available.
    """

    geometry_quality: float | None
    data_quality: float | None
    overall_quality: float | None
    experiment_data_available: bool | str


async def _fetch_summary_quality_scores(
    batch: Iterable[str],
    session: RetryClient,
    url: str = "https://www.ebi.ac.uk/pdbe/api/v2/validation/summary_quality_scores/entry",
) -> dict[str, Scores]:
    request = converter.dumps(",".join(batch))
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    async with session.post(url, data=request, headers=headers) as response:
        response.raise_for_status()
        raw_response = await response.read()
        return converter.loads(raw_response, dict[str, Scores])


async def fetch_summary_quality_scores_in_batches(
    pdb_ids: set[str],
    timeout: int = 1800,
    batch_size: int = 20,
    sleep_between_batches: float = 0.1,
) -> dict[str, Scores]:
    """Fetches validation summary quality scores for a set of PDB IDs in batches.

    Args:
        pdb_ids: A set of PDB IDs to fetch quality scores for.
        timeout: Total timeout for fetching all batches, in seconds.
        batch_size: Number of PDB IDs to include in each batch request.
        sleep_between_batches: Time to sleep between batch requests, in seconds.

    Returns:
        A dict mapping PDB IDs to their quality scores.

    Raises:
        aiohttp.ClientResponseError: If any batch request fails with a non-2xx status code.
    """
    total = len(pdb_ids)
    scores = {}
    async with friendly_session(total_timeout=timeout) as session:
        with tqdm(
            total=total,
            desc="Searching for quality scores for PDBe entries",
            disable=total < batch_size,
            unit="acc",
        ) as pbar:
            for batch in batched(pdb_ids, batch_size, strict=False):
                batch_summries = await _fetch_summary_quality_scores(batch, session)
                scores.update(batch_summries)
                pbar.update(len(batch))
                await sleep(sleep_between_batches)
    return scores
