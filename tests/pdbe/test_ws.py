import pytest
from aiohttp import ClientResponseError

from protein_quest.pdbe.ws import Scores, fetch_summary_quality_scores_in_batches


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_fetch_summary_quality_scores_in_batches_happypath():
    # Ids from fixtures in conftest.py
    pdb_ids = {
        "3JRS",
        "2Y29",
        "1AMB",
        "6O5I",
        "1A02",
        "1UN5",
    }
    scores = await fetch_summary_quality_scores_in_batches(pdb_ids, batch_size=4)

    expected = {
        "1a02": Scores(
            geometry_quality=6.15,
            data_quality=7.84,
            overall_quality=6.73,
            experiment_data_available=True,
        ),
        "1amb": Scores(
            geometry_quality=None,
            data_quality=None,
            overall_quality=None,
            experiment_data_available="unknown",
        ),
        "1un5": Scores(
            geometry_quality=45.23,
            data_quality=None,
            overall_quality=45.23,
            experiment_data_available=False,
        ),
        "2y29": Scores(
            geometry_quality=55.9,
            data_quality=81.09,
            overall_quality=63.83,
            experiment_data_available=True,
        ),
        "3jrs": Scores(
            geometry_quality=31.58,
            data_quality=22.4,
            overall_quality=27.13,
            experiment_data_available=True,
        ),
        "6o5i": Scores(
            geometry_quality=64.38,
            data_quality=61.45,
            overall_quality=63.17,
            experiment_data_available=True,
        ),
    }
    assert scores == expected


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_fetch_summary_quality_scores_in_batches_badid():
    pdb_ids = {"BADID"}
    with pytest.raises(ClientResponseError) as excinfo:
        await fetch_summary_quality_scores_in_batches(pdb_ids, batch_size=4)

    assert excinfo.value.status == 404
