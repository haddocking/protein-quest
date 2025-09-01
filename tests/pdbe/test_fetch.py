from pathlib import Path

import pytest

from protein_quest.pdbe.fetch import fetch, sync_fetch


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_fetch(tmp_path: Path):
    theid = "2Y29"
    ids = [theid]

    results = await fetch(ids, tmp_path)

    expected = {theid: tmp_path / f"{theid.lower()}.cif.gz"}
    assert results == expected


@pytest.mark.default_cassette("test_fetch.yaml")  # pyrefly: ignore[not-callable]
@pytest.mark.vcr
def test_sync_fetch(tmp_path: Path):
    theid = "2Y29"
    ids = [theid]

    results = sync_fetch(ids, tmp_path)

    expected = {theid: tmp_path / f"{theid.lower()}.cif.gz"}
    assert results == expected
