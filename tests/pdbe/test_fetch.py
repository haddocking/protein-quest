import tarfile
from pathlib import Path

import pytest

from protein_quest.pdbe.fetch import fetch, fetch_to_tar, sync_fetch


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


@pytest.mark.default_cassette("test_fetch.yaml")  # pyrefly: ignore[not-callable]
@pytest.mark.asyncio
@pytest.mark.vcr
async def test_fetch_to_tar(tmp_path: Path):
    theid = "2Y29"
    tar_path = tmp_path / "downloads.tar"

    results, failures = await fetch_to_tar([theid], tar_path)

    assert failures == {}
    assert results == {theid: f"{theid.lower()}.cif.gz"}
    assert tar_path.exists()
    assert not (tmp_path / f"{theid.lower()}.cif.gz").exists()
    with tarfile.open(tar_path, "r") as tar:
        assert tar.getmember(f"{theid.lower()}.cif.gz")
