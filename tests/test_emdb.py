from pathlib import Path

import pytest

from protein_quest.emdb import fetch


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_fetch(tmp_path: Path):
    # use small emdb entry
    emdb_ids = ["EMD-1470"]

    results = await fetch(emdb_ids, tmp_path)
    expected = {"EMD-1470": tmp_path / "emd_1470.map.gz"}
    assert results == expected
    assert all(path.exists() for path in results.values())
