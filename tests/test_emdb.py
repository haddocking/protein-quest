from pathlib import Path

import pytest

from protein_quest.emdb import fetch


@pytest.mark.vcr
def test_fetch(tmp_path: Path):
    # use small emdb entry
    emdb_ids = ["EMD-1470"]

    results = fetch(emdb_ids, tmp_path)
    expected = {"EMD-1470": tmp_path / "emd_1470.map.gz"}
    assert results == expected
