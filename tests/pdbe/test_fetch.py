from pathlib import Path

import pytest

from protein_quest.pdbe.fetch import fetch


@pytest.mark.vcr
def test_fetch(tmp_path: Path):
    theid = "2Y29"
    ids = [theid]

    results = fetch(ids, tmp_path)

    expected = {theid: tmp_path / f"{theid.lower()}.cif.gz"}
    assert results == expected
