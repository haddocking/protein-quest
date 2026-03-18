from io import StringIO
from pathlib import Path

import pytest

from protein_quest.pdbe.fetch import fetch, read_pdb_ids_from_csv, sync_fetch


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


class TestReadPdbIdsFromCsv:
    def test_reads_pdb_id_column(self):
        csv_data = StringIO("pdb_id\n2Y29\n8WAS\n2Y29\n")

        ids = read_pdb_ids_from_csv(csv_data)

        assert ids == {"2Y29", "8WAS"}

    def test_reads_model_identifier_for_pdbe_provider(self):
        csv_data = StringIO("model_provider,model_identifier\npdbe,8WAS\npdbe,2Y29\n")

        ids = read_pdb_ids_from_csv(csv_data)

        assert ids == {"8WAS", "2Y29"}
