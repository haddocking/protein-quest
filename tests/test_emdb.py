from pathlib import Path

import pytest

from protein_quest.emdb import fetch, read_emdb_ids_from_csv


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_fetch(tmp_path: Path):
    # use small emdb entry
    emdb_ids = ["EMD-1470"]

    results = await fetch(emdb_ids, tmp_path)
    expected = {"EMD-1470": tmp_path / "emd_1470.map.gz"}
    assert results == expected
    assert all(path.exists() for path in results.values())


class TestReadEMDBIdsFromCSV:
    def test_with_header(self, tmp_path: Path):
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("emdb_id,uniprot_acc\nEMD-1470,ABC123\nEMD-1480,DEF456\n")
        expected_ids = {"EMD-1470", "EMD-1480"}
        ids = read_emdb_ids_from_csv(csv_file)
        assert ids == expected_ids

    def test_without_header(self, tmp_path: Path):
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("EMD-1470\nEMD-1480\n")
        expected_ids = {"EMD-1470", "EMD-1480"}
        ids = read_emdb_ids_from_csv(csv_file)
        assert ids == expected_ids
