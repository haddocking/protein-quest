from io import StringIO
from pathlib import Path
from textwrap import dedent

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
    def test_with_header(self):
        csv_data = dedent("""\
            emdb_id,uniprot_acc
            EMD-1470,ABC123
            EMD-1480,DEF456
            """)
        file = StringIO(csv_data)
        expected_ids = {"EMD-1470", "EMD-1480"}
        ids = read_emdb_ids_from_csv(file)
        assert ids == expected_ids

    def test_without_header(self):
        csv_data = dedent("""\
            EMD-1470
            EMD-1480
            """)
        file = StringIO(csv_data)
        expected_ids = {"EMD-1470", "EMD-1480"}
        ids = read_emdb_ids_from_csv(file)
        assert ids == expected_ids
