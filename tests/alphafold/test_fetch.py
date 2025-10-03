from pathlib import Path

import pytest

from protein_quest.alphafold.fetch import fetch_many


@pytest.mark.vcr
def test_fetch_many(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    results = fetch_many(ids, tmp_path, {"summary", "pdb"})

    assert len(results) == 1
    fresult = results[0]
    assert fresult.uniprot_acc == theid
    assert fresult.summary is not None
    assert (tmp_path / f"{theid}.json").exists()
    assert fresult.pdb_file and fresult.pdb_file.exists()


@pytest.mark.vcr
def test_fetch_many_gzipped(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    results = fetch_many(ids, tmp_path, {"summary", "pdb", "cif"}, gzip_files=True)

    assert len(results) == 1
    fresult = results[0]
    assert fresult.uniprot_acc == theid
    assert fresult.summary is not None
    assert (tmp_path / f"{theid}.json").exists()
    assert fresult.pdb_file and fresult.pdb_file.exists()
    assert fresult.pdb_file.suffix == ".gz"
    assert fresult.cif_file and fresult.cif_file.exists()
    assert fresult.cif_file.suffix == ".gz"
    assert fresult.bcif_file is None
