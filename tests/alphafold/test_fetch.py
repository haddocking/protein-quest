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
    assert fresult.uniprot_accession == theid
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
    assert fresult.uniprot_accession == theid
    assert fresult.summary is not None
    assert (tmp_path / f"{theid}.json").exists()
    assert fresult.pdb_file and fresult.pdb_file.exists()
    assert fresult.pdb_file.suffix == ".gz"
    assert fresult.cif_file and fresult.cif_file.exists()
    assert fresult.cif_file.suffix == ".gz"
    assert fresult.bcif_file is None


@pytest.mark.vcr
def test_fetch_many_all_isoforms(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    results = fetch_many(ids, tmp_path, {"summary"}, all_isoforms=True)

    # On https://www.uniprot.org/uniprotkb/P05067/entry#sequences
    # there are 11 isoforms.
    # Its P05067-3 isoform is on https://alphafold.ebi.ac.uk/entry/AF-P05067-3-F1
    # , but is not returned by the prediction API endpoint, so we expect 10 results here
    assert len(results) == 10
    assert all(result.uniprot_accession and result.uniprot_accession.startswith(theid) for result in results)
    canonical_results = [r for r in results if r.summary.uniprotAccession == theid]
    assert len(canonical_results) == 1
