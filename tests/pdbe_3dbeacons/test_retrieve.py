import gzip
from io import StringIO
from pathlib import Path

import pytest
from cattrs import ClassValidationError

from protein_quest.pdbe_3dbeacons.retrieve import (
    RetrieveStructureRow,
    RetrieveStructureSummary,
    read_retrieve_structure_rows,
    retrieve_structures,
)
from protein_quest.utils import DirectoryCacher


def assert_is_gzipped(path: Path) -> None:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        assert handle.read(1)


@pytest.mark.vcr
@pytest.mark.asyncio
async def test_retrieve_structures_defaults_to_cif_gz(tmp_path: Path):
    rows = [
        RetrieveStructureRow(
            provider="swissmodel",
            model_identifier="Q9NTW7_329-603:5v3m.1.C",
            model_url="https://swissmodel.expasy.org/3d-beacons/uniprot/Q9NTW7.cif?range=329-603&template=5v3m.1.C&provider=swissmodel",
            model_format="MMCIF",
        )
    ]

    summary = await retrieve_structures(rows, tmp_path)

    output_file = tmp_path / "swissmodel~Q9NTW7_329-603:5v3m.1.C.cif.gz"
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    assert_is_gzipped(output_file)

    assert summary == RetrieveStructureSummary(requested=1, downloaded=1, skipped=0, converted=0, final=1)


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_retrieve_structures_native_pdb_convert_to_cif_gz(tmp_path: Path):
    rows = [
        RetrieveStructureRow(
            provider="isoformio",
            model_identifier="CHS.58955.1",
            model_url="https://storage.googleapis.com/isoform.io/pdb/CHS.58955.1.pdb",
            model_format="PDB",
        )
    ]

    summary = await retrieve_structures(rows, tmp_path)

    output_file = tmp_path / "isoformio~CHS.58955.1.cif.gz"
    source_file = tmp_path / "isoformio~CHS.58955.1.pdb"
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    assert_is_gzipped(output_file)
    assert not source_file.exists()
    assert summary == RetrieveStructureSummary(requested=1, downloaded=1, skipped=0, converted=1, final=1)


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_retrieve_structures_skips_not_found(tmp_path: Path):
    rows = [
        RetrieveStructureRow(
            provider="pdbe",
            model_identifier="not-a-real-entry",
            model_url="https://www.ebi.ac.uk/pdbe/static/entry/not-a-real-entry_updated.cif",
            model_format="MMCIF",
        )
    ]

    summary = await retrieve_structures(rows, tmp_path)

    assert summary == RetrieveStructureSummary(requested=1, downloaded=0, skipped=1, converted=0, final=0)
    assert not any(tmp_path.iterdir())


@pytest.mark.vcr
@pytest.mark.default_cassette("test_retrieve_structures_native_pdb_convert_to_cif_gz.yaml")
@pytest.mark.asyncio
async def test_retrieve_structures_raw_keeps_native_extension(tmp_path: Path):
    rows = [
        RetrieveStructureRow(
            provider="isoformio",
            model_identifier="CHS.58955.1",
            model_url="https://storage.googleapis.com/isoform.io/pdb/CHS.58955.1.pdb",
            model_format="PDB",
        )
    ]

    summary = await retrieve_structures(rows, tmp_path, raw=True)

    output_file = tmp_path / "isoformio~CHS.58955.1.pdb"
    converted_target = tmp_path / "isoformio~CHS.58955.1.cif.gz"
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    assert not converted_target.exists()
    assert summary == RetrieveStructureSummary(requested=1, downloaded=1, skipped=0, converted=0, final=1)


@pytest.mark.vcr
@pytest.mark.default_cassette("test_retrieve_structures_native_pdb_convert_to_cif_gz.yaml")
@pytest.mark.asyncio
async def test_retrieve_structures_uses_cacher_and_skips_download(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    rows = [
        RetrieveStructureRow(
            provider="isoformio",
            model_identifier="CHS.58955.1",
            model_url="https://storage.googleapis.com/isoform.io/pdb/CHS.58955.1.pdb",
            model_format="PDB",
        )
    ]

    cacher = DirectoryCacher(tmp_path / "cache", copy_method="copy")
    cached_file_name = "isoformio~CHS.58955.1.cif.gz"
    cached_source = tmp_path / "seed" / cached_file_name
    cached_source.parent.mkdir(parents=True, exist_ok=True)
    await cacher.write_bytes(cached_source, b"cached-gz-content")

    output_dir = tmp_path
    with caplog.at_level("DEBUG", logger="protein_quest.pdbe_3dbeacons.retrieve"):
        summary = await retrieve_structures(rows, output_dir, cacher=cacher)

    output_file = output_dir / cached_file_name
    assert output_file.exists()
    assert output_file.read_bytes() == b"cached-gz-content"
    assert summary == RetrieveStructureSummary(
        requested=0,
        downloaded=0,
        skipped=0,
        converted=0,
        final=0,
        cached=1,
    )
    assert f"File {cached_file_name} already exists in cache, skipping download." in caplog.text


@pytest.mark.vcr
@pytest.mark.default_cassette("test_retrieve_structures_defaults_to_cif_gz.yaml")
@pytest.mark.asyncio
async def test_retrieve_structures_skips_download_when_output_exists(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    rows = [
        RetrieveStructureRow(
            provider="swissmodel",
            model_identifier="Q9NTW7_329-603:5v3m.1.C",
            model_url="https://swissmodel.expasy.org/3d-beacons/uniprot/Q9NTW7.cif?range=329-603&template=5v3m.1.C&provider=swissmodel",
            model_format="MMCIF",
        )
    ]

    output_file = tmp_path / "swissmodel~Q9NTW7_329-603:5v3m.1.C.cif.gz"
    output_file.write_bytes(b"existing-content")

    with caplog.at_level("DEBUG", logger="protein_quest.pdbe_3dbeacons.retrieve"):
        summary = await retrieve_structures(rows, tmp_path)

    assert output_file.read_bytes() == b"existing-content"
    assert summary == RetrieveStructureSummary(
        requested=0,
        downloaded=0,
        skipped=0,
        converted=0,
        final=0,
        cached=0,
    )
    assert f"File {output_file.name} already exists in {tmp_path}, skipping download." in caplog.text


class TestReadRetrieveStructureRows:
    def test_allows_extra_columns_and_strips_values(self):
        csv_content = StringIO(
            "provider,model_identifier,model_url,model_format,uniprot_accession\n"
            "swissmodel,Q9NTW7_329-603:5v3m.1.C,https://example.org/model.cif,MMCIF,Q9NTW7\n"
        )

        rows = read_retrieve_structure_rows(csv_content)

        assert rows == [
            RetrieveStructureRow(
                provider="swissmodel",
                model_identifier="Q9NTW7_329-603:5v3m.1.C",
                model_url="https://example.org/model.cif",
                model_format="MMCIF",
            )
        ]

    def test_fails_on_missing_provider_field(self):
        csv_content = StringIO(
            "provider,model_identifier,model_url,model_format,uniprot_accession\n"
            ",Q9NTW7_329-603:5v3m.1.C,https://example.org/model.cif,MMCIF,Q9NTW7\n"
        )

        with pytest.raises(ClassValidationError):
            read_retrieve_structure_rows(csv_content)
