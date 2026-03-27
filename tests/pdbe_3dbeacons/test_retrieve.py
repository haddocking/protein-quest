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


class TestRetrieveStructures:
    @pytest.mark.vcr
    @pytest.mark.asyncio
    async def test_defaults_to_cif_gz(self, tmp_path: Path):
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
    async def test_native_pdb_convert_to_cif_gz(self, tmp_path: Path):
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
    async def test_skips_not_found(self, tmp_path: Path):
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
    @pytest.mark.default_cassette("TestRetrieveStructures.test_native_pdb_convert_to_cif_gz.yaml")
    @pytest.mark.asyncio
    async def test_raw_keeps_native_extension(self, tmp_path: Path):
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
    @pytest.mark.asyncio
    async def test_uses_cacher_and_skips_download(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
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

    @pytest.mark.asyncio
    @pytest.mark.vcr
    @pytest.mark.default_cassette("TestRetrieveStructures.test_native_pdb_convert_to_cif_gz.yaml")
    async def test_caches_converted_output_for_later_run(self, tmp_path: Path):
        rows = [
            RetrieveStructureRow(
                provider="isoformio",
                model_identifier="CHS.58955.1",
                model_url="https://storage.googleapis.com/isoform.io/pdb/CHS.58955.1.pdb",
                model_format="PDB",
            )
        ]

        cacher = DirectoryCacher(tmp_path / "cache", copy_method="copy")
        first_output_dir = tmp_path / "first"
        second_output_dir = tmp_path / "second"
        # populate cache with first run
        expected_name = "isoformio~CHS.58955.1.cif.gz"
        await retrieve_structures(rows, first_output_dir, cacher=cacher)

        second_summary = await retrieve_structures(rows, second_output_dir, cacher=cacher)

        first_output_file = first_output_dir / expected_name
        second_output_file = second_output_dir / expected_name
        assert second_output_file.read_bytes() == first_output_file.read_bytes()
        assert second_summary == RetrieveStructureSummary(
            requested=0,
            downloaded=0,
            skipped=0,
            converted=0,
            final=0,
            cached=1,
        )

    @pytest.mark.asyncio
    @pytest.mark.vcr
    async def test_cif_to_gzip_cached(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        rows = [
            RetrieveStructureRow(
                provider="ped",
                model_identifier="PED03502e001",
                model_url="https://deposition.proteinensemble.org/api/v1/entries/PED03502/ensembles/e001/ensemble-sample/",
                model_format="MMCIF",
            )
        ]

        cacher = DirectoryCacher(tmp_path / "cache", copy_method="copy")
        output_dir = tmp_path / "second"
        # populate cache with first run
        expected_name = "ped~PED03502e001.cif.gz"

        with caplog.at_level("DEBUG", logger="protein_quest.pdbe_3dbeacons.retrieve"):
            summary = await retrieve_structures(rows, output_dir, cacher=cacher)

        output_file = output_dir / expected_name
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        assert_is_gzipped(output_file)
        assert summary == RetrieveStructureSummary(
            requested=1,
            downloaded=1,
            skipped=0,
            converted=1,
            final=1,
            cached=0,
        )

    @pytest.mark.vcr
    @pytest.mark.asyncio
    async def test_skips_download_when_output_exists(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
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

    @pytest.mark.asyncio
    async def test_raises_on_invalid_format(self, tmp_path: Path):
        rows = [
            RetrieveStructureRow(
                provider="swissmodel",
                model_identifier="42",
                model_url="http://example.org/model.42",
                model_format="INVALID_FORMAT",  # type: ignore[reportArgumentType]
            )
        ]

        with pytest.raises(ValueError, match="Unsupported model format: INVALID_FORMAT"):
            await retrieve_structures(rows, tmp_path)

    @pytest.mark.asyncio
    async def test_raises_on_invalid_provider(self, tmp_path: Path):
        rows = [
            RetrieveStructureRow(
                provider="invalid_provider",  # type: ignore[reportArgumentType]
                model_identifier="42",
                model_url="http://example.org/model.cif",
                model_format="MMCIF",
            )
        ]

        with pytest.raises(ValueError, match="Unsupported provider: invalid_provider"):
            await retrieve_structures(rows, tmp_path)

    @pytest.mark.asyncio
    async def test_raw_exists(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        rows = [
            RetrieveStructureRow(
                provider="swissmodel",
                model_identifier="Q9NTW7_329-603:5v3m.1.C",
                model_url="https://swissmodel.expasy.org/3d-beacons/uniprot/Q9NTW7.cif?range=329-603&template=5v3m.1.C&provider=swissmodel",
                model_format="MMCIF",
            )
        ]

        output_file = tmp_path / "swissmodel~Q9NTW7_329-603:5v3m.1.C.cif"
        output_file.write_bytes(b"existing-content")

        with caplog.at_level("DEBUG", logger="protein_quest.pdbe_3dbeacons.retrieve"):
            summary = await retrieve_structures(rows, tmp_path, raw=True)

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

    @pytest.mark.asyncio
    async def test_raw_exists_in_cache(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        rows = [
            RetrieveStructureRow(
                provider="swissmodel",
                model_identifier="Q9NTW7_329-603:5v3m.1.C",
                model_url="https://swissmodel.expasy.org/3d-beacons/uniprot/Q9NTW7.cif?range=329-603&template=5v3m.1.C&provider=swissmodel",
                model_format="MMCIF",
            )
        ]
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cacher = DirectoryCacher(cache_dir)
        filename = "swissmodel~Q9NTW7_329-603:5v3m.1.C.cif"
        cache_source_file = tmp_path / filename
        await cacher.write_bytes(cache_source_file, b"cached-content")

        output_dir = tmp_path / "output"

        with caplog.at_level("DEBUG", logger="protein_quest.pdbe_3dbeacons.retrieve"):
            summary = await retrieve_structures(rows, output_dir, cacher=cacher, raw=True)
        output_file = output_dir / filename
        assert output_file.exists()
        assert output_file.read_bytes() == b"cached-content"
        assert summary == RetrieveStructureSummary(
            requested=0,
            downloaded=0,
            skipped=0,
            converted=0,
            final=0,
            cached=1,
        )
        assert f"File {filename} already exists in cache, skipping download." in caplog.text


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
