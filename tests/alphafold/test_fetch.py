import re
from pathlib import Path

import pytest
from yarl import URL

from protein_quest.alphafold.entry_summary import EntrySummary
from protein_quest.alphafold.fetch import (
    AlphaFoldEntry,
    DownloadableFormat,
    downloadable_formats,
    fetch_alphafold_db_version,
    fetch_many,
    files_for_alphafold_entries,
)
from protein_quest.converter import converter


@pytest.fixture
def sample_entry_summary() -> EntrySummary:
    summary_json = '{"toolUsed":"AlphaFold Monomer v2.0 pipeline","providerId":"GDM","entityType":"protein","isUniProt":true,"modelEntityId":"AF-P50613-F1","modelCreatedDate":"2025-08-01T00:00:00Z","sequenceVersionDate":"1996-10-01T00:00:00Z","globalMetricValue":82.0,"fractionPlddtVeryLow":0.162,"fractionPlddtLow":0.092,"fractionPlddtConfident":0.127,"fractionPlddtVeryHigh":0.618,"latestVersion":6,"allVersions":[1,2,3,4,5,6],"sequence":"MALDVKSRAKRYEKLDFLGEGQFATVYKARDKNTNQIVAIKKIKLGHRSEAKDGINRTALREIKLLQELSHPNIIGLLDAFGHKSNISLVFDFMETDLEVIIKDNSLVLTPSHIKAYMLMTLQGLEYLHQHWILHRDLKPNNLLLDENGVLKLADFGLAKSFGSPNRAYTHQVVTRWYRAPELLFGARMYGVGVDMWAVGCILAELLLRVPFLPGDSDLDQLTRIFETLGTPTEEQWPDMCSLPDYVTFKSFPGIPLHHIFSAAGDDLLDLIQGLFLFNPCARITATQALKMKYFSNRPGPTPGCQLPRPNCPVETLKEQSNPALAIKRKRTEALEQGGLPKKLIF","sequenceStart":1,"sequenceEnd":346,"sequenceChecksum":"0A94BFA7DD416CEB","isUniProtReviewed":true,"gene":"CDK7","uniprotAccession":"P50613","uniprotId":"CDK7_HUMAN","uniprotDescription":"Cyclin-dependent kinase 7","taxId":9606,"organismScientificName":"Homo sapiens","isUniProtReferenceProteome":true,"bcifUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-model_v6.bcif","cifUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-model_v6.cif","pdbUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-model_v6.pdb","paeImageUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-predicted_aligned_error_v6.png","msaUrl":"https://alphafold.ebi.ac.uk/files/msa/AF-P50613-F1-msa_v6.a3m","plddtDocUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-confidence_v6.json","paeDocUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-predicted_aligned_error_v6.json","amAnnotationsUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-aa-substitutions.csv","amAnnotationsHg19Url":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-hg19.csv","amAnnotationsHg38Url":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-hg38.csv","entryId":"AF-P50613-F1","uniprotSequence":"MALDVKSRAKRYEKLDFLGEGQFATVYKARDKNTNQIVAIKKIKLGHRSEAKDGINRTALREIKLLQELSHPNIIGLLDAFGHKSNISLVFDFMETDLEVIIKDNSLVLTPSHIKAYMLMTLQGLEYLHQHWILHRDLKPNNLLLDENGVLKLADFGLAKSFGSPNRAYTHQVVTRWYRAPELLFGARMYGVGVDMWAVGCILAELLLRVPFLPGDSDLDQLTRIFETLGTPTEEQWPDMCSLPDYVTFKSFPGIPLHHIFSAAGDDLLDLIQGLFLFNPCARITATQALKMKYFSNRPGPTPGCQLPRPNCPVETLKEQSNPALAIKRKRTEALEQGGLPKKLIF","uniprotStart":1,"uniprotEnd":346,"isReferenceProteome":true,"isReviewed":true}'

    return converter.loads(summary_json, EntrySummary)


@pytest.fixture
def sample_entry(sample_entry_summary: EntrySummary) -> AlphaFoldEntry:
    return AlphaFoldEntry(
        uniprot_accession="P05067",
        summary=sample_entry_summary,
        summary_file=Path("P05067.json"),
        pdb_file=Path("P05067.pdb"),
        cif_file=Path("P05067.cif"),
        bcif_file=Path("P05067.bcif"),
        pae_doc_file=Path("P05067.pae.json"),
        msa_file=Path("P05067.msa"),
        am_annotations_file=Path("P05067.am_annotations.csv"),
        am_annotations_hg19_file=Path("P05067.am_annotations_hg19.csv"),
        am_annotations_hg38_file=Path("P05067.am_annotations_hg38.csv"),
        plddt_doc_file=Path("P05067.plddt.json"),
    )


class TestAlphaFoldEntry:
    @pytest.mark.parametrize("af_format", downloadable_formats)
    def test_by_format(self, af_format: DownloadableFormat, sample_entry: AlphaFoldEntry):
        file_path = sample_entry.by_format(af_format)
        assert file_path is not None

    def test_by_format_bad_format(self, sample_entry: AlphaFoldEntry):
        with pytest.raises(ValueError, match="Invalid format"):
            sample_entry.by_format("bad_format")  # type: ignore  # noqa: PGH003

    def test_nr_of_files_none(self):
        entry = AlphaFoldEntry(
            uniprot_accession="P05067",
        )
        assert entry.nr_of_files() == 0

    def test_nr_of_files_all(self, sample_entry: AlphaFoldEntry):
        assert sample_entry.nr_of_files() == 10

    def test_relative_to(self, sample_entry: AlphaFoldEntry, tmp_path: Path):
        session_dir = tmp_path / "session1"
        session_dir.mkdir()
        abs_entry = AlphaFoldEntry(
            uniprot_accession="P05067",
            summary=sample_entry.summary,
            summary_file=session_dir / "P05067.json",
            pdb_file=session_dir / "P05067.pdb",
            cif_file=session_dir / "P05067.cif",
            bcif_file=session_dir / "P05067.bcif",
            pae_doc_file=session_dir / "P05067.pae.json",
            msa_file=session_dir / "P05067.msa",
            am_annotations_file=session_dir / "P05067.am_annotations.csv",
            am_annotations_hg19_file=session_dir / "P05067.am_annotations_hg19.csv",
            am_annotations_hg38_file=session_dir / "P05067.am_annotations_hg38.csv",
            plddt_doc_file=session_dir / "P05067.plddt.json",
        )

        relative_entry = abs_entry.relative_to(session_dir)

        assert relative_entry == sample_entry


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
    canonical_results = [r for r in results if r.summary is not None and r.summary.uniprotAccession == theid]
    assert len(canonical_results) == 1


@pytest.mark.vcr
def test_fetch_many_no_summary(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    results = fetch_many(ids, tmp_path, {"cif"}, gzip_files=False)

    assert len(results) == 1
    fresult = results[0]
    fexpected = AlphaFoldEntry(
        uniprot_accession=theid,
        cif_file=tmp_path / "AF-P05067-F1-model_v6.cif",
    )
    assert fresult == fexpected


@pytest.mark.vcr
def test_fetch_many_no_summary_with_version(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    results = fetch_many(ids, tmp_path, {"cif"}, db_version="6", gzip_files=False)

    assert len(results) == 1
    fresult = results[0]
    fexpected = AlphaFoldEntry(
        uniprot_accession=theid,
        cif_file=tmp_path / "AF-P05067-F1-model_v6.cif",
    )
    assert fresult == fexpected


def test_fetch_many_no_formats(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    with pytest.raises(ValueError, match="At least one format must be specified"):
        fetch_many(ids, tmp_path, set())


def test_fetch_many_all_isoforms_no_summary(tmp_path: Path):
    theid = "P05067"
    ids = [theid]

    with pytest.raises(
        ValueError, match=re.escape("Cannot fetch all isoforms when 'summary' is not in 'formats' set.")
    ):
        fetch_many(ids, tmp_path, {"cif"}, all_isoforms=True)


@pytest.mark.vcr
@pytest.mark.asyncio
async def test_fetch_alphafold_db_version():
    result = await fetch_alphafold_db_version()

    expected = "6"
    assert result == expected


def test_files_for_alphafold_entries():
    db_version = "6"
    uniprot_accessions = ["P05067", "Q9H9K5"]
    formats = downloadable_formats
    result = files_for_alphafold_entries(uniprot_accessions, formats, db_version, gzip_files=False)

    expected = {
        "P05067": {
            "amAnnotations": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-aa-substitutions.csv"),
                "AF-P05067-F1-aa-substitutions.csv",
            ),
            "amAnnotationsHg19": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-hg19.csv"),
                "AF-P05067-F1-hg19.csv",
            ),
            "amAnnotationsHg38": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-hg38.csv"),
                "AF-P05067-F1-hg38.csv",
            ),
            "bcif": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-model_v6.bcif"),
                "AF-P05067-F1-model_v6.bcif",
            ),
            "cif": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-model_v6.cif"),
                "AF-P05067-F1-model_v6.cif",
            ),
            "msa": (
                URL("https://alphafold.ebi.ac.uk/files/msa/AF-P05067-F1-msa_v6.a3m"),
                "AF-P05067-F1-msa_v6.a3m",
            ),
            "paeDoc": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-predicted_aligned_error_v6.json"),
                "AF-P05067-F1-predicted_aligned_error_v6.json",
            ),
            "pdb": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-model_v6.pdb"),
                "AF-P05067-F1-model_v6.pdb",
            ),
            "plddtDoc": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-confidence_v6.json"),
                "AF-P05067-F1-confidence_v6.json",
            ),
        },
        "Q9H9K5": {
            "amAnnotations": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-aa-substitutions.csv"),
                "AF-Q9H9K5-F1-aa-substitutions.csv",
            ),
            "amAnnotationsHg19": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-hg19.csv"),
                "AF-Q9H9K5-F1-hg19.csv",
            ),
            "amAnnotationsHg38": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-hg38.csv"),
                "AF-Q9H9K5-F1-hg38.csv",
            ),
            "bcif": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-model_v6.bcif"),
                "AF-Q9H9K5-F1-model_v6.bcif",
            ),
            "cif": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-model_v6.cif"),
                "AF-Q9H9K5-F1-model_v6.cif",
            ),
            "msa": (
                URL("https://alphafold.ebi.ac.uk/files/msa/AF-Q9H9K5-F1-msa_v6.a3m"),
                "AF-Q9H9K5-F1-msa_v6.a3m",
            ),
            "paeDoc": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-predicted_aligned_error_v6.json"),
                "AF-Q9H9K5-F1-predicted_aligned_error_v6.json",
            ),
            "pdb": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-model_v6.pdb"),
                "AF-Q9H9K5-F1-model_v6.pdb",
            ),
            "plddtDoc": (
                URL("https://alphafold.ebi.ac.uk/files/AF-Q9H9K5-F1-confidence_v6.json"),
                "AF-Q9H9K5-F1-confidence_v6.json",
            ),
        },
    }
    assert result == expected


def test_files_for_alphafold_entries_someformats_gzipped():
    db_version = "6"
    uniprot_accessions = ["P05067"]
    formats: set[DownloadableFormat] = {"pdb", "cif"}
    result = files_for_alphafold_entries(uniprot_accessions, formats, db_version, gzip_files=True)

    expected = {
        "P05067": {
            "cif": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-model_v6.cif"),
                "AF-P05067-F1-model_v6.cif.gz",
            ),
            "pdb": (
                URL("https://alphafold.ebi.ac.uk/files/AF-P05067-F1-model_v6.pdb"),
                "AF-P05067-F1-model_v6.pdb.gz",
            ),
        },
    }
    assert result == expected
