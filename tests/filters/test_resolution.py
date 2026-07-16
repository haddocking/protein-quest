import csv
import gzip
from dataclasses import replace
from pathlib import Path

import pytest
from cyclopts.types import StdioPath

from protein_quest.errors import ResolutionUnsetError
from protein_quest.filters.resolution import (
    NoUniProtAccessionError,
    OutsideTopError,
    ResolutionFilterStatistics,
    SequenceIdentityBelowThresholdError,
    copy_resolution_statistics,
    filter_files_on_resolution,
    filter_on_sequence_identity,
    load_resolution_statistics,
    sort_resolution_statistics,
    write_resolution_stats,
)


def _make_stats(
    filename: str,
    accession: str | None,
    resolution: float,
    pdb_id: str | None = None,
    total_residue_count: int = 100,
    is_alphafold: bool = False,
    uniprot_start: int = 0,
    uniprot_end: int = 0,
    sequence_identity: float = 0.0,
    chain_length: int = 0,
    discard_reason: Exception | None = None,
    passed: bool = False,
    output_file: Path | None = None,
) -> ResolutionFilterStatistics:
    if pdb_id is None:
        pdb_id = filename
    return ResolutionFilterStatistics(
        id=pdb_id,
        input_file=Path(f"/fake/{filename}"),
        uniprot_accession=accession,
        chain_id=None,
        resolution=resolution,
        total_residue_count=total_residue_count,
        is_alphafold=is_alphafold,
        uniprot_start=uniprot_start,
        uniprot_end=uniprot_end,
        sequence_identity=sequence_identity,
        chain_length=chain_length,
        passed=passed,
        output_file=output_file,
        discard_reason=discard_reason,
    )


class TestResolutionFilterStatistics:
    def test_is_hashable(self):
        stats = _make_stats("a.cif.gz", "P12345", resolution=1.0)
        assert isinstance(hash(stats), int)

    def test_passed_equality(self):
        stats1 = ResolutionFilterStatistics(
            id="a.cif.gz",
            input_file=Path("/fake/a.cif.gz"),
            uniprot_accession="P12345",
            chain_id=None,
            resolution=1.0,
            total_residue_count=100,
            is_alphafold=False,
            uniprot_start=1,
            uniprot_end=10,
            sequence_identity=1.0,
            chain_length=10,
            passed=True,
            output_file=Path("/fake/output/a.cif.gz"),
            discard_reason=None,
        )
        stats2 = ResolutionFilterStatistics(
            id="a.cif.gz",
            input_file=Path("/fake/a.cif.gz"),
            uniprot_accession="P12345",
            chain_id=None,
            resolution=1.0,
            total_residue_count=100,
            is_alphafold=False,
            uniprot_start=1,
            uniprot_end=10,
            sequence_identity=1.0,
            chain_length=10,
            passed=True,
            output_file=Path("/fake/output/a.cif.gz"),
            discard_reason=None,
        )
        assert stats1 == stats2
        assert hash(stats1) == hash(stats2)

    def test_discarded_equality(self):
        stats1 = ResolutionFilterStatistics(
            id="a.cif.gz",
            input_file=Path("/fake/a.cif.gz"),
            uniprot_accession="P12345",
            chain_id=None,
            resolution=1.0,
            total_residue_count=100,
            is_alphafold=True,
            uniprot_start=1,
            uniprot_end=10,
            sequence_identity=1.0,
            chain_length=10,
            passed=False,
            output_file=Path("/fake/output/a.cif.gz"),
            discard_reason=ValueError("Some reason"),
        )
        stats2 = ResolutionFilterStatistics(
            id="a.cif.gz",
            input_file=Path("/fake/a.cif.gz"),
            uniprot_accession="P12345",
            chain_id=None,
            resolution=1.0,
            total_residue_count=100,
            is_alphafold=True,
            uniprot_start=1,
            uniprot_end=10,
            sequence_identity=1.0,
            chain_length=10,
            passed=False,
            output_file=Path("/fake/output/a.cif.gz"),
            discard_reason=ValueError("Some reason"),
        )
        assert stats1 == stats2
        assert hash(stats1) == hash(stats2)


class TestYieldResolutionStatistics:
    def test_metadata_in_order(self, sample_cif: Path, sample2_cif: Path, af_cif: Path, nmr_cif: Path):
        input_files = [sample_cif, sample2_cif, af_cif, nmr_cif]
        results = load_resolution_statistics(input_files, scheduler_address="sequential")

        expected = [
            ResolutionFilterStatistics(
                id="3JRS",
                input_file=sample_cif,
                uniprot_accession="Q8VZS8",
                chain_id='A',
                resolution=2.05,
                total_residue_count=173,
                is_alphafold=False,
                uniprot_start=8,
                uniprot_end=211,
                sequence_identity=1.0,
                chain_length=173,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                id="2Y29",
                input_file=sample2_cif,
                uniprot_accession="P05067",
                chain_id='A',
                resolution=2.3,
                total_residue_count=8,
                is_alphafold=False,
                uniprot_start=687,
                uniprot_end=692,
                sequence_identity=1.0,
                chain_length=8,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                id="AF-A0A0C5B5G6-F1",
                input_file=af_cif,
                uniprot_accession="A0A0C5B5G6",
                chain_id='A',
                resolution=0.0,
                total_residue_count=16,
                is_alphafold=True,
                uniprot_start=1,
                uniprot_end=16,
                sequence_identity=1.0,
                chain_length=16,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                id="1AMB",
                input_file=nmr_cif,
                uniprot_accession="P05067",
                chain_id='A',
                resolution=0.0,
                total_residue_count=28,
                is_alphafold=False,
                uniprot_start=672,
                uniprot_end=699,
                sequence_identity=1.0,
                chain_length=28,
                passed=False,
                output_file=None,
            ),
        ]
        assert results == expected

    def test_empty_input(self):
        assert load_resolution_statistics([], scheduler_address="sequential") == []


class TestLoadResolutionStatistics:
    @pytest.mark.parametrize("scheduler_address", [None, "sequential"])
    def test_metadata_in_order(
        self,
        sample_cif: Path,
        sample2_cif: Path,
        af_cif: Path,
        nmr_cif: Path,
        scheduler_address: str | None,
    ):
        input_files = [sample_cif, sample2_cif, af_cif, nmr_cif]
        results = load_resolution_statistics(input_files, scheduler_address=scheduler_address)

        expected = [
            ResolutionFilterStatistics(
                id="3JRS",
                input_file=sample_cif,
                uniprot_accession="Q8VZS8",
                chain_id='A',
                resolution=2.05,
                total_residue_count=173,
                is_alphafold=False,
                uniprot_start=8,
                uniprot_end=211,
                sequence_identity=1.0,
                chain_length=173,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                id="2Y29",
                input_file=sample2_cif,
                uniprot_accession="P05067",
                chain_id=None,
                resolution=2.3,
                total_residue_count=8,
                is_alphafold=False,
                uniprot_start=687,
                uniprot_end=692,
                sequence_identity=1.0,
                chain_length=8,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                id="AF-A0A0C5B5G6-F1",
                input_file=af_cif,
                uniprot_accession="A0A0C5B5G6",
                chain_id=None,
                resolution=0.0,
                total_residue_count=16,
                is_alphafold=True,
                uniprot_start=1,
                uniprot_end=16,
                sequence_identity=1.0,
                chain_length=16,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                id="1AMB",
                input_file=nmr_cif,
                uniprot_accession="P05067",
                chain_id=None,
                resolution=0.0,
                total_residue_count=28,
                is_alphafold=False,
                uniprot_start=672,
                uniprot_end=699,
                sequence_identity=1.0,
                chain_length=28,
                passed=False,
                output_file=None,
            ),
        ]

        assert results == expected

    def test_empty_input(self):
        assert load_resolution_statistics([], scheduler_address="sequential") == []

    def test_multiple_accessions_uniprotlessresult(self, multi_accession_cif: Path):
        result = load_resolution_statistics([multi_accession_cif], scheduler_address="sequential")

        expected = ResolutionFilterStatistics(
            id="1A02",
            input_file=multi_accession_cif,
            uniprot_accession=None,
            chain_id='A',
            resolution=2.7,
            total_residue_count=513,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=40,
            passed=False,
            output_file=None,
        )
        assert result == [expected]

    def test_invalid_file_sets_discard_reason(self, tmp_path: Path):
        invalid_file = tmp_path / "invalid.cif.gz"
        invalid_file.write_text("not a valid structure file")

        results = load_resolution_statistics([invalid_file], scheduler_address="sequential")

        expected = ResolutionFilterStatistics(
            id="invalid.cif",
            input_file=invalid_file,
            uniprot_accession=None,
            chain_id=None,
            resolution=0.0,
            total_residue_count=0,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=0,
            passed=False,
            output_file=None,
            discard_reason=RuntimeError(f"File not in the gzip format: {invalid_file}"),
        )
        assert results == [expected]


class TestFilterFilesOnResolution:
    def test_strict_no_uniprot_does_not_pass(self, sample2_cif: Path, tmp_path: Path):
        no_uniprot = tmp_path / "no-uniprot.cif.gz"
        with gzip.open(sample2_cif, "rt", encoding="utf-8") as handle:
            text = handle.read()
        # Remove UniProt database entries while keeping a valid mmCIF structure.
        no_uniprot.write_bytes(gzip.compress(text.replace("UNP", "XXX").encode("utf-8")))

        output_dir = tmp_path / "output"
        results = list(
            filter_files_on_resolution(
                input_files=[no_uniprot],
                output_dir=output_dir,
                top=1,
                min_sequence_identity=0.0,
                scheduler_address="sequential",
                lax=False,
            )
        )

        expected = ResolutionFilterStatistics(
            id="2Y29",
            input_file=no_uniprot,
            uniprot_accession=None,
            chain_id='A',
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=8,
            passed=False,
            output_file=None,
            discard_reason=NoUniProtAccessionError(no_uniprot),
        )
        assert results == [expected]
        assert output_dir.exists()
        assert list(output_dir.iterdir()) == []

    def test_no_uniprot_can_pass_when_grouping_disabled(self, no_uniprot_cif: Path, tmp_path: Path):
        output_dir = tmp_path / "output"
        results = list(
            filter_files_on_resolution(
                input_files=[no_uniprot_cif], output_dir=output_dir, top=1, group_by=False, min_sequence_identity=0.0
            )
        )

        expected = ResolutionFilterStatistics(
            id="2Y29",
            input_file=no_uniprot_cif,
            uniprot_accession=None,
            chain_id='A',
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=8,
            passed=True,
            output_file=output_dir / no_uniprot_cif.name,
        )
        assert results == [expected]
        assert results[0].output_file is not None and results[0].output_file.exists()

    def test_lax_mode_passes_structure_without_resolution(self, af_cif: Path, tmp_path: Path):
        output_dir = tmp_path / "output"
        results = list(
            filter_files_on_resolution(
                input_files=[af_cif], output_dir=output_dir, top=1, lax=True, copy_method="symlink"
            )
        )

        expected = ResolutionFilterStatistics(
            id="AF-A0A0C5B5G6-F1",
            input_file=af_cif,
            uniprot_accession="A0A0C5B5G6",
            chain_id='A',
            resolution=0.0,
            total_residue_count=16,
            is_alphafold=True,
            uniprot_start=1,
            uniprot_end=16,
            sequence_identity=1.0,
            chain_length=16,
            passed=True,
            output_file=output_dir / af_cif.name,
            discard_reason=ResolutionUnsetError(af_cif),
        )
        assert results == [expected]
        assert results[0].output_file is not None and results[0].output_file.exists()

    def test_lax_mode_passes_structure_without_uniprot_accession(self, no_uniprot_cif: Path):
        output_dir = no_uniprot_cif.parent / "output"
        results = list(
            filter_files_on_resolution(
                input_files=[no_uniprot_cif], output_dir=output_dir, top=1, lax=True, min_sequence_identity=0.0
            )
        )

        expected = ResolutionFilterStatistics(
            id="2Y29",
            input_file=no_uniprot_cif,
            uniprot_accession=None,
            chain_id='A',
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=8,
            passed=True,
            output_file=output_dir / no_uniprot_cif.name,
            discard_reason=NoUniProtAccessionError(no_uniprot_cif),
        )
        assert results == [expected]
        assert results[0].output_file is not None and results[0].output_file.exists()


class TestSortResolutionStatisticsYesGroupByNoCoverage:
    def test_groupby_accession(self):
        a1 = _make_stats("a1.cif.gz", "P11111", resolution=1.0)
        a2 = _make_stats("a2.cif.gz", "P11111", resolution=2.0)
        b1 = _make_stats("b1.cif.gz", "P22222", resolution=1.5)
        b2 = _make_stats("b2.cif.gz", "P22222", resolution=3.0)
        expected = [
            replace(a1, passed=True),
            replace(a2, passed=False, discard_reason=OutsideTopError(top=1, rank=2)),
            replace(b1, passed=True),
            replace(b2, passed=False, discard_reason=OutsideTopError(top=1, rank=2)),
        ]

        results = sort_resolution_statistics([a2, a1, b2, b1], top=1, coverage=False, group_by=True)

        assert results == expected

    def test_top_limit(self):
        a = _make_stats("a.cif.gz", "P12345", resolution=1.0)
        b = _make_stats("b.cif.gz", "P12345", resolution=2.0)
        c = _make_stats("c.cif.gz", "P12345", resolution=3.0)

        results = sort_resolution_statistics([a, b, c], top=2, coverage=False, group_by=True)

        passed_names = {r.input_file.name for r in results if r.passed}
        assert passed_names == {"a.cif.gz", "b.cif.gz"}

    def test_results_sorted_by_resolution_within_group(self):
        a = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        b = _make_stats("b.cif.gz", "P12345", resolution=3.0)
        c = _make_stats("c.cif.gz", "P12345", resolution=1.0)

        results = sort_resolution_statistics([a, b, c], top=10, coverage=False, group_by=True)

        assert [r.input_file.name for r in results] == ["c.cif.gz", "a.cif.gz", "b.cif.gz"]


class TestSortResolutionStatisticsNoGroupbyNoCoverage:
    def test_sorts_by_resolution(self):
        a = _make_stats("a.cif.gz", "P11111", resolution=1.0)  # best
        b = _make_stats("b.cif.gz", "P22222", resolution=1.5)
        c = _make_stats("c.cif.gz", "P33333", resolution=2.0)  # worst
        expected = [
            replace(a, passed=True),
            replace(b, passed=False, discard_reason=OutsideTopError(top=1, rank=2)),
            replace(c, passed=False, discard_reason=OutsideTopError(top=1, rank=3)),
        ]

        results = sort_resolution_statistics([c, a, b], top=1, coverage=False, group_by=False)

        assert results == expected

    def test_missing_accession_ignored(self):
        with_acc = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        no_acc = _make_stats("z.cif.gz", None, resolution=1.0)  # best
        expected = [
            replace(no_acc, passed=True),
            replace(with_acc, passed=False, discard_reason=OutsideTopError(top=1, rank=2)),
        ]

        results = sort_resolution_statistics([with_acc, no_acc], top=1, coverage=False, group_by=False)

        assert results == expected


class TestSortResolutionStatisticsNoGroupbyYesCoverage:
    def test_kitchensink(self):
        p1a = _make_stats(
            "p1a.cif.gz",
            "P11111",
            resolution=1.0,
            sequence_identity=0.95,
            chain_length=100,
            uniprot_start=1,
            uniprot_end=10,
        )
        p1b = _make_stats(
            "p1b.cif.gz",
            "P11111",
            resolution=1.5,
            sequence_identity=0.80,
            chain_length=80,
            uniprot_start=20,
            uniprot_end=30,
        )
        p2a = _make_stats(
            "p2a.cif.gz",
            "P22222",
            resolution=1.2,
            sequence_identity=0.90,
            chain_length=90,
            uniprot_start=1,
            uniprot_end=10,
        )
        p2b = _make_stats(
            "p2b.cif.gz",
            "P22222",
            resolution=1.7,
            sequence_identity=0.70,
            chain_length=70,
            uniprot_start=20,
            uniprot_end=30,
        )
        expected = [
            replace(p1a, passed=True),
            replace(p2a, passed=True),
            replace(p1b, passed=False, discard_reason=OutsideTopError(top=2, rank=3)),
            replace(p2b, passed=False, discard_reason=OutsideTopError(top=2, rank=4)),
        ]

        results = sort_resolution_statistics([p2b, p1b, p2a, p1a], top=2, coverage=True, group_by=False)

        assert results == expected


class TestSortResolutionStatisticsYesGroupbyYesCoverage:
    def test_kitchensink(self):
        p1a = _make_stats(
            "p1a.cif.gz",
            "P11111",
            resolution=1.0,
            sequence_identity=1.00,
            chain_length=100,
            uniprot_start=1,
            uniprot_end=100,
        )
        p1b = _make_stats(
            "p1b.cif.gz",
            "P11111",
            resolution=1.5,
            sequence_identity=1.00,
            chain_length=80,
            uniprot_start=1,
            uniprot_end=100,
        )
        p1c = _make_stats(
            "p1c.cif.gz",
            "P11111",
            resolution=1.8,  # worst resolution of P11111, last in group
            sequence_identity=1.00,
            chain_length=100,
            uniprot_start=1,
            uniprot_end=100,
        )
        p2a = _make_stats(
            "p2a.cif.gz",
            "P22222",
            resolution=1.2,
            sequence_identity=1.00,
            chain_length=100,
            uniprot_start=1,
            uniprot_end=100,
        )
        p2b = _make_stats(
            "p2b.cif.gz",
            "P22222",
            resolution=1.7,
            sequence_identity=1.00,
            chain_length=100,
            uniprot_start=1,
            uniprot_end=100,
        )
        expected = [
            replace(p1a, passed=True),
            replace(p1b, passed=True),
            replace(p1c, passed=False, discard_reason=OutsideTopError(top=2, rank=3)),
            replace(p2a, passed=True),
            replace(p2b, passed=True),
        ]

        results = sort_resolution_statistics([p2b, p1b, p1c, p2a, p1a], top=2, coverage=True, group_by=True)

        assert results == expected


class TestCopyResolutionStatistics:
    def test_passed_is_copied_and_output_set(self, sample2_cif: Path, tmp_path: Path):
        stats = ResolutionFilterStatistics(
            id="2Y29",
            input_file=sample2_cif,
            uniprot_accession="P05067",
            chain_id=None,
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=687,
            uniprot_end=692,
            sequence_identity=1.0,
            chain_length=8,
            passed=True,
            output_file=None,
        )
        output_dir = tmp_path / "output"
        expected = replace(stats, output_file=output_dir / sample2_cif.name)

        results = list(copy_resolution_statistics([stats], output_dir))

        assert results == [expected]
        assert results[0].output_file is not None and results[0].output_file.exists()

    def test_not_passed_is_unchanged(self, sample2_cif: Path, tmp_path: Path):
        stats = ResolutionFilterStatistics(
            id="2Y29",
            input_file=sample2_cif,
            uniprot_accession="P05067",
            chain_id=None,
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=687,
            uniprot_end=692,
            sequence_identity=1.0,
            chain_length=8,
            passed=False,
            output_file=None,
        )
        output_dir = tmp_path / "output"

        results = list(copy_resolution_statistics([stats], output_dir))

        assert results[0].output_file is None
        assert not any(f for f in output_dir.iterdir() if f.is_file())


def test_filter_on_sequence_identity(caplog: pytest.LogCaptureFixture):
    good = _make_stats("good.cif.gz", "P12345", resolution=1.0, sequence_identity=0.95)
    bad = _make_stats("bad.cif.gz", "P12345", resolution=1.0, sequence_identity=0.5)
    expected = [
        replace(good, passed=False, discard_reason=None),
        replace(
            bad,
            passed=False,
            discard_reason=SequenceIdentityBelowThresholdError(bad.input_file, bad.sequence_identity, 0.9),
        ),
    ]

    results = list(filter_on_sequence_identity(0.9, [good, bad]))

    assert results == expected
    assert (
        "Discarding /fake/bad.cif.gz due to sequence identity 0.500 below minimal sequence identity 0.900"
        in caplog.text
    )


def test_write_stats_with_edge_cases(tmp_path: StdioPath):
    """Test writing stats with various edge cases in discard_reason."""
    stats = [
        _make_stats("a.cif.gz", "P12345", resolution=1.0, passed=True),
        _make_stats(
            "b.cif.gz",
            "P12345",
            resolution=2.0,
            passed=False,
            discard_reason=ValueError('Message with "quotes" inside'),
        ),
        _make_stats(
            "c.cif.gz",
            "P12345",
            resolution=3.0,
            passed=False,
            discard_reason=ValueError("Line 1\nLine 2"),
        ),
        _make_stats(
            "d.cif.gz",
            "P12345",
            resolution=4.0,
            passed=False,
            discard_reason=ValueError("value1, value2, value3"),
        ),
        _make_stats(
            "e.cif.gz",
            "P12345",
            resolution=5.0,
            passed=False,
            discard_reason=ValueError('Error: "invalid" data\nat line 42'),
        ),
    ]
    output = tmp_path / "stats.csv"

    write_resolution_stats(stats, output)

    with output.open() as f:
        rows = list(csv.DictReader(f))

    expected = [
        {
            "uniprot_accession": "P12345",
            "chain_id": "",
            "pdb_id": "a.cif.gz",
            "input_file": "/fake/a.cif.gz",
            "resolution": "1.0",
            "total_residue_count": "100",
            "is_alphafold": "False",
            "uniprot_start": "0",
            "uniprot_end": "0",
            "sequence_identity": "0.000",
            "chain_length": "0",
            "passed": "True",
            "output_file": "",
            "discard_reason": "",
            "discard_reason_type": "",
        },
        {
            "uniprot_accession": "P12345",
            "chain_id": "",
            "pdb_id": "b.cif.gz",
            "input_file": "/fake/b.cif.gz",
            "resolution": "2.0",
            "total_residue_count": "100",
            "is_alphafold": "False",
            "uniprot_start": "0",
            "uniprot_end": "0",
            "sequence_identity": "0.000",
            "chain_length": "0",
            "passed": "False",
            "output_file": "",
            "discard_reason": 'Message with "quotes" inside',
            "discard_reason_type": "ValueError",
        },
        {
            "uniprot_accession": "P12345",
            "chain_id": "",
            "pdb_id": "c.cif.gz",
            "input_file": "/fake/c.cif.gz",
            "resolution": "3.0",
            "total_residue_count": "100",
            "is_alphafold": "False",
            "uniprot_start": "0",
            "uniprot_end": "0",
            "sequence_identity": "0.000",
            "chain_length": "0",
            "passed": "False",
            "output_file": "",
            "discard_reason": "Line 1\nLine 2",
            "discard_reason_type": "ValueError",
        },
        {
            "uniprot_accession": "P12345",
            "chain_id": "",
            "pdb_id": "d.cif.gz",
            "input_file": "/fake/d.cif.gz",
            "resolution": "4.0",
            "total_residue_count": "100",
            "is_alphafold": "False",
            "uniprot_start": "0",
            "uniprot_end": "0",
            "sequence_identity": "0.000",
            "chain_length": "0",
            "passed": "False",
            "output_file": "",
            "discard_reason": "value1, value2, value3",
            "discard_reason_type": "ValueError",
        },
        {
            "uniprot_accession": "P12345",
            "chain_id": "",
            "pdb_id": "e.cif.gz",
            "input_file": "/fake/e.cif.gz",
            "resolution": "5.0",
            "total_residue_count": "100",
            "is_alphafold": "False",
            "uniprot_start": "0",
            "uniprot_end": "0",
            "sequence_identity": "0.000",
            "chain_length": "0",
            "passed": "False",
            "output_file": "",
            "discard_reason": 'Error: "invalid" data\nat line 42',
            "discard_reason_type": "ValueError",
        },
    ]
    assert rows == expected
