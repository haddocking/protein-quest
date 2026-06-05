import gzip
from dataclasses import replace
from pathlib import Path

import pytest

from protein_quest.filters.resolution import (
    ResolutionFilterStatistics,
    copy_resolution_statistics,
    filter_files_on_resolution,
    iter_resolution_statistics,
    load_resolution_statistics,
    sort_resolution_statistics,
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
) -> ResolutionFilterStatistics:
    if pdb_id is None:
        pdb_id = filename
    return ResolutionFilterStatistics(
        id=pdb_id,
        input_file=Path(f"/fake/{filename}"),
        uniprot_accession=accession,
        resolution=resolution,
        total_residue_count=total_residue_count,
        is_alphafold=is_alphafold,
        uniprot_start=uniprot_start,
        uniprot_end=uniprot_end,
        sequence_identity=sequence_identity,
        chain_length=chain_length,
        passed=False,
        output_file=None,
        discard_reason=discard_reason,
    )


def assert_resolution_filter_statistics(
    result: ResolutionFilterStatistics,
    expected: ResolutionFilterStatistics,
) -> None:
    assert result.input_file == expected.input_file
    assert result.id == expected.id
    assert result.uniprot_accession == expected.uniprot_accession
    assert result.resolution == expected.resolution
    assert result.total_residue_count == expected.total_residue_count
    assert result.is_alphafold == expected.is_alphafold
    assert result.uniprot_start == expected.uniprot_start
    assert result.uniprot_end == expected.uniprot_end
    assert result.sequence_identity == pytest.approx(expected.sequence_identity, rel=1e-3, abs=0.0)
    assert result.chain_length == expected.chain_length
    assert result.passed == expected.passed
    assert result.output_file == expected.output_file


def test_resolution_filter_statistics_is_hashable():
    stats = _make_stats("a.cif.gz", "P12345", resolution=1.0)
    assert isinstance(hash(stats), int)


class TestIterResolutionStatistics:
    def test_metadata_in_order(self, sample_cif: Path, sample2_cif: Path, af_cif: Path, nmr_cif: Path):
        input_files = [sample_cif, sample2_cif, af_cif, nmr_cif]
        results = list(iter_resolution_statistics(input_files))

        expected = [
            ResolutionFilterStatistics(
                id="3JRSB2A",
                input_file=sample_cif,
                uniprot_accession="Q8VZS8",
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
        for result, expected_result in zip(results, expected, strict=True):
            assert_resolution_filter_statistics(result, expected_result)

    def test_empty_input(self):
        assert list(iter_resolution_statistics([])) == []


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
                id="3JRSB2A",
                input_file=sample_cif,
                uniprot_accession="Q8VZS8",
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

        for result, expected_result in zip(results, expected, strict=True):
            assert_resolution_filter_statistics(result, expected_result)

    def test_empty_input(self):
        assert load_resolution_statistics([], scheduler_address="sequential") == []

    def test_multiple_accessions_uniprotlessresult(self, multi_accession_cif: Path):
        result = load_resolution_statistics([multi_accession_cif], scheduler_address="sequential")

        expected = ResolutionFilterStatistics(
            id="1A02",
            input_file=multi_accession_cif,
            uniprot_accession=None,
            resolution=2.7,
            total_residue_count=513,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=513,
            passed=False,
            output_file=None,
        )
        assert len(result) == 1
        assert_resolution_filter_statistics(result[0], expected)


class TestFilterFilesOnResolution:
    def test_no_uniprot_does_not_pass(self, sample2_cif: Path, tmp_path: Path):
        no_uniprot = tmp_path / "no-uniprot.cif.gz"
        with gzip.open(sample2_cif, "rt", encoding="utf-8") as handle:
            text = handle.read()
        # Remove UniProt database entries while keeping a valid mmCIF structure.
        no_uniprot.write_bytes(gzip.compress(text.replace("UNP", "XXX").encode("utf-8")))

        output_dir = tmp_path / "output"
        results = list(
            filter_files_on_resolution(
                input_files=[no_uniprot], output_dir=output_dir, top=1, min_sequence_identity=0.0
            )
        )

        expected = ResolutionFilterStatistics(
            id="2Y29",
            input_file=no_uniprot,
            uniprot_accession=None,
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=8,
            passed=False,
            output_file=None,
        )
        assert len(results) == 1
        assert_resolution_filter_statistics(results[0], expected)
        assert output_dir.exists()
        assert list(output_dir.iterdir()) == []

    def test_no_uniprot_can_pass_when_grouping_disabled(self, sample2_cif: Path, tmp_path: Path):
        no_uniprot = tmp_path / "no-uniprot.cif.gz"
        with gzip.open(sample2_cif, "rt", encoding="utf-8") as handle:
            text = handle.read()
        # Remove UniProt database entries while keeping a valid mmCIF structure.
        no_uniprot.write_bytes(gzip.compress(text.replace("UNP", "XXX").encode("utf-8")))

        output_dir = tmp_path / "output"
        results = list(
            filter_files_on_resolution(
                input_files=[no_uniprot], output_dir=output_dir, top=1, group_by=False, min_sequence_identity=0.0
            )
        )

        expected = ResolutionFilterStatistics(
            id="2Y29",
            input_file=no_uniprot,
            uniprot_accession=None,
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=8,
            passed=True,
            output_file=output_dir / no_uniprot.name,
        )
        assert len(results) == 1
        assert_resolution_filter_statistics(results[0], expected)
        assert results[0].output_file is not None and results[0].output_file.exists()


@pytest.mark.parametrize("group_by, coverage", [(False, False), (True, False), (False, True), (True, True)])
def test_sort_resolution_statistics_order_alphafold_xray_nmr(group_by: bool, coverage: bool):
    shared = {
        "accession": "P12345",
        "sequence_identity": 1.0,
        "uniprot_start": 1,
        "uniprot_end": 100,
        "chain_length": 100,
    }
    af = _make_stats("af.cif.gz", resolution=0.0, is_alphafold=True, **shared)
    xray = _make_stats("xray.cif.gz", resolution=1.0, **shared)
    nmr = _make_stats("nmr.cif.gz", resolution=0.0, **shared)

    results = sort_resolution_statistics([xray, nmr, af], top=3, coverage=coverage, group_by=group_by)

    # any coverage/groupby combi should give same order
    expected = [af, xray, nmr]
    assert results == expected


class TestSortResolutionStatisticsYesGroupByNoCoverage:
    def test_groupby_accession(self):
        a1 = _make_stats("a1.cif.gz", "P11111", resolution=1.0)
        a2 = _make_stats("a2.cif.gz", "P11111", resolution=2.0)
        b1 = _make_stats("b1.cif.gz", "P22222", resolution=1.5)
        b2 = _make_stats("b2.cif.gz", "P22222", resolution=3.0)
        # Structures without uniprot always at bottom
        c1 = _make_stats("c1.cif.gz", None, resolution=0.5)
        c2 = _make_stats("c2.cif.gz", None, resolution=0.1)

        results = sort_resolution_statistics([a2, a1, b2, b1, c1, c2], top=1, coverage=False, group_by=True)

        expected = [
            replace(a1, passed=True),
            replace(a2, passed=False),
            replace(b1, passed=True),
            replace(b2, passed=False),
            replace(c2, passed=False),
            replace(c1, passed=False),
        ]
        assert results == expected

    def test_top_limit(self):
        a = _make_stats("a.cif.gz", "P12345", resolution=1.0)
        b = _make_stats("b.cif.gz", "P12345", resolution=2.0)
        c = _make_stats("c.cif.gz", "P12345", resolution=3.0)

        results = sort_resolution_statistics([a, b, c], top=2, coverage=False, group_by=True)

        passed_names = {r.input_file.name for r in results if r.passed}
        assert passed_names == {"a.cif.gz", "b.cif.gz"}

    def test_none_accession_is_skipped(self):
        good = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        no_acc = _make_stats("z.cif.gz", None, resolution=1.0)

        results = sort_resolution_statistics([good, no_acc], top=10, coverage=False, group_by=True)

        expected = [
            replace(good, passed=True),
            replace(no_acc, passed=False),
        ]
        assert results == expected

    def test_results_sorted_by_resolution_within_group(self):
        a = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        b = _make_stats("b.cif.gz", "P12345", resolution=3.0)
        c = _make_stats("c.cif.gz", "P12345", resolution=1.0)

        results = sort_resolution_statistics([a, b, c], top=10, coverage=False, group_by=True)

        assert [r.input_file.name for r in results] == ["c.cif.gz", "a.cif.gz", "b.cif.gz"]


class TestSortResolutionStatisticsNoGroupbyNoCoverage:
    def test_groupby_none_uses_global_top(self):
        a = _make_stats("a.cif.gz", "P11111", resolution=1.0)  # best
        b = _make_stats("b.cif.gz", "P22222", resolution=1.5)
        c = _make_stats("c.cif.gz", "P33333", resolution=2.0)  # worst

        results = sort_resolution_statistics([c, a, b], top=1, coverage=False, group_by=False)

        expected = [
            replace(a, passed=True),
            replace(b, passed=False),
            replace(c, passed=False),
        ]
        assert results == expected

    def test_missing_accession_ignored(self):
        with_acc = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        no_acc = _make_stats("z.cif.gz", None, resolution=1.0)  # best

        results = sort_resolution_statistics([with_acc, no_acc], top=1, coverage=False, group_by=False)

        expected = [
            replace(no_acc, passed=True),
            replace(with_acc, passed=False),
        ]
        assert results == expected

    def test_discarded_last(self):
        good = _make_stats("a.cif.gz", "P12345", resolution=1.0)
        discarded = _make_stats("b.cif.gz", "P12345", resolution=0.5, discard_reason=ValueError("reason1"))

        results = sort_resolution_statistics([discarded, good], top=1, coverage=False, group_by=False)

        expected = [
            replace(good, passed=True),
            replace(discarded, passed=False),
        ]
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
        n1 = _make_stats(
            "n1.cif.gz",
            None,
            resolution=0.5,
            chain_length=50,
        )  # without accession, between with accession and discarded
        d1 = _make_stats(
            "d1.cif.gz",
            "P33333",
            resolution=0.8,
            discard_reason=ValueError("reason1"),
        )  # baddest

        results = sort_resolution_statistics([p2b, n1, p1b, d1, p2a, p1a], top=2, coverage=True, group_by=False)

        expected = [
            replace(p1a, passed=True),
            replace(p2a, passed=True),
            replace(p1b, passed=False),
            replace(p2b, passed=False),
            replace(n1, passed=False),
            replace(d1, passed=False),
        ]

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
        n1 = _make_stats(
            "n1.cif.gz",
            None,
            resolution=0.5,
            chain_length=100,
        )  # without accession, between with accession and discarded
        d1 = _make_stats(
            "d1.cif.gz",
            "P33333",
            resolution=0.8,
            discard_reason=ValueError("reason1"),
        )  # baddest

        results = sort_resolution_statistics([p2b, n1, p1b, p1c, d1, p2a, p1a], top=2, coverage=True, group_by=True)

        expected = [
            replace(p1a, passed=True),
            replace(p1b, passed=True),
            replace(p1c, passed=False),
            replace(p2a, passed=True),
            replace(p2b, passed=True),
            replace(n1, passed=False),  # not passed because top filled by uniprot groups
            replace(d1, passed=False),
        ]

        assert results == expected

    def test_top_filled_with_accessionless(self):
        p1a = _make_stats(
            "p1a.cif.gz",
            "P11111",
            resolution=1.0,
            sequence_identity=1.00,
            chain_length=100,
            uniprot_start=1,
            uniprot_end=100,
        )
        n1 = _make_stats(
            "n1.cif.gz",
            None,
            resolution=0.5,
            chain_length=100,
        )  # without accession, best resolution
        n2 = _make_stats(
            "n2.cif.gz",
            None,
            resolution=0.6,
            chain_length=100,
        )  # without accession, second best resolution

        results = sort_resolution_statistics([n1, n2, p1a], top=2, coverage=True, group_by=True)

        expected = [
            replace(p1a, passed=True),
            replace(n1, passed=True),
            replace(n2, passed=False),
        ]

        assert results == expected


class TestCopyResolutionStatistics:
    def test_passed_is_copied_and_output_set(self, sample2_cif: Path, tmp_path: Path):
        stats = ResolutionFilterStatistics(
            id="2Y29",
            input_file=sample2_cif,
            uniprot_accession="P05067",
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

        results = list(copy_resolution_statistics([stats], output_dir))

        assert len(results) == 1
        assert results[0].output_file == output_dir / sample2_cif.name
        assert results[0].output_file is not None and results[0].output_file.exists()

    def test_not_passed_is_unchanged(self, sample2_cif: Path, tmp_path: Path):
        stats = ResolutionFilterStatistics(
            id="2Y29",
            input_file=sample2_cif,
            uniprot_accession="P05067",
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


class TestDiscardReasons:
    def test_load_invalid_file_sets_discard_reason(self, tmp_path: Path):
        invalid_file = tmp_path / "invalid.cif.gz"
        invalid_file.write_text("not a valid structure file")

        results = load_resolution_statistics([invalid_file], scheduler_address="sequential")

        assert len(results) == 1
        assert results[0].discard_reason is not None
        assert results[0].passed is False
        assert results[0].uniprot_accession is None

    def test_sequence_identity_filter_sets_discard_reason(self, sample2_cif: Path, tmp_path: Path):
        output_dir = tmp_path / "output"

        results = list(
            filter_files_on_resolution(
                input_files=[sample2_cif],
                output_dir=output_dir,
                top=10,
                min_sequence_identity=1.5,  # impossible value
            )
        )

        assert len(results) == 1
        assert results[0].discard_reason is not None
        assert results[0].passed is False
        assert "Sequence identity" in str(results[0].discard_reason)
        assert results[0].output_file is None

    def test_grouping_excludes_discarded_from_ranking(self):
        good1 = _make_stats("good1.cif.gz", "P12345", resolution=1.0)
        good2 = _make_stats("good2.cif.gz", "P12345", resolution=2.0)
        discarded = replace(_make_stats("bad.cif.gz", "P12345", resolution=0.5), discard_reason=ValueError("test"))

        results = sort_resolution_statistics([discarded, good2, good1], top=1, coverage=False, group_by=True)

        # good1 should pass (best resolution), discarded should not be considered for passing
        assert results[0].input_file.name == "good1.cif.gz"
        assert results[0].passed is True
        assert results[1].input_file.name == "good2.cif.gz"
        assert results[1].passed is False
        assert results[2].input_file.name == "bad.cif.gz"
        assert results[2].passed is False

    def test_coverage_grouping_excludes_discarded_from_ranking(self):
        good1 = _make_stats("good1.cif.gz", "P12345", resolution=1.0, uniprot_start=1, uniprot_end=10)
        good2 = _make_stats("good2.cif.gz", "P12345", resolution=2.0, uniprot_start=20, uniprot_end=30)
        discarded = replace(
            _make_stats("bad.cif.gz", "P12345", resolution=0.5, uniprot_start=1, uniprot_end=10),
            discard_reason=ValueError("test"),
        )

        results = sort_resolution_statistics([good2, good1, discarded], top=1, coverage=True, group_by=False)

        # good1 should pass (best resolution in its cluster), discarded should not be considered
        passed = [r for r in results if r.passed]
        assert len(passed) == 1
        assert passed[0].input_file.name == "good1.cif.gz"
        # discarded should still be in results
        assert any(r.input_file.name == "bad.cif.gz" for r in results)

    def test_coverage_grouping_appends_accessionless_in_grouped_mode(self):
        good1 = _make_stats("good1.cif.gz", "P12345", resolution=1.0, uniprot_start=1, uniprot_end=10)
        good2 = _make_stats("good2.cif.gz", "P12345", resolution=1.5, uniprot_start=20, uniprot_end=30)
        no_acc = _make_stats("no_acc.cif.gz", None, resolution=0.5)

        # With top=3, both good structures pass (top per group), leaving room for 1 accessionless
        results = sort_resolution_statistics([good1, good2, no_acc], top=3, coverage=True, group_by=False)

        # Both good structures should pass (limited by available structures in group)
        assert results[0].input_file.name == "good1.cif.gz"
        assert results[0].passed is True
        assert results[1].input_file.name == "good2.cif.gz"
        assert results[1].passed is True
        # no_acc should be appended and marked as passed to fill quota (2 passed + 1 accessionless = 3)
        assert results[2].input_file.name == "no_acc.cif.gz"
        assert results[2].passed is True

    def test_coverage_grouping_appends_accessionless_without_passing_when_quota_full(self):
        good1 = _make_stats("good1.cif.gz", "P12345", resolution=1.0, uniprot_start=1, uniprot_end=10)
        good2 = _make_stats("good2.cif.gz", "P12345", resolution=1.5, uniprot_start=20, uniprot_end=30)
        no_acc = _make_stats("no_acc.cif.gz", None, resolution=0.5)

        results = sort_resolution_statistics([good1, good2, no_acc], top=2, coverage=True, group_by=False)

        # both good structures should pass
        passed = [r for r in results if r.passed]
        assert len(passed) == 2
        # no_acc should still be in results but not passed
        no_acc_result = [r for r in results if r.input_file.name == "no_acc.cif.gz"]
        assert len(no_acc_result) == 1
        assert no_acc_result[0].passed is False
