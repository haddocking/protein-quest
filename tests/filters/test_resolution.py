import gzip
import logging
from dataclasses import replace
from pathlib import Path

import pytest

from protein_quest.filters.resolution import (
    ResolutionFilterStatistics,
    copy_resolution_statistics,
    filter_files_on_resolution,
    group_resolution_statistics,
    iter_resolution_statistics,
    resolution_sort_key,
)


def _make_stats(
    filename: str,
    accession: str | None,
    resolution: float,
    total_residue_count: int = 100,
    is_alphafold: bool = False,
) -> ResolutionFilterStatistics:
    return ResolutionFilterStatistics(
        input_file=Path(f"/fake/{filename}"),
        uniprot_accession=accession,
        resolution=resolution,
        total_residue_count=total_residue_count,
        is_alphafold=is_alphafold,
        passed=False,
        output_file=None,
    )


def test_resolution_sort_key():
    a = _make_stats("a.cif.gz", "P12345", resolution=1.0)
    b = _make_stats("b.cif.gz", "P12345", resolution=2.0)
    c = _make_stats("c.cif.gz", "P12345", resolution=0.0) # undesirable
    af = _make_stats("af.cif.gz", "P12345", resolution=0.0, is_alphafold=True) # best
    nu = _make_stats("nu.cif.gz", None, resolution=3.5) # uniprot is ignored

    assert sorted([b, c, a, af, nu], key=resolution_sort_key) == [af, a, b, nu, c]


class TestIterResolutionStatistics:
    def test_metadata_in_order(self, sample_cif: Path, sample2_cif: Path, af_cif: Path, nmr_cif: Path):
        input_files = [sample_cif, sample2_cif, af_cif, nmr_cif]
        results = list(iter_resolution_statistics(input_files))

        expected = [
            ResolutionFilterStatistics(
                input_file=sample_cif,
                uniprot_accession="Q8VZS8",
                resolution=2.05,
                total_residue_count=173,
                is_alphafold=False,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                input_file=sample2_cif,
                uniprot_accession="P05067",
                resolution=2.3,
                total_residue_count=8,
                is_alphafold=False,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                input_file=af_cif,
                uniprot_accession="A0A0C5B5G6",
                resolution=0.0,
                total_residue_count=16,
                is_alphafold=True,
                passed=False,
                output_file=None,
            ),
            ResolutionFilterStatistics(
                input_file=nmr_cif,
                uniprot_accession="P05067",
                resolution=0.0,
                total_residue_count=28,
                is_alphafold=False,
                passed=False,
                output_file=None,
            ),
        ]
        assert results == expected

    def test_empty_input(self):
        assert list(iter_resolution_statistics([])) == []


def test_filter_files_on_resolution_no_uniprot_does_not_pass(sample2_cif: Path, tmp_path: Path):
    no_uniprot = tmp_path / "no-uniprot.cif.gz"
    with gzip.open(sample2_cif, "rt", encoding="utf-8") as handle:
        text = handle.read()
    # Remove UniProt database entries while keeping a valid mmCIF structure.
    no_uniprot.write_bytes(gzip.compress(text.replace("UNP", "XXX").encode("utf-8")))

    output_dir = tmp_path / "output"
    results = list(filter_files_on_resolution(input_files=[no_uniprot], output_dir=output_dir, top=1))

    expected = ResolutionFilterStatistics(
        input_file=no_uniprot,
        uniprot_accession=None,
        resolution=2.3,
        total_residue_count=8,
        is_alphafold=False,
        passed=False,
        output_file=None,
    )
    assert results == [expected]
    assert output_dir.exists()
    assert list(output_dir.iterdir()) == []


def test_filter_files_on_resolution_no_uniprot_can_pass_when_grouping_disabled(sample2_cif: Path, tmp_path: Path):
    no_uniprot = tmp_path / "no-uniprot.cif.gz"
    with gzip.open(sample2_cif, "rt", encoding="utf-8") as handle:
        text = handle.read()
    # Remove UniProt database entries while keeping a valid mmCIF structure.
    no_uniprot.write_bytes(gzip.compress(text.replace("UNP", "XXX").encode("utf-8")))

    output_dir = tmp_path / "output"
    results = list(filter_files_on_resolution(input_files=[no_uniprot], output_dir=output_dir, top=1, group_by=None))

    expected = ResolutionFilterStatistics(
        input_file=no_uniprot,
        uniprot_accession=None,
        resolution=2.3,
        total_residue_count=8,
        is_alphafold=False,
        passed=True,
        output_file=output_dir / no_uniprot.name,
    )
    assert results == [expected]
    assert results[0].output_file is not None and results[0].output_file.exists()


class TestGroupResolutionStatistics:
    def test_groupby_accession(self):
        a1 = _make_stats("a1.cif.gz", "P11111", resolution=1.0)
        a2 = _make_stats("a2.cif.gz", "P11111", resolution=2.0)
        b1 = _make_stats("b1.cif.gz", "P22222", resolution=1.5)
        b2 = _make_stats("b2.cif.gz", "P22222", resolution=3.0)
        # Structures without uniprot are never passed
        c1 = _make_stats("c1.cif.gz", None, resolution=0.5)
        c2 = _make_stats("c2.cif.gz", None, resolution=0.1)

        results = group_resolution_statistics([a2, a1, b2, b1, c1, c2], top=1, group_by="uniprot_accession")

        expected = [
            replace(a1, passed=True),
            replace(a2, passed=False),
            replace(b1, passed=True),
            replace(b2, passed=False),
            replace(c1, passed=False),
            replace(c2, passed=False),
        ]
        assert results == expected

    def test_groupby_none_uses_global_top(self):
        # Same accession-based ranking key as production code; global top=1 should keep only a.cif.gz.
        best = _make_stats("a.cif.gz", "P11111", resolution=1.0)
        second = _make_stats("b.cif.gz", "P22222", resolution=1.5)
        third = _make_stats("c.cif.gz", "P33333", resolution=2.0)

        results = group_resolution_statistics([third, best, second], top=1, group_by=None)

        expected = [
            replace(best, passed=True),
            replace(second, passed=False),
            replace(third, passed=False),
        ]
        assert results == expected

    def test_top_limit(self):
        a = _make_stats("a.cif.gz", "P12345", resolution=1.0)
        b = _make_stats("b.cif.gz", "P12345", resolution=2.0)
        c = _make_stats("c.cif.gz", "P12345", resolution=3.0)

        results = group_resolution_statistics([a, b, c], top=2)

        passed_names = {r.input_file.name for r in results if r.passed}
        assert passed_names == {"a.cif.gz", "b.cif.gz"}

    def test_none_accession_is_skipped(self, caplog: pytest.LogCaptureFixture):
        good = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        no_acc = _make_stats("z.cif.gz", None, resolution=1.0)

        with caplog.at_level(logging.WARNING):
            results = group_resolution_statistics([good, no_acc], top=10)

        assert results[-1] is no_acc
        assert not results[-1].passed
        assert "z.cif.gz" in caplog.text

    def test_results_sorted_by_filename_within_group(self):
        c = _make_stats("c.cif.gz", "P12345", resolution=1.0)
        a = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        b = _make_stats("b.cif.gz", "P12345", resolution=3.0)

        results = group_resolution_statistics([c, a, b], top=10)

        assert [r.input_file.name for r in results] == ["a.cif.gz", "b.cif.gz", "c.cif.gz"]

    def test_groups_are_independent(self):
        p1 = _make_stats("p1.cif.gz", "P11111", resolution=1.0)
        p2 = _make_stats("p2.cif.gz", "P22222", resolution=2.0)

        results = group_resolution_statistics([p1, p2], top=1)

        assert all(r.passed for r in results)

    def test_none_grouping_uses_global_top(self):
        # Same accession-based ranking key as production code; global top=1 should keep only a.cif.gz.
        best = _make_stats("a.cif.gz", "P11111", resolution=1.0)
        second = _make_stats("b.cif.gz", "P22222", resolution=1.5)
        third = _make_stats("c.cif.gz", "P33333", resolution=2.0)

        results = group_resolution_statistics([third, best, second], top=1, group_by=None)

        passed_names = {r.input_file.name for r in results if r.passed}
        assert passed_names == {"a.cif.gz"}
        assert [r.input_file.name for r in results] == ["a.cif.gz", "b.cif.gz", "c.cif.gz"]

    def test_none_grouping_does_not_warn_for_missing_accession(self, caplog: pytest.LogCaptureFixture):
        with_acc = _make_stats("a.cif.gz", "P12345", resolution=2.0)
        no_acc = _make_stats("z.cif.gz", None, resolution=1.0)

        with caplog.at_level(logging.WARNING):
            results = group_resolution_statistics([with_acc, no_acc], top=1, group_by=None)

        assert all("No UniProt accession found" not in record.message for record in caplog.records)
        passed_names = {r.input_file.name for r in results if r.passed}
        assert passed_names == {"z.cif.gz"}


class TestCopyResolutionStatistics:
    def test_passed_is_copied_and_output_set(self, sample2_cif: Path, tmp_path: Path):
        stats = ResolutionFilterStatistics(
            input_file=sample2_cif,
            uniprot_accession="P05067",
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
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
            input_file=sample2_cif,
            uniprot_accession="P05067",
            resolution=2.3,
            total_residue_count=8,
            is_alphafold=False,
            passed=False,
            output_file=None,
        )
        output_dir = tmp_path / "output"

        results = list(copy_resolution_statistics([stats], output_dir))

        assert results[0].output_file is None
        assert not any(f for f in output_dir.iterdir() if f.is_file())
