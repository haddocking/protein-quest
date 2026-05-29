"""Convert CLI tests for protein-quest."""

import csv
from pathlib import Path

import pytest

from protein_quest.cli import main


def test_convert_structures_to_cifgz(sample_cif: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    input_file = input_dir / sample_cif.name
    input_file.write_bytes(sample_cif.read_bytes())
    output_dir = tmp_path / "output"

    main(
        [
            "convert",
            "structures",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--output-format",
            ".cif.gz",
        ]
    )

    assert (output_dir / "3JRS_B2A.cif.gz").exists()
    captured = capsys.readouterr()
    assert ".cif.gz" in captured.err


def test_convert_clusters_writes_default_outputs(
    sample2_cif: Path,
    af_cif: Path,
    nmr_cif: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    """convert clusters always writes clusters.csv and stats.csv."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / sample2_cif.name).symlink_to(sample2_cif)
    (input_dir / af_cif.name).symlink_to(af_cif)
    (input_dir / nmr_cif.name).symlink_to(nmr_cif)

    output_dir = tmp_path / "output"

    main(["convert", "clusters", str(input_dir), str(output_dir)])

    clusters_fn = output_dir / "clusters.csv"
    stats_fn = output_dir / "stats.csv"
    assert clusters_fn.exists()
    assert stats_fn.exists()

    with clusters_fn.open() as handle:
        cluster_rows = list(csv.DictReader(handle))
    with stats_fn.open() as handle:
        stats_rows = list(csv.DictReader(handle))

    assert len(cluster_rows) == 3
    assert {row["uniprot_accession"] for row in cluster_rows} == {"P05067", "A0A0C5B5G6"}
    assert len(stats_rows) == 2

    captured = capsys.readouterr()
    assert "Wrote clusters for 2 accessions" in captured.err


def test_convert_clusters_optional_outputs(
    sample2_cif: Path,
    nmr_cif: Path,
    tmp_path: Path,
):
    """convert clusters writes optional per-accession artifacts when requested."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / sample2_cif.name).symlink_to(sample2_cif)
    (input_dir / nmr_cif.name).symlink_to(nmr_cif)

    output_dir = tmp_path / "output"

    main(
        [
            "convert",
            "clusters",
            str(input_dir),
            str(output_dir),
            "--condensed-distances",
            "--linkage-matrix",
            "--dendrogram",
        ]
    )

    condensed_fn = output_dir / "P05067_condensed.csv"
    linkage_fn = output_dir / "P05067_linkage.csv"
    dendrogram_fn = output_dir / "P05067_dendrogram.nwk"

    assert condensed_fn.exists()
    assert linkage_fn.exists()
    assert dendrogram_fn.exists()

    with condensed_fn.open() as handle:
        condensed_rows = list(csv.DictReader(handle))
    assert len(condensed_rows) == 1

    with linkage_fn.open() as handle:
        linkage_rows = list(csv.DictReader(handle))
    assert len(linkage_rows) == 1

    dendrogram_text = dendrogram_fn.read_text(encoding="utf-8")
    assert dendrogram_text.endswith(";")
