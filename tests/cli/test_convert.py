"""Convert CLI tests for protein-quest."""

import csv
from pathlib import Path

import pytest

from protein_quest.cli import main
from protein_quest.structure.formats import read_structure
from protein_quest.structure.uniprot import structure_to_uniprot


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

    assert (output_dir / sample_cif.name).exists()
    captured = capsys.readouterr()
    assert ".cif.gz" in captured.err


def test_convert_structures_with_injected_uniprot(no_uniprot_cif: Path, tmp_path: Path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / no_uniprot_cif.name).symlink_to(no_uniprot_cif)
    output_dir = tmp_path / "output"
    pdb2uniprotcsv = tmp_path / "pdb2uniprot.csv"
    with pdb2uniprotcsv.open("w", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "chain", "uniprot_accession"])
        writer.writeheader()
        writer.writerows(
            [
                {"pdb_id": "2Y29", "chain": "A", "uniprot_accession": "P01100"},
            ]
        )

    main(
        [
            "convert",
            "structures",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--output-format",
            ".cif.gz",
            "--uniprots",
            str(pdb2uniprotcsv),
        ]
    )

    output_file = output_dir / (no_uniprot_cif.name + ".gz")
    assert output_file.exists()
    structure = read_structure(output_file)
    injected_uniprot = structure_to_uniprot(structure)
    expected = {"2Y29": {("A", "P01100")}}
    assert injected_uniprot == expected


@pytest.mark.parametrize(
    ("chain_id", "extra_args"),
    [
        pytest.param("B", [], id="defaults-to-auth"),
        pytest.param("A", ["--chain-system", "label"], id="explicit-label-system"),
    ],
)
def test_convert_structures_with_injected_uniprot_chain_system(
    cif_8rw8: Path,
    tmp_path: Path,
    chain_id: str,
    extra_args: list[str],
    caplog: pytest.LogCaptureFixture,
):
    caplog.set_level("WARNING")
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    local_cif = input_dir / cif_8rw8.name
    local_cif.symlink_to(cif_8rw8)

    output_dir = tmp_path / "output"
    pdb2uniprotcsv = tmp_path / "pdb2uniprot.csv"
    pdb_id = read_structure(cif_8rw8).name
    with pdb2uniprotcsv.open("w", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "chain", "uniprot_accession"])
        writer.writeheader()
        writer.writerow({"pdb_id": pdb_id, "chain": chain_id, "uniprot_accession": "P12345"})

    main(
        [
            "convert",
            "structures",
            str(input_dir),
            "--output-dir",
            str(output_dir),
            "--output-format",
            ".cif.gz",
            "--uniprots",
            str(pdb2uniprotcsv),
            *extra_args,
        ]
    )

    output_file = output_dir / cif_8rw8.name
    assert output_file.exists()

    structure = read_structure(output_file)
    injected_uniprot = structure_to_uniprot(structure)
    assert any(uniprot == "P12345" for _, uniprot in injected_uniprot[pdb_id])
    assert "Expected: {('B', 'P12345')}" in caplog.text


def test_convert_clusters_writes_clusters_output(
    sample2_cif: Path,
    af_cif: Path,
    nmr_cif: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
):
    """convert clusters always writes the main clusters output file."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / sample2_cif.name).symlink_to(sample2_cif)
    (input_dir / af_cif.name).symlink_to(af_cif)
    (input_dir / nmr_cif.name).symlink_to(nmr_cif)

    clusters_fn = tmp_path / "clusters.csv"
    stats_fn = tmp_path / "stats.csv"

    main(["convert", "clusters", str(input_dir), str(clusters_fn)])

    assert clusters_fn.exists()
    assert not stats_fn.exists()

    with clusters_fn.open() as handle:
        cluster_rows = list(csv.DictReader(handle))

    assert len(cluster_rows) == 3
    assert {row["uniprot_accession"] for row in cluster_rows} == {"P05067", "A0A0C5B5G6"}

    captured = capsys.readouterr()
    assert "Wrote clusters for 2 uniprot accessions" in captured.err


def test_convert_clusters_optional_outputs(
    sample2_cif: Path,
    nmr_cif: Path,
    tmp_path: Path,
):
    """convert clusters writes optional stats and clustered artifacts when requested."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / sample2_cif.name).symlink_to(sample2_cif)
    (input_dir / nmr_cif.name).symlink_to(nmr_cif)

    clusters_fn = tmp_path / "clusters.csv"
    stats_fn = tmp_path / "stats.csv"
    condensed_fn = tmp_path / "distances.csv"
    linkage_fn = tmp_path / "linkage.csv"
    dendrogram_dir = tmp_path / "dendrograms"

    main(
        [
            "convert",
            "clusters",
            str(input_dir),
            str(clusters_fn),
            "--stats",
            str(stats_fn),
            "--condensed-distances",
            str(condensed_fn),
            "--linkage-matrix",
            str(linkage_fn),
            "--dendrogram",
            str(dendrogram_dir),
        ]
    )

    dendrogram_fn = dendrogram_dir / "P05067_dendrogram.nwk"

    assert clusters_fn.exists()
    assert stats_fn.exists()
    assert condensed_fn.exists()
    assert linkage_fn.exists()
    assert dendrogram_fn.exists()

    with condensed_fn.open() as handle:
        condensed_rows = list(csv.DictReader(handle))
    assert list(condensed_rows[0].keys()) == ["uniprot_accession", "structure_i", "structure_j", "distance"]
    assert len(condensed_rows) == 1
    assert condensed_rows[0]["uniprot_accession"] == "P05067"

    with linkage_fn.open() as handle:
        linkage_rows = list(csv.DictReader(handle))
    assert next(iter(linkage_rows[0].keys())) == "uniprot_accession"
    assert len(linkage_rows) == 1
    assert linkage_rows[0]["uniprot_accession"] == "P05067"

    dendrogram_text = dendrogram_fn.read_text(encoding="utf-8")
    assert dendrogram_text.endswith(";")


@pytest.mark.parametrize(
    "option",
    [
        "--condensed-distances",
        "--linkage-matrix",
        "--dendrogram",
    ],
)
def test_convert_clusters_optional_output_requires_value(
    sample2_cif: Path,
    tmp_path: Path,
    option: str,
):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / sample2_cif.name).symlink_to(sample2_cif)

    output_file = tmp_path / "clusters.csv"

    with pytest.raises(SystemExit) as exc_info:
        main(["convert", "clusters", str(input_dir), str(output_file), option])

    assert exc_info.value.code == 1


def test_convert_uniprot_ungrouped_multiple_accessions(
    sample2_cif: Path,
    multi_accession_cif: Path,
    tmp_path: Path,
):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / sample2_cif.name).symlink_to(sample2_cif)
    (input_dir / multi_accession_cif.name).symlink_to(multi_accession_cif)
    output_file = tmp_path / "uniprot.txt"

    main(["convert", "uniprot", str(input_dir), str(output_file)])

    lines = output_file.read_text(encoding="utf-8").splitlines()
    assert lines == ["P01100", "P05067", "P05412", "Q13469"]


def test_convert_uniprot_grouped_multiple_accessions(
    sample2_cif: Path,
    multi_accession_cif: Path,
    tmp_path: Path,
):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    local_sample2 = input_dir / sample2_cif.name
    local_sample2.symlink_to(sample2_cif)
    local_multi = input_dir / multi_accession_cif.name
    local_multi.symlink_to(multi_accession_cif)
    output_file = tmp_path / "uniprot_grouped.txt"

    main(["convert", "uniprot", str(input_dir), str(output_file), "--grouped"])

    lines = output_file.read_text(encoding="utf-8").splitlines()
    assert lines == [
        f"{local_multi},P01100",
        f"{local_multi},P05412",
        f"{local_multi},Q13469",
        f"{local_sample2},P05067",
    ]
