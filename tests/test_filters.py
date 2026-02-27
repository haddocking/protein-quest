import tarfile
from pathlib import Path

import pytest

from protein_quest.filters import (
    ChainFilterStatistics,
    ResidueFilterStatistics,
    filter_chain_rows_mixed_io,
    filter_files_on_chain,
    filter_files_on_residues,
)
from protein_quest.structure import ChainNotFoundError


@pytest.mark.parametrize(
    "scheduler_address,expected_progress_bar",
    [
        (None, "Completed"),  # creates a local cluster
        ("sequential", "file/s"),
    ],
)
def test_filter_files_on_chain_local_cluster(
    sample2_cif: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    scheduler_address: str | None,
    expected_progress_bar: str,
):
    file2chains = [
        (sample2_cif, "A"),  # should pass
        (sample2_cif, "B"),  # should be discarded
    ]

    results = filter_files_on_chain(file2chains, tmp_path, scheduler_address=scheduler_address)

    expected_passed = ChainFilterStatistics(
        input_file=sample2_cif,
        chain_id="A",
        passed=True,
        output_file=tmp_path / "2Y29_A2A.cif.gz",
    )
    assert expected_passed.output_file and expected_passed.output_file.exists()
    expected_discarded = ChainFilterStatistics(
        input_file=sample2_cif,
        chain_id="B",
        passed=False,
        output_file=None,
        discard_reason=ChainNotFoundError("B", sample2_cif, {"A"}),
    )
    assert results == [expected_passed, expected_discarded]

    _, stderr = capsys.readouterr()
    assert expected_progress_bar in stderr


def test_filter_files_on_residues(sample_cif: Path, sample2_cif: Path, tmp_path: Path):
    results = list(
        filter_files_on_residues(
            input_files=[sample_cif, sample2_cif],
            output_dir=tmp_path,
            min_residues=100,
            max_residues=200,
        )
    )
    expected_passed = ResidueFilterStatistics(
        input_file=sample_cif,
        residue_count=173,
        passed=True,
        output_file=tmp_path / sample_cif.name,
    )
    assert expected_passed.output_file and expected_passed.output_file.exists()
    expected_discarded = ResidueFilterStatistics(
        input_file=sample2_cif,
        residue_count=8,
        passed=False,
        output_file=None,
    )

    assert results == [expected_passed, expected_discarded]


def test_filter_chain_rows_mixed_io_tar_to_dir(sample2_cif: Path, tmp_path: Path):
    rows = [{"pdb_id": "2Y29", "chain": "A"}]
    input_tar = tmp_path / "downloads-pdbe.tar"
    with tarfile.open(input_tar, "w") as tar:
        tar.add(sample2_cif, arcname=sample2_cif.name)

    output_dir = tmp_path / "filtered"
    nr_written, errors, discards = filter_chain_rows_mixed_io(
        rows,
        input_tar,
        output_dir,
        scheduler_address="sequential",
    )

    assert nr_written == 1
    assert errors == []
    assert discards == []
    assert (output_dir / "2Y29_A2A.cif.gz").exists()


def test_filter_chain_rows_mixed_io_dir_to_tar(sample2_cif: Path, tmp_path: Path):
    rows = [{"pdb_id": "2Y29", "chain": "A"}]
    output_tar = tmp_path / "filtered-chains.tar"

    nr_written, errors, discards = filter_chain_rows_mixed_io(
        rows,
        sample2_cif.parent,
        output_tar,
        scheduler_address="sequential",
    )

    assert nr_written == 1
    assert errors == []
    assert discards == []
    assert output_tar.exists()
    with tarfile.open(output_tar, "r") as tar:
        assert "2Y29_A2A.cif.gz" in tar.getnames()


def test_filter_chain_rows_mixed_io_tar_to_tar(sample2_cif: Path, tmp_path: Path):
    rows = [{"pdb_id": "2Y29", "chain": "A"}]
    input_tar = tmp_path / "downloads-pdbe.tar"
    with tarfile.open(input_tar, "w") as tar:
        tar.add(sample2_cif, arcname=sample2_cif.name)
    output_tar = tmp_path / "filtered-chains.tar"

    nr_written, errors, discards = filter_chain_rows_mixed_io(
        rows,
        input_tar,
        output_tar,
        scheduler_address="sequential",
    )

    assert nr_written == 1
    assert errors == []
    assert discards == []
    assert output_tar.exists()
    with tarfile.open(output_tar, "r") as tar:
        assert "2Y29_A2A.cif.gz" in tar.getnames()


def test_filter_chain_rows_mixed_io_input_tar_notfound(tmp_path: Path):
    rows = [{"pdb_id": "2Y29", "chain": "A"}]
    input_tar = tmp_path / "downloads-pdbe.tar"
    with tarfile.open(input_tar, "w"):
        pass
    output_tar = tmp_path / "filtered-chains.tar"

    nr_written, errors, discards = filter_chain_rows_mixed_io(
        rows,
        input_tar,
        output_tar,
        scheduler_address="sequential",
    )

    assert nr_written == 0
    assert len(errors) == 1
    assert "No structure file found for 2Y29" in str(errors[0])
    assert discards == []
    assert not output_tar.exists()


def test_filter_chain_rows_mixed_io_input_tar_partial_notfound(sample2_cif: Path, tmp_path: Path):
    rows = [{"pdb_id": "2Y29", "chain": "A"}, {"pdb_id": "1ABC", "chain": "A"}]
    input_tar = tmp_path / "downloads-pdbe.tar"
    with tarfile.open(input_tar, "w") as tar:
        tar.add(sample2_cif, arcname=sample2_cif.name)
    output_tar = tmp_path / "filtered-chains.tar"

    nr_written, errors, discards = filter_chain_rows_mixed_io(
        rows,
        input_tar,
        output_tar,
        scheduler_address="sequential",
    )

    assert nr_written == 1
    assert len(errors) == 1
    assert "No structure file found for 1ABC" in str(errors[0])
    assert discards == []
    assert output_tar.exists()
    with tarfile.open(output_tar, "r") as tar:
        assert "2Y29_A2A.cif.gz" in tar.getnames()


def test_filter_chain_rows_mixed_io_tar_to_tar_with_dask_scheduler(sample2_cif: Path, tmp_path: Path):
    rows = [{"pdb_id": "2Y29", "chain": "A"}, {"pdb_id": "2Y29", "chain": "A"}]
    input_tar = tmp_path / "downloads-pdbe.tar"
    with tarfile.open(input_tar, "w") as tar:
        tar.add(sample2_cif, arcname=sample2_cif.name)
    output_tar = tmp_path / "filtered-chains.tar"

    nr_written, errors, discards = filter_chain_rows_mixed_io(rows, input_tar, output_tar, scheduler_address=None)

    assert nr_written == 1
    assert errors == []
    assert discards == []
    assert output_tar.exists()
    with tarfile.open(output_tar, "r") as tar:
        assert "2Y29_A2A.cif.gz" in tar.getnames()
    assert not list(output_tar.parent.glob(".filtered-chains.worker-*.tar"))
