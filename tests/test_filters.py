from pathlib import Path

import pytest

from protein_quest.filters import ChainFilterStatistics, filter_files_on_chain
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
