from pathlib import Path

import pytest

from protein_quest.filters.chain import ChainFilterStatistics, filter_file_on_chain, filter_files_on_chain
from protein_quest.structure.chains import ChainIdSystem
from protein_quest.structure.errors import ChainNotFoundError


class TestFilterFilesOnChain:
    @pytest.mark.parametrize(
        "scheduler_address,expected_progress_bar",
        [
            (None, "Completed"),  # creates a local cluster
            ("sequential", "file/s"),
        ],
    )
    def test_local_cluster(
        self,
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

    @pytest.mark.parametrize(
        ("chain_id", "chain_system"),
        [
            pytest.param("ZZ", "auth", id="invalid-auth-chain"),
            pytest.param("ZZ", "label", id="invalid-label-chain"),
            pytest.param("A", "auth", id="wrong-system-A-as-auth"),
            pytest.param("B", "label", id="wrong-system-B-as-label"),
        ],
    )
    def test_invalid_chain(
        self,
        cif_8rw8: Path,
        tmp_path: Path,
        chain_id: str,
        chain_system: ChainIdSystem,
    ):
        results = filter_files_on_chain(
            [(cif_8rw8, chain_id)],
            tmp_path,
            chain_system=chain_system,
            scheduler_address="sequential",
        )

        assert len(results) == 1
        result = results[0]
        assert not result.passed
        assert result.output_file is None
        assert isinstance(result.discard_reason, ChainNotFoundError)


@pytest.mark.parametrize(
    ("chain_id", "chain_system"),
    [
        ("B", "auth"),
        ("A", "label"),
    ],
)
def test_filter_file_on_chain_system_handling(
    cif_8rw8: Path,
    tmp_path: Path,
    chain_id: str,
    chain_system: ChainIdSystem,
):
    result = filter_file_on_chain((cif_8rw8, chain_id), tmp_path, chain_system=chain_system)

    assert result.passed is True
    assert result.output_file is not None
    assert result.output_file.exists()
    assert result.output_file.name.endswith("_B2A.cif.gz")


def test_filter_file_on_chain_force_overwrites_existing_output(sample2_cif: Path, tmp_path: Path):
    output_file = tmp_path / "2Y29_A2A.cif.gz"
    output_file.write_text("fake content")

    result = filter_file_on_chain((sample2_cif, "A"), tmp_path, force=True)

    assert result.passed is True
    assert result.output_file == output_file
    assert result.output_file.exists()
    assert result.output_file.read_bytes() != b"fake content"
