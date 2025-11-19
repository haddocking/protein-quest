import pytest
from distributed import Client

from protein_quest.parallel import MyProgressBar, dask_map_with_progress


def test_MyProgressBar_interval_env(monkeypatch):
    monkeypatch.setenv("TQDM_MININTERVAL", "1234")

    with Client():
        progress_bar = MyProgressBar([])
        assert progress_bar.interval == 1234


def run_dask_map_with_progress():
    def square(x: int) -> int:
        return x**2

    with Client() as client:
        result = dask_map_with_progress(
            client,
            square,
            range(5),
        )
    assert result == [0, 1, 4, 9, 16]


def test_dask_map_with_progress(capsys: pytest.CaptureFixture, caplog: pytest.LogCaptureFixture):
    caplog.set_level("INFO")

    run_dask_map_with_progress()

    captured = capsys.readouterr()
    assert "Completed" in captured.err

    assert "Follow progress on dask dashboard at" in caplog.text


def test_dask_map_with_progress_disabled(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
    monkeypatch.setenv("TQDM_DISABLE", "1")

    run_dask_map_with_progress()

    captured = capsys.readouterr()
    assert "Completed" not in captured.err
