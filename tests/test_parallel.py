import pytest
from distributed import Client

from protein_quest.parallel import MyProgressBar, configure_dask_scheduler, dask_map_with_progress, map_with_progress


class TestMapWithProgress:
    def test_tqdm(self, capsys: pytest.CaptureFixture[str]):
        result = map_with_progress(
            scheduler_address="sequential",
            func=lambda x: x * 2,
            iterable=[0, 1, 2, 3, 4],
            map_with_progress_options={"tqdm_desc": "tqdm-desc-test", "tqdm_unit": "tqdm-unit-item"},
        )

        assert result == [0, 2, 4, 6, 8]
        captured = capsys.readouterr()
        assert "tqdm-desc-test" in captured.err and "tqdm-unit-item" in captured.err

    def test_dask_local(self, capsys: pytest.CaptureFixture[str]):
        result = map_with_progress(
            scheduler_address=None,
            func=lambda x: x * 2,
            iterable=[0, 1, 2, 3, 4],
            map_with_progress_options={"tqdm_desc": "test parallel", "tqdm_unit": "item"},
        )

        assert result == [0, 2, 4, 6, 8]
        captured = capsys.readouterr()
        assert "Completed" in captured.err

    def test_dask_with_adress(self):
        with configure_dask_scheduler(None, name="running-cluster") as cluster:
            scheduler_address = cluster if isinstance(cluster, str) else cluster.scheduler_address

            result = map_with_progress(
                scheduler_address=scheduler_address,
                func=lambda x: x * 2,
                iterable=[0, 1, 2, 3, 4],
                map_with_progress_options={"tqdm_desc": "test parallel", "tqdm_unit": "item"},
            )

            assert result == [0, 2, 4, 6, 8]

    def test_dask_with_name(self):
        with configure_dask_scheduler(None, name="running-cluster") as cluster:
            name = cluster if isinstance(cluster, str) else cluster.name

            result = map_with_progress(
                scheduler_address=None,
                func=lambda x: x * 2,
                iterable=[0, 1, 2, 3, 4],
                map_with_progress_options={
                    "tqdm_desc": "test parallel",
                    "tqdm_unit": "item",
                    "dask_scheduler_name": name,
                },
            )

            assert result == [0, 2, 4, 6, 8]

    def test_dask_with_cluster(self):
        with configure_dask_scheduler(None, name="running-cluster") as cluster:
            if isinstance(cluster, str):
                msg = "Expected a Cluster instance, got a string."
                raise TypeError(msg)

            result = map_with_progress(
                scheduler_address=cluster,
                func=lambda x: x * 2,
                iterable=[0, 1, 2, 3, 4],
            )

            assert result == [0, 2, 4, 6, 8]

    def test_tqdm_disabled(self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]):
        monkeypatch.setenv("TQDM_DISABLE", "1")

        result = map_with_progress(
            scheduler_address="sequential",
            func=lambda x: x * 2,
            iterable=[0, 1, 2, 3, 4],
        )

        assert result == [0, 2, 4, 6, 8]
        captured = capsys.readouterr()
        assert "Completed" not in captured.err

    def test_posargs(self):
        result = map_with_progress(
            "sequential",
            lambda x, offset: x + offset,
            [0, 1, 2, 3, 4],
            None,
            10,
        )

        assert result == [10, 11, 12, 13, 14]

    def test_kwargs(self):
        result = map_with_progress(
            scheduler_address="sequential",
            func=lambda x, offset=0: x + offset,
            iterable=[0, 1, 2, 3, 4],
            offset=10,
        )

        assert result == [10, 11, 12, 13, 14]


def test_MyProgressBar_interval_env(monkeypatch: pytest.MonkeyPatch):
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


def test_dask_map_with_progress(capsys: pytest.CaptureFixture[str], caplog: pytest.LogCaptureFixture):
    caplog.set_level("INFO")

    run_dask_map_with_progress()

    captured = capsys.readouterr()
    assert "Completed" in captured.err

    assert "Follow progress on dask dashboard at" in caplog.text


def test_dask_map_with_progress_disabled(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]):
    monkeypatch.setenv("TQDM_DISABLE", "1")

    run_dask_map_with_progress()

    captured = capsys.readouterr()
    assert "Completed" not in captured.err
