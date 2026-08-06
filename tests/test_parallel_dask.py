import pytest

try:
    from protein_quest.parallel_dask import Client, MyProgressBar, configure_dask_scheduler, dask_map_with_progress
except (AttributeError, ImportError):
    pytest.skip("distributed not installed", allow_module_level=True)

from protein_quest.parallel import map_with_progress
from protein_quest.pdbe.ws import Scores


class TestMapWithProgress:
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

    def test_dask_with_address(self):
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

    def test_dask_with_big_kwarg(self):
        score_count = 1190
        item_count = 2382
        scores = {
            f"id{i:05d}": Scores(
                geometry_quality=1.0,
                data_quality=1.0,
                overall_quality=1.0,
                experiment_data_available=True,
            )
            for i in range(score_count)
        }
        pdb_ids = [f"id{i % score_count:05d}" for i in range(item_count)]

        def myfunc(pdb_id: str, scores: dict[str, Scores]) -> float | None:
            score = scores.get(pdb_id)
            if score is None:
                return 0.0
            return score.geometry_quality

        result = map_with_progress(scheduler_address=None, func=myfunc, iterable=pdb_ids, scores=scores)

        assert len(result) == item_count
        assert set(result) == {1.0}


def test_MyProgressBar_interval_env(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("TQDM_MININTERVAL", "1234")
    assert Client is not None

    with Client():
        progress_bar = MyProgressBar([])
        assert progress_bar.interval == 1234


def run_dask_map_with_progress():
    def square(x: int) -> int:
        return x**2

    assert Client is not None

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
