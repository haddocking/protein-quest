import pytest

from protein_quest.parallel import map_with_progress


class TestMapWithProgressUsingSequential:
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
