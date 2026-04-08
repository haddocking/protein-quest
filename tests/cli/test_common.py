import logging
from pathlib import Path

import pytest

from protein_quest.cli.common import CacheParameter, StdioPathValidator, console, setup_logging, to_cacher
from protein_quest.utils import Cacher, DirectoryCacher, PassthroughCacher


class TestStdioPathValidator:
    def test_valid_path(self, tmp_path: Path):
        valid_file = tmp_path / "valid.txt"
        valid_file.touch()
        validator = StdioPathValidator(exists=True)

        assert validator(Path, Path(valid_file)) is None

    def test_stdin_stdout(self):
        validator = StdioPathValidator(exists=True)
        assert validator(Path, Path("-")) is None

    def test_nonexistent_path(self, tmp_path: Path):
        nonexistent_file = tmp_path / "nonexistent.txt"
        validator = StdioPathValidator(exists=True)

        with pytest.raises(ValueError, match="does not exist"):
            validator(Path, Path(nonexistent_file))

    def test_stdin_stdout_conflict(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.chdir(tmp_path)
        conflict_file = tmp_path / "-"
        conflict_file.touch()
        validator = StdioPathValidator(exists=True)

        with pytest.raises(ValueError, match='"-" is reserved for stdin/stdout'):
            validator(Path, Path("-"))


class TestSetupLogging:
    @pytest.fixture(autouse=True)
    def isolate_logging_state(self):
        old_level = logging.root.level
        old_handlers = list(logging.root.handlers)
        old_console_quiet = console.quiet
        yield
        logging.root.handlers = old_handlers
        logging.root.setLevel(old_level)
        console.quiet = old_console_quiet

    @pytest.mark.parametrize(
        "verbose, quiet, expected_level, expected_console_quiet",
        [
            (0, 0, logging.WARNING, False),
            (1, 0, logging.INFO, False),
            (2, 0, logging.DEBUG, False),
            (3, 0, logging.DEBUG, False),
            (0, 1, logging.ERROR, False),
            (0, 2, logging.CRITICAL, True),
            (0, 3, logging.CRITICAL, True),
        ],
    )
    def test_levels_and_console_quiet(
        self,
        verbose: int,
        quiet: int,
        expected_level: int,
        expected_console_quiet: bool,
    ):
        logging.root.handlers = []
        console.quiet = False

        setup_logging(verbose=verbose, quiet=quiet)

        assert logging.root.level == expected_level
        assert console.quiet == expected_console_quiet

    def test_raises_when_verbose_and_quiet_are_both_set(self):
        with pytest.raises(ValueError, match="cannot be used together"):
            setup_logging(verbose=1, quiet=1)

    @pytest.mark.parametrize(
        "verbose, quiet, expected_level",
        [
            (1, 0, logging.INFO),
            (0, 1, logging.ERROR),
        ],
    )
    def test_updates_all_existing_handlers(self, verbose: int, quiet: int, expected_level: int):
        handler1 = logging.StreamHandler()
        handler2 = logging.StreamHandler()
        handler1.setLevel(logging.CRITICAL)
        handler2.setLevel(logging.DEBUG)
        logging.root.handlers = [handler1, handler2]

        setup_logging(verbose=verbose, quiet=quiet)

        assert logging.root.level == expected_level
        assert len(logging.root.handlers) == 2
        assert all(handler.level == expected_level for handler in logging.root.handlers)

    def test_creates_stream_handler_when_no_handlers_exist(self):
        logging.root.handlers = []

        setup_logging(verbose=1, quiet=0)

        assert len(logging.root.handlers) == 1
        assert isinstance(logging.root.handlers[0], logging.StreamHandler)
        assert logging.root.handlers[0].level == logging.INFO


@pytest.mark.parametrize(
    "param, expected",
    [
        (None, PassthroughCacher()),
        (CacheParameter(no_cache=True), PassthroughCacher()),
        (CacheParameter(), DirectoryCacher()),
        (
            CacheParameter(cache_dir=Path("/tmp/cache"), copy_method="symlink"),
            DirectoryCacher(cache_dir=Path("/tmp/cache"), copy_method="symlink"),
        ),
    ],
)
def test_to_cacher(param: CacheParameter | None, expected: Cacher):
    cacher = to_cacher(param)
    assert cacher == expected
