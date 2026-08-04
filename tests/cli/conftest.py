from pathlib import Path

import pytest


@pytest.fixture
def config_file_path(tmp_path: Path) -> Path:
    return tmp_path / "config.toml"


@pytest.fixture(autouse=True)
def never_use_user_config(config_file_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Make cli tests ignore any user config file
    monkeypatch.setenv("PROTEIN_QUEST_CONFIG", str(config_file_path))
