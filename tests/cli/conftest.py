import os
from pathlib import Path


def pytest_configure():
    # Make cli tests ignore any user or site config file
    cfg = Path("/dummy-config.toml")
    os.environ["PROTEIN_QUEST_CONFIG"] = str(cfg)
