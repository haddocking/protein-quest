"""Module for cli parsers and handlers."""

import logging
import os
import sys
from collections.abc import Sequence
from pathlib import Path

import cyclopts
from cyclopts import App
from platformdirs import PlatformDirs
from rich.console import Console
from rich.traceback import install as install_rich_traceback
from rocrate_action_recorder.adapters.cyclopts import (
    record_cyclopts,
)

from protein_quest.__version__ import __version__
from protein_quest.cli.convert import convert_app
from protein_quest.cli.filter import filter_app
from protein_quest.cli.mcp import mcp_app
from protein_quest.cli.retrieve import retrieve_app
from protein_quest.cli.search import search_app

console = Console(stderr=True)
rprint = console.print
logger = logging.getLogger(__name__)


def make_app() -> App:
    """Make CLI App instance.

    The optional user configuration file is located at `~/.config/protein-quest/config.toml`.
    Or wherever the `PROTEIN_QUEST_CONFIG` environment variable points to.
    """
    pd = PlatformDirs("protein-quest")
    config_path = pd.user_config_path / "config.toml"
    env_path = os.environ.get("PROTEIN_QUEST_CONFIG")
    if env_path:
        config_path = Path(env_path)
    config = cyclopts.config.Toml(config_path)

    app = App(
        name="protein-quest",
        version=__version__,
        help="Protein Quest CLI",
        config=config,
    )

    app.register_install_completion_command()
    install_rich_traceback(console=console, suppress=[cyclopts])

    app.command(search_app)
    app.command(retrieve_app)
    app.command(filter_app)
    app.command(convert_app)
    app.command(mcp_app)
    return app


app4docs = make_app()
"""Only used for generating documentation, do not use for any other purpose.

As cyclopts mkdocs plugin does not support an App generator function, only an App instance.
"""


def main(argv: Sequence[str] | None = None):
    """Main entry point for the CLI.

    Args:
        argv: List of command line arguments. If None, uses sys.argv.
    """
    actual_argv = argv or sys.argv[1:]

    app = make_app()
    with record_cyclopts(app, tokens=actual_argv, dataset_license="CC BY 4.0"):
        app(actual_argv)
