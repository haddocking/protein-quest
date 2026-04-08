"""Common CLI options for protein-quest."""

import logging
from collections.abc import Iterable
from dataclasses import dataclass
from os import linesep
from pathlib import Path
from typing import Annotated, Any, Literal

from cyclopts import Group, Parameter, validators
from cyclopts.types import NonNegativeFloat, PositiveInt, StdioPath
from rich.console import Console
from rocrate_action_recorder.adapters.cyclopts import INPUT_DIR, INPUT_FILE, OUTPUT_DIR, OUTPUT_FILE, RECORD_TRIGGER

from protein_quest.utils import Cacher, DirectoryCacher, PassthroughCacher, user_cache_root_dir

# Custom annotated types for common CLI parameters
Limit = PositiveInt
"""Type for limit parameters (positive integers > 0)."""

Timeout = PositiveInt
"""Type for timeout parameters (positive integers > 0)."""

MinResidues = PositiveInt
"""Type for minimum residue parameters (positive integers > 0)."""

MaxResidues = PositiveInt
"""Type for maximum residue parameters (positive integers > 0)."""

BatchSize = PositiveInt
"""Type for batch size parameters (positive integers > 0)."""

MinSequenceLength = PositiveInt
"""Type for minimum sequence length parameters (positive integers > 0)."""

MaxSequenceLength = PositiveInt
"""Type for maximum sequence length parameters (positive integers > 0)."""

ConfidenceThreshold = NonNegativeFloat
"""Type for confidence threshold parameters (non-negative floats >= 0)."""

console = Console(stderr=True)

cache_group = Group("Cache")
common_group = Group("Common")


def setup_logging(verbose: int, quiet: int) -> None:
    """Configure root logging and console quiet mode

    Args:
        verbose: Verbosity level
        quiet: Quietness level.

    Raises:
        ValueError: If both verbose and quiet are used together.
    """
    if verbose > 0 and quiet > 0:
        msg = "-v/--verbose and -q/--quiet cannot be used together"
        raise ValueError(msg)

    base_level = logging.WARNING
    adjustment = (quiet - verbose) * 10
    level = max(logging.DEBUG, min(logging.CRITICAL, base_level + adjustment))

    logging.root.setLevel(level)
    for handler in logging.root.handlers:
        handler.setLevel(level)
    if not logging.root.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(level)
        logging.root.addHandler(handler)
    if quiet > 1:
        # -qq silences console completely
        console.quiet = True


@Parameter(name="*", group=cache_group)
@dataclass
class CacheParameter:
    """Cache-related CLI options.

    Args:
        no_cache: Disable caching of files to central location.
        cache_dir: Directory to use as cache for files.
        copy_method: How to make target file be same file as source file.
            By default uses hardlinks to save disk space.
            Note that hardlinks only work within the same filesystem and are harder to track.
            If you want to track cached files easily then use 'symlink'.
            On Windows you need developer mode or admin privileges to create symlinks.
    """

    no_cache: Annotated[
        bool,
        Parameter(
            negative="",
        ),
    ] = False
    cache_dir: Path = user_cache_root_dir()  # noqa: RUF009
    copy_method: Literal["copy", "symlink", "hardlink"] = "hardlink"


def to_cacher(cache_params: CacheParameter | None) -> Cacher:
    """Initialize cacher based on parameters.

    Args:
        cache_params: Cache parameters from CLI.

    Returns:
        Initialized cacher instance.
    """
    if cache_params is None or cache_params.no_cache:
        return PassthroughCacher()
    return DirectoryCacher(cache_dir=cache_params.cache_dir, copy_method=cache_params.copy_method)


@Parameter(name="*", group=common_group)
@dataclass
class Common:
    """Common CLI options shared across all commands.

    Args:
        verbose: Increase verbosity (use multiple times for more detail).
        quiet: Decrease verbosity (use multiple times for less output).
        prov: Whether to write provenance information about the
            command execution to ro-crate-metadata.json file.
    """

    verbose: Annotated[
        int,
        Parameter(
            name=("-v", "--verbose"),
            count=True,
        ),
    ] = 0

    quiet: Annotated[
        int,
        Parameter(
            name=("-q", "--quiet"),
            count=True,
        ),
    ] = 0

    prov: Annotated[
        bool,
        Parameter(
            negative="",
        ),
        RECORD_TRIGGER,
    ] = False

    def __post_init__(self):
        """Automatically configure logging when Common instance is created."""
        setup_logging(verbose=self.verbose, quiet=self.quiet)


def write_lines(file: StdioPath, lines: Iterable[str]):
    """Write lines to a file

    If file is "-", writes to stdout.

    Creates parent directories if they do not exist.

    Args:
        file: Path to output file.
        lines: Iterable of lines to write.
    """
    if str(file) != "-":
        file.parent.mkdir(parents=True, exist_ok=True)
    with file.open("w", encoding="utf-8") as f:
        for line in lines:
            f.write(line + linesep)


class StdioPathValidator(validators.Path):
    """Custom Path validator that allows "-" for stdin/stdout.

    This validator checks if the path is "-", and if so, it bypasses Path validation.

    If given '-', it also checks if a file named '-' already exists, which would conflict with stdin/stdout usage.
    """

    def __call__(self, type_: Any, path: Any):
        if str(path) == "-":
            if Path("-").exists():
                msg = '"-" is reserved for stdin/stdout but a file named "-" already exists.'
                raise ValueError(msg)
            return None
        return super().__call__(type_, path)


InputFile = Annotated[StdioPath, Parameter(validator=StdioPathValidator(exists=True, dir_okay=False)), INPUT_FILE]
"""Type for input file parameters (file paths that can also be "-" for stdin)."""
InputDir = Annotated[Path, Parameter(validator=validators.Path(exists=True, file_okay=False)), INPUT_DIR]
"""Type for input directory parameters (directory paths)."""
OutputFile = Annotated[StdioPath, Parameter(validator=StdioPathValidator(dir_okay=False)), OUTPUT_FILE]
"""Type for output file parameters (file paths that can also be "-" for stdout)."""
OutputDir = Annotated[Path, Parameter(validator=validators.Path(file_okay=False)), OUTPUT_DIR]
"""Type for output directory parameters (directory paths)."""
