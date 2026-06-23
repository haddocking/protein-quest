"""Errors used by multiple modules in protein_quest."""

from pathlib import Path


class ResolutionUnsetError(ValueError):
    """Indicates that a structure has no resolution set."""

    def __init__(self, what: Path | str) -> None:
        super().__init__(f"Resolution is unset for {what}")
        self.what = what

    def __hash__(self) -> int:
        return hash((self.what,))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ResolutionUnsetError):
            return NotImplemented
        return self.what == other.what
