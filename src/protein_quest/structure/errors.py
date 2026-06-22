"""Exceptions for structure module operations."""

from pathlib import Path


class ChainNotFoundError(IndexError):
    """Exception raised when a chain is not found in a structure."""

    def __init__(self, chain_id: str, file: Path | str | None, available_chains: set[str]):
        super().__init__(f"Chain {chain_id} not found in {file}. Available chains are: {available_chains}")
        self.available_chains = available_chains
        self.chain_id = chain_id
        self.file = file

    def __reduce__(self):
        """Helper for pickling the exception."""
        return (self.__class__, (self.chain_id, self.file, self.available_chains))

    def __eq__(self, other):
        if not isinstance(other, ChainNotFoundError):
            return NotImplemented
        return (
            self.chain_id == other.chain_id
            and self.file == other.file
            and self.available_chains == other.available_chains
        )

    def __hash__(self):
        return hash((self.chain_id, str(self.file), frozenset(self.available_chains)))
