"""Module for reformatting/filtering PDB results from UniProtKB SPARQL."""

import logging
from dataclasses import dataclass
from functools import cached_property

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PdbResult:
    """Result of a PDB search in UniProtKB.

    Attributes:
        id: PDB ID (for example "1H3O").
        method: Method used for the PDB entry (for example "X-ray diffraction").
        uniprot_chains: Chains in UniProt format (for example "A/B=1-42,A/B=50-99").
            The chain ids used in string are in 'auth' [chain id system][protein_quest.structure.chains.ChainIdSystem].
        resolution: Resolution of the PDB entry (for example "2.0" for 2.0 Å). Optional.
    """

    id: str
    method: str
    uniprot_chains: str
    resolution: str | None = None

    @cached_property
    def chain(self) -> str:
        """The first chain from the UniProt chains aka self.uniprot_chains.

        The chain ids used in string are in 'auth' [chain id system][protein_quest.structure.chains.ChainIdSystem].
        """
        return _first_chain_from_uniprot_chains(self.uniprot_chains)

    @cached_property
    def chain_length(self) -> int:
        """The length of the chain from the UniProt chains aka self.uniprot_chains.

        This is different from `self.uniprot_end - self.uniprot_start + 1`
        when there are multiple ranges in the UniProt chains string.

        This and end/start can be used to determine sequence identity
        (see 3D beacons API definition).
        """
        try:
            return _chain_length_from_uniprot_chains(self.uniprot_chains)
        except ValueError as e:
            raise PdbChainLengthError(self.id, self.uniprot_chains) from e

    @cached_property
    def uniprot_start(self) -> int:
        """Lowest UniProt residue position covered by this PDB result."""
        try:
            start, _ = _uniprot_start_end_from_uniprot_chains(self.uniprot_chains)
        except ValueError as e:
            raise PdbChainLengthError(self.id, self.uniprot_chains) from e
        else:
            return start

    @cached_property
    def uniprot_end(self) -> int:
        """Highest UniProt residue position covered by this PDB result."""
        try:
            _, end = _uniprot_start_end_from_uniprot_chains(self.uniprot_chains)
        except ValueError as e:
            raise PdbChainLengthError(self.id, self.uniprot_chains) from e
        else:
            return end

    @cached_property
    def sequence_identity(self) -> float:
        """Sequence identity of the PDB entry to the UniProt sequence.

        Calculated as (number of residues in the chain) / (length of the covered UniProt sequence).
        """
        try:
            return self.chain_length / (self.uniprot_end - self.uniprot_start + 1)
        except PdbChainLengthError:
            # If chain length cannot be determined, we cannot calculate sequence identity, return 0.0
            return 0.0

    @cached_property
    def resolution_value(self) -> float:
        """Resolution as a float for ordering.

        Missing or non-numeric resolutions become ``0.0`` so they rank as
        undesirable by resolution-based sorting.
        """
        value = self.resolution
        if value is None:
            return 0.0
        try:
            return float(value)
        except ValueError:
            return 0.0

    @cached_property
    def geometry_quality(self) -> float | None:
        """Geometry quality score (``0.0`` - ``100.0``) or ``None`` if unavailable."""
        return None


type PdbResults = dict[str, set[PdbResult]]
"""Dictionary with uniprot accessions as keys and sets of PDB results as values."""


def _chain_length_or_zero(entry: PdbResult) -> int:
    """Return chain length or 0 when the chain length cannot be determined."""
    try:
        return entry.chain_length
    except PdbChainLengthError:
        return 0


def _first_chain_from_uniprot_chains(uniprot_chains: str) -> str:
    """Extracts the first chain identifier from a UniProt chains string.

    The UniProt chains string is formatted (with EBNF notation) as follows:

        chain_group=range(,chain_group=range)*

    where:
        chain_group := chain_id(/chain_id)*
        chain_id    := [A-Za-z0-9]+
        range       := start-end
        start, end  := integer

    Args:
        uniprot_chains: A string representing UniProt chains, For example "B/D=1-81".

    Returns:
        The first chain identifier from the UniProt chain string. For example "B".
    """
    chains = uniprot_chains.split("=")
    parts = chains[0].split("/")
    chain = parts[0]
    try:
        # Workaround for Q9Y2Q5 │ 5YK3 │ 1/B/G=1-124, 1 does not exist but B does
        int(chain)
        if len(parts) > 1:
            return parts[1]
    except ValueError:
        # A letter
        pass
    return chain


def _chain_length_from_uniprot_chains(uniprot_chains: str) -> int:
    """Calculates the total length of chain from a UniProt chains string.

    See `_first_chain_from_uniprot_chains` for the format of the UniProt chains string.

    Args:
        uniprot_chains: A string representing UniProt chains, For example "B/D=1-81".

    Returns:
        The length of the chain in the UniProt chain string. For example 81 for "B/D=1-81".
    """
    total_length = 0
    chains = uniprot_chains.split(",")
    for chain in chains:
        _, rangestr = chain.split("=")
        start, stop = rangestr.split("-")
        # Residue positions are 1-based so + 1
        total_length += int(stop) - int(start) + 1
    return total_length


def _uniprot_start_end_from_uniprot_chains(uniprot_chains: str) -> tuple[int, int]:
    """Extract minimum start and maximum end residue positions from UniProt chains."""
    starts: list[int] = []
    ends: list[int] = []
    for chain in uniprot_chains.split(","):
        _, rangestr = chain.split("=")
        start, stop = rangestr.split("-")
        starts.append(int(start))
        ends.append(int(stop))

    return min(starts), max(ends)


class PdbChainLengthError(ValueError):
    """Raised when a UniProt chain description does not yield a chain length."""

    def __init__(self, pdb_id: str, uniprot_chains: str):
        msg = f"Could not determine chain length of '{pdb_id}' from '{uniprot_chains}'"
        super().__init__(msg)


def _sort_resolution_key(entry: PdbResult) -> tuple[int, float, int, str]:
    """Build a deterministic sort key for PDB resolution ranking.

    Lower resolution ranks first. Entries with missing or invalid resolution are
    undesirable and rank after entries with a real resolution. When resolution is
    tied, entries with more residues rank first, then PDB ID breaks ties.
    """
    resolution = entry.resolution_value
    chain_length = _chain_length_or_zero(entry)
    if resolution != 0.0:
        return (1, resolution, -chain_length, entry.id)
    return (2, 0.0, -chain_length, entry.id)


def filter_pdb_results_on_resolution(
    pdb_results: PdbResults,
    top: int,
) -> PdbResults:
    """Filter PDB results to top entries per UniProt accession by resolution.

    Entries are ranked by lower resolution first, then higher chain length,
    and finally deterministic PDB ID ordering.

    Args:
        pdb_results: Dictionary with UniProt accessions mapped to PDB entries.
        top: Maximum number of PDB entries to keep for each accession.

    Returns:
        Filtered dictionary with top-ranked entries per accession.
    """
    if top <= 0:
        msg = f"Top must be a positive integer, got {top}"
        raise ValueError(msg)
    results: PdbResults = {}
    for uniprot_accession, pdb_entries in pdb_results.items():
        ranked = sorted(pdb_entries, key=_sort_resolution_key)
        top_entries = set(ranked[:top])
        if top_entries:
            results[uniprot_accession] = top_entries

    return results


def filter_pdb_results_on_chain_length(
    pdb_results: PdbResults,
    min_residues: int | None,
    max_residues: int | None,
    keep_invalid: bool = False,
) -> PdbResults:
    """Filter PDB results based on chain length.

    Args:
        pdb_results: Dictionary with protein IDs as keys and sets of PDB results as values.
        min_residues: Minimum number of residues required in the chain mapped to the UniProt accession.
            If None, no minimum is applied.
        max_residues: Maximum number of residues allowed in chain mapped to the UniProt accession.
            If None, no maximum is applied.
        keep_invalid: If True, PDB results with invalid chain length (could not be determined) are kept.
            If False, PDB results with invalid chain length are filtered out.
            Warnings are logged when length can not be determined.

    Returns:
        Filtered dictionary with protein IDs as keys and sets of PDB results as values.
    """
    if min_residues is None and max_residues is None:
        # No filtering needed
        return pdb_results
    if min_residues is not None and max_residues is not None and max_residues <= min_residues:
        msg = f"Maximum number of residues ({max_residues}) must be > minimum number of residues ({min_residues})"
        raise ValueError(msg)
    results: PdbResults = {}
    for uniprot_accession, pdb_entries in pdb_results.items():
        filtered_pdb_entries = set()
        for pdb_entry in pdb_entries:
            try:
                if (min_residues is None or pdb_entry.chain_length >= min_residues) and (
                    max_residues is None or pdb_entry.chain_length <= max_residues
                ):
                    filtered_pdb_entries.add(pdb_entry)
            except PdbChainLengthError:
                if keep_invalid:
                    logger.warning(
                        f"Could not determine chain length of '{pdb_entry.id}' from '{pdb_entry.uniprot_chains}' "
                        f"belonging to uniprot accession '{uniprot_accession}', "
                        "for completeness not filtering it out"
                    )
                    filtered_pdb_entries.add(pdb_entry)
                else:
                    logger.warning(
                        f"Filtering out PDB entry '{pdb_entry.id}' belonging to uniprot accession "
                        f"'{uniprot_accession}' due to invalid chain length from '{pdb_entry.uniprot_chains}'"
                    )
        if filtered_pdb_entries:
            # Only include uniprot_accession if there are any pdb entries left after filtering
            results[uniprot_accession] = filtered_pdb_entries
    return results
