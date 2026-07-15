"""Helpers for parsing and formatting UniProt chain strings.

For example the Uniprot Sparql endpoint returns a string like ``A=1-100``
to show in which chain of a PDB entry the Uniprot accession is located.
"""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class UniprotChainRange:
    """One chain-group to UniProt residue-range assignment.

    Attributes:
        chain_ids: Chain ids in the group, in the order they appeared in the source string.
        start: Inclusive UniProt residue start.
        end: Inclusive UniProt residue end.
    """

    chain_ids: tuple[str, ...]
    start: int
    end: int

    def __post_init__(self):
        if not self.chain_ids:
            msg = "UniProt chain range must contain at least one chain id"
            raise ValueError(msg)
        if self.end < self.start:
            msg = f"UniProt chain range end {self.end} must be >= start {self.start}"
            raise ValueError(msg)

    def __len__(self) -> int:
        """Inclusive length of this UniProt residue range."""
        return self.end - self.start + 1

    @property
    def preferred_chain_id(self) -> str:
        """Preferred chain id for legacy single-chain consumers.

        Numeric leading chain ids are skipped when an alternative chain id exists.
        """
        first_chain = self.chain_ids[0]
        try:
            int(first_chain)
        except ValueError:
            return first_chain
        if len(self.chain_ids) > 1:
            return self.chain_ids[1]
        return first_chain


type UniprotChains = tuple[UniprotChainRange, ...]
"""Parsed UniProt chain-group/range assignments in source order."""


@dataclass(frozen=True)
class UniprotChainMapping:
    """UniProt accession plus parsed chain/range assignments for injection.

    Attributes:
        uniprot_accession: UniProt accession to inject.
        chain_ranges: Parsed chain/range assignments in either auth or label system,
            depending on the caller context.
    """

    uniprot_accession: str
    chain_ranges: UniprotChains


type Pdb2RangeMappings = dict[str, set[UniprotChainMapping]]
"""Dictionary mapping PDB ID to structured UniProt chain mappings used for injection.

Each mapping carries residue-range information for _struct_ref_seq injection.
"""


def parse_chain_ids(uniprot_chains: str) -> tuple[str, ...]:
    """Parse all unique chain ids from a UniProt chains string.

    Unlike ``parse_uniprot_chains()``, this helper only requires valid chain
    groups and tolerates invalid or missing ranges such as ``A=-``.

    Args:
        uniprot_chains: UniProt chains string such as ``A/B=1-100``.

    Returns:
        Chain ids in first-seen order with duplicates removed.

    Raises:
        ValueError: If the string does not contain at least one valid chain
            group or contains an empty chain id.
    """
    ordered_chain_ids: list[str] = []
    for chain_group in uniprot_chains.split(","):
        raw_chain_ids, _separator, _raw_range = chain_group.partition("=")
        if not raw_chain_ids:
            msg = "UniProt chains must contain at least one chain group"
            raise ValueError(msg)
        for chain_id in raw_chain_ids.split("/"):
            if not chain_id:
                msg = "UniProt chain group must contain at least one chain id"
                raise ValueError(msg)
            if chain_id not in ordered_chain_ids:
                ordered_chain_ids.append(chain_id)
    if not ordered_chain_ids:
        msg = "UniProt chains must contain at least one chain group"
        raise ValueError(msg)
    return tuple(ordered_chain_ids)


def preferred_chain_id(chain_ids: tuple[str, ...]) -> str:
    """Return the preferred chain id for legacy single-chain consumers.

    Numeric leading chain ids are skipped when a non-numeric alternative is
    available.

    Args:
        chain_ids: Chain ids in source order.

    Returns:
        Preferred chain id for consumers that only support a single chain.

    Raises:
        ValueError: If ``chain_ids`` is empty.
    """
    if not chain_ids:
        msg = "UniProt chains must contain at least one chain id"
        raise ValueError(msg)
    first_chain = chain_ids[0]
    try:
        int(first_chain)
    except ValueError:
        return first_chain
    if len(chain_ids) > 1:
        return chain_ids[1]
    return first_chain


def parse_uniprot_chains(uniprot_chains: str) -> UniprotChains:
    """Parse a UniProt chains string into chain-group/range assignments.

    The input format is ``chain_group=range(,chain_group=range)*`` where
    ``chain_group := chain_id(/chain_id)*`` and ``range := start-end``.

    Args:
        uniprot_chains: UniProt chains string such as
            ``A/B=2-459,A/B=520-610``.

    Returns:
        Parsed chain-group/range assignments in source order.

    Raises:
        ValueError: If the input contains an invalid or incomplete range, or if
            a parsed range is otherwise invalid.
    """
    parsed: list[UniprotChainRange] = []
    for chain_group in uniprot_chains.split(","):
        raw_chain_ids, raw_range = chain_group.split("=")
        start_str, end_str = raw_range.split("-")
        parsed.append(
            UniprotChainRange(
                chain_ids=tuple(raw_chain_ids.split("/")),
                start=int(start_str),
                end=int(end_str),
            )
        )
    return tuple(parsed)


def format_uniprot_chains(chain_ranges: UniprotChains) -> str:
    """Format parsed UniProt chain assignments back into a string.

    Args:
        chain_ranges: Parsed UniProt chain assignments.

    Returns:
        String representation such as ``A/B=2-459,C=520-610``.
    """
    return ",".join(
        f"{'/'.join(chain_range.chain_ids)}={chain_range.start}-{chain_range.end}" for chain_range in chain_ranges
    )


def first_chain_id(chain_ranges: UniprotChains) -> str:
    """Return the preferred chain id from parsed UniProt chains.

    Args:
        chain_ranges: Parsed UniProt chain assignments.

    Returns:
        Preferred chain id from the first chain group.

    Raises:
        ValueError: If ``chain_ranges`` is empty.
    """
    if not chain_ranges:
        msg = "UniProt chains must contain at least one chain group"
        raise ValueError(msg)
    return preferred_chain_id(chain_ranges[0].chain_ids)


def all_chain_ids(chain_ranges: UniprotChains) -> tuple[str, ...]:
    """Return all unique chain ids in first-seen order.

    Args:
        chain_ranges: Parsed UniProt chain assignments.

    Returns:
        Unique chain ids in the order they first appear.
    """
    return parse_chain_ids(format_uniprot_chains(chain_ranges))


def chain_length(chain_ranges: UniprotChains) -> int:
    """Return the total aligned UniProt residue count across all ranges.

    Args:
        chain_ranges: Parsed UniProt chain assignments.

    Returns:
        Sum of inclusive lengths across all ranges.
    """
    return sum(len(chain_range) for chain_range in chain_ranges)


def covered_uniprot_range(chain_ranges: UniprotChains) -> tuple[int, int]:
    """Return the minimum start and maximum end residue positions.

    Args:
        chain_ranges: Parsed UniProt chain assignments.

    Returns:
        Tuple of ``(minimum_start, maximum_end)``.

    Raises:
        ValueError: If ``chain_ranges`` is empty.
    """
    if not chain_ranges:
        msg = "UniProt chains must contain at least one chain group"
        raise ValueError(msg)
    starts = [chain_range.start for chain_range in chain_ranges]
    ends = [chain_range.end for chain_range in chain_ranges]
    return min(starts), max(ends)
