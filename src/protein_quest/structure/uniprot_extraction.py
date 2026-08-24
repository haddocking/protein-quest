"""UniProt extraction and mapping helpers for structures."""

import logging
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import gemmi

from protein_quest.structure.chains import _label_mappings_to_auth_system
from protein_quest.structure.formats import read_structure, read_structure_as_cif_block
from protein_quest.structure.sifts import uniprot_chain_mappings_from_cif
from protein_quest.uniprot_chains import (
    UniprotChainMapping,
    UniprotChainMappings,
    UniprotChainRange,
)

logger = logging.getLogger(__name__)


def uniprot_chain_mappings_from_struct_ref_seq(structure: gemmi.Structure) -> UniprotChainMappings:
    """Extract UniProt chain mappings from `_struct_ref_seq` rows.

    Args:
        structure: The structure containing ``_struct_ref`` and ``_struct_ref_seq`` records.

    Returns:
        Set of UniProt chain mappings with ranges per chain. Empty if no UNP data found."""
    block = structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
    struct_ref = block.get_mmcif_category("_struct_ref.")
    struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")

    if not struct_ref or not struct_ref_seq:
        return set()

    unp_indices = [i for i, db_name in enumerate(struct_ref["db_name"]) if db_name == "UNP"]
    unp_accessions = {struct_ref["pdbx_db_accession"][i] for i in unp_indices}
    if not unp_accessions:
        return set()

    acc_to_ranges: dict[str, list[UniprotChainRange]] = defaultdict(list)
    for i, acc in enumerate(struct_ref_seq["pdbx_db_accession"]):
        if acc in unp_accessions:
            chain_id = struct_ref_seq["pdbx_strand_id"][i]
            try:
                beg = int(struct_ref_seq["db_align_beg"][i])
                end = int(struct_ref_seq["db_align_end"][i])
            except (ValueError, TypeError):
                logger.info(
                    "Skipping struct_ref_seq row with align_id %s of %s due to non-numeric db_align_beg/db_align_end",
                    struct_ref_seq["align_id"][i],
                    acc,
                )
                continue
            acc_to_ranges[acc].append(
                UniprotChainRange(
                    chain_ids=(chain_id,),
                    start=beg,
                    end=end,
                )
            )

    return {
        UniprotChainMapping(uniprot_accession=acc, chain_ranges=tuple(ranges)) for acc, ranges in acc_to_ranges.items()
    }


ChainUniprotPair = namedtuple("ChainUniprotPair", ["chain_id", "uniprot_accession"])
"""Pair of chain id and UniProt accession for mapping purposes."""


@dataclass(frozen=True, slots=True)
class FlattenedUniprotChainMapping:
    """Collapsed `_struct_ref_seq` like alignment information for one chain.

    Attributes:

        uniprot_accession: The UniProt accession.
        uniprot_start: The start position of the alignment on the UniProt sequence.
        uniprot_end: The end position of the alignment on the UniProt sequence.
        chain_id: The chain ID in the 'auth' [chain ID system][protein_quest.structure.chains.ChainIdSystem].
        sequence_identity: The sequence identity of the alignment.
        aligned_residue_count: The number of aligned residues in the alignment.
    """

    uniprot_accession: str
    uniprot_start: int
    uniprot_end: int
    chain_id: str
    sequence_identity: float
    aligned_residue_count: int


def flatten_uniprot_chain_mappings(mappings: UniprotChainMappings) -> set[FlattenedUniprotChainMapping]:
    """Flatten a set of UniprotChainMapping.

    Each (accession, chain) group is collapsed into one record with merged start/end
    and ``aligned_residue_count`` summed across all ranges.

    Args:
        mappings: Set of UniprotChainMapping.

    Returns:
        Set of flattened per-(accession, chain) records with merged start/end,
        summed aligned residue counts, and computed sequence identity.
    """
    groups: dict[tuple[str, str], list[UniprotChainRange]] = defaultdict(list)
    for mapping in mappings:
        for chain_range in mapping.chain_ranges:
            chain_id = chain_range.preferred_chain_id
            groups[(mapping.uniprot_accession, chain_id)].append(chain_range)

    result: set[FlattenedUniprotChainMapping] = set()
    for (accession, chain_id), ranges in groups.items():
        starts = [r.start for r in ranges]
        ends = [r.end for r in ranges]
        uniprot_start = min(starts)
        uniprot_end = max(ends)
        aligned_residue_count = sum(len(r) for r in ranges)
        reference_span = uniprot_end - uniprot_start + 1
        sequence_identity = aligned_residue_count / reference_span
        result.add(
            FlattenedUniprotChainMapping(
                uniprot_accession=accession,
                uniprot_start=uniprot_start,
                uniprot_end=uniprot_end,
                chain_id=chain_id,
                sequence_identity=sequence_identity,
                aligned_residue_count=aligned_residue_count,
            )
        )
    return result


def best_uniprot_per_chain(mappings: set[FlattenedUniprotChainMapping]) -> set[FlattenedUniprotChainMapping]:
    groups: dict[str, list[FlattenedUniprotChainMapping]] = defaultdict(list)
    for mapping in mappings:
        groups[mapping.chain_id].append(mapping)

    result: set[FlattenedUniprotChainMapping] = set()
    for chain_seqs in groups.values():
        if len(chain_seqs) == 1:
            result.add(chain_seqs[0])
        else:
            result.add(min(chain_seqs, key=_best_uniprot_sort_key))
    return result


def _best_uniprot_sort_key(mapping: FlattenedUniprotChainMapping) -> tuple[int, str]:
    return (-mapping.aligned_residue_count, mapping.uniprot_accession)


UniprotSource = Literal["both", "sifts", "struct_ref_seq", "fallback"]
"""From which source to extract UniProt accessions from a structure."""


def structure_to_uniprot(
    structure_file: Path,
    source: UniprotSource = "both",
    one_uniprot_per_chain: bool = True,
    structure: gemmi.Structure | None = None,
) -> set[FlattenedUniprotChainMapping]:
    """Extract UniProt chain mappings from a structure.

    Args:
        structure_file: Source mmCIF path. When provided for ``.cif`` or
            ``.cif.gz`` files, SIFTS mappings are extracted from raw
            ``_pdbx_sifts_unp_segments`` rows as structure does not contain that info.
        source: UniProt source to read from.

            - ``sifts``: Read from entity ``sifts_unp_acc`` values.
            - ``struct_ref_seq``: Read from ``_struct_ref_seq`` filtered by
                ``_struct_ref`` records with ``db_name=UNP``.
            - ``both``: Merge SIFTS and struct_ref_seq results.
            - ``fallback``: Return SIFTS when available, otherwise ``struct_ref_seq``.
        one_uniprot_per_chain: If True, return only the best UniProt per chain,
            based on highest aligned residue count, with ties broken alphabetically by accession.
            Otherwise, return all UniProt mappings for each chain.
        structure: The structure containing SIFTS and/or ``_struct_ref_seq`` data.
            Can be passed if caller already read structure.
            If not passed will use [protein_quest.structure.formats.read_structure][] to read `structure_file`.

    Returns:
        Set of flattened per-(accession, chain) records with merged start/end,
        summed aligned residue counts, and computed sequence identity.
    """
    if structure is None:
        structure = read_structure(structure_file)

    sift_mappings: set[FlattenedUniprotChainMapping] = set()
    if source in {"sifts", "both", "fallback"}:
        block = read_structure_as_cif_block(structure_file)
        raw_sift_mappings = uniprot_chain_mappings_from_cif(block) if block is not None else set()
        sift_mappings = flatten_uniprot_chain_mappings(_label_mappings_to_auth_system(structure, raw_sift_mappings))

    struct_ref_mappings: set[FlattenedUniprotChainMapping] = set()
    if source in {"struct_ref_seq", "both", "fallback"}:
        struct_ref_mappings = flatten_uniprot_chain_mappings(uniprot_chain_mappings_from_struct_ref_seq(structure))

    mappings: set[FlattenedUniprotChainMapping] = set()
    if source == "sifts":
        mappings = sift_mappings
    elif source == "struct_ref_seq":
        mappings = struct_ref_mappings
    elif source == "both":
        mappings = sift_mappings | struct_ref_mappings
    elif source == "fallback":
        mappings = sift_mappings or struct_ref_mappings
    else:
        msg = f"Invalid source '{source}', must be one of 'sifts', 'struct_ref_seq', 'both', or 'fallback'"
        raise ValueError(msg)

    if one_uniprot_per_chain:
        return best_uniprot_per_chain(mappings)
    return mappings


def structure2uniprot_accessions(structure_file: Path) -> set[str]:
    """Extract UniProt accessions from a gemmi Structure object.

    Logs a warning and returns an empty set if no accessions are found in structure.

    Args:
        structure_file: Source mmCIF path to preserve raw SIFTS segment data.

    Returns:
        A set of UniProt accessions found in the structure.
    """
    mappings = structure_to_uniprot(structure_file, one_uniprot_per_chain=False)
    return {mapping.uniprot_accession for mapping in mappings}
