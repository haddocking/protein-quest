"""UniProt extraction and injection helpers for structures."""

import logging
from collections import defaultdict, namedtuple
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

import gemmi

from protein_quest.structure.chains import (
    ChainExtractionProvenance,
    ChainIdSystem,
    get_label2auth_chains,
    retrieve_chain_extraction_provenance,
)
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.uniprot_chains import (
    Pdb2UniprotChainsMapping,
    UniprotChainMapping,
    UniprotChainMappings,
    UniprotChainRange,
    all_chain_ids,
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


def _positions_to_ranges(positions: Iterable[int]) -> list[tuple[int, int]]:
    if not positions:
        return []

    sorted_positions = sorted(positions)
    ranges = []
    start = sorted_positions[0]
    end = sorted_positions[0]

    for pos in sorted_positions[1:]:
        if pos == end + 1:
            end = pos
        else:
            ranges.append((start, end))
            start = pos
            end = pos

    ranges.append((start, end))
    return ranges


def _subchains2sifts_unp_acc(structure: gemmi.Structure) -> dict[str, list[str]]:
    sc2ua = {}
    for entity in structure.entities:
        entity_uniprots = entity.sifts_unp_acc
        subchains = entity.subchains
        for subchain in subchains:
            sc2ua[subchain] = entity_uniprots
    return sc2ua


def uniprot_chain_mappings_from_sifts(structure: gemmi.Structure) -> UniprotChainMappings:
    """Extract UniProt chain mappings from SIFTS data.

    Args:
        structure: The structure containing SIFTS data.

    Returns:
        Set of UniProt chain mappings with ranges per chain. Empty if no SIFTS data found."""
    sc2ua = _subchains2sifts_unp_acc(structure)

    chain_ua2up: dict[ChainUniprotPair, set[int]] = defaultdict(set)
    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            subchain_id = polymer.subchain_id()
            entity_uniprots = sc2ua.get(subchain_id)
            if not entity_uniprots:
                continue
            for residue in polymer:
                sifts_unp: tuple[str, int, int] = residue.sifts_unp
                uniprot_pos = sifts_unp[2]
                if uniprot_pos == 0:
                    continue
                uniprot_accession = entity_uniprots[sifts_unp[1]]
                key = ChainUniprotPair(chain.name, uniprot_accession)
                chain_ua2up[key].add(uniprot_pos)

    acc_to_ranges: dict[str, list[UniprotChainRange]] = defaultdict(list)
    for (chain_id, uniprot_accession), positions in chain_ua2up.items():
        ranges = _positions_to_ranges(positions)
        for start, end in ranges:
            acc_to_ranges[uniprot_accession].append(
                UniprotChainRange(
                    chain_ids=(chain_id,),
                    start=start,
                    end=end,
                )
            )

    return {
        UniprotChainMapping(uniprot_accession=acc, chain_ranges=tuple(ranges)) for acc, ranges in acc_to_ranges.items()
    }


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
    for s in mappings:
        groups[s.chain_id].append(s)

    result: set[FlattenedUniprotChainMapping] = set()
    for chain_seqs in groups.values():
        if len(chain_seqs) == 1:
            result.add(chain_seqs[0])
        else:
            result.add(min(chain_seqs, key=lambda s: (-s.aligned_residue_count, s.uniprot_accession)))
    return result


UniprotSource = Literal["both", "sifts", "struct_ref_seq", "fallback"]
"""From which source to extract UniProt accessions from a structure."""


def structure_to_uniprot(
    structure: gemmi.Structure, source: UniprotSource = "both", one_uniprot_per_chain: bool = True
) -> set[FlattenedUniprotChainMapping]:
    """Extract UniProt chain mappings from a structure.

    Args:
        structure: The structure containing SIFTS and/or ``_struct_ref_seq`` data.
        source: UniProt source to read from.

            - ``sifts``: Read from entity ``sifts_unp_acc`` values.
            - ``struct_ref_seq``: Read from ``_struct_ref_seq`` filtered by
                ``_struct_ref`` records with ``db_name=UNP``.
            - ``both``: Merge SIFTS and struct_ref_seq results.
            - ``fallback``: Return SIFTS when available, otherwise ``struct_ref_seq``.
        one_uniprot_per_chain: If True, return only the best UniProt per chain,
            based on highest aligned residue count, with ties broken alphabetically by accession.
            Otherwise, return all UniProt mappings for each chain.

    Returns:
        Set of flattened per-(accession, chain) records with merged start/end,
        summed aligned residue counts, and computed sequence identity.
    """
    sift_mappings = flatten_uniprot_chain_mappings(uniprot_chain_mappings_from_sifts(structure))
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


def structure2uniprot_accessions(structure: gemmi.Structure) -> set[str]:
    """Extract UniProt accessions from a gemmi Structure object.

    Logs a warning and returns an empty set if no accessions are found in structure.

    Args:
        structure: The gemmi Structure object to extract UniProt accessions from.

    Returns:
        A set of UniProt accessions found in the structure."""
    mappings = structure_to_uniprot(structure, one_uniprot_per_chain=False)
    return {m.uniprot_accession for m in mappings}


def _append_uniprot_to_structure(structure: gemmi.Structure, chain_mappings: UniprotChainMappings) -> gemmi.Structure:
    block = structure.make_mmcif_block()
    struct_ref = block.get_mmcif_category("_struct_ref.")
    struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
    label2auth = get_label2auth_chains(structure)
    auth2label = {auth: label for label, auth in label2auth.items()}
    chain2entity_id: dict[str, str] = {
        chain: entity.name for entity in structure.entities for chain in entity.subchains
    }
    fillable_struct_ref_seq_cols = {
        "align_id",
        "ref_id",
        "pdbx_strand_id",
        "pdbx_db_accession",
        "db_align_beg",
        "db_align_end",
    }
    for mapping in chain_mappings:
        for chain_range in mapping.chain_ranges:
            for auth_chain in chain_range.chain_ids:
                # Chain in the mapping is of auth system, but entity subchains are label system
                # so we need to convert to label system to get entity id.
                try:
                    label_chain = auth2label[auth_chain]
                except KeyError:
                    raise ChainNotFoundError(auth_chain, structure.name, set(auth2label.keys())) from None
                entity_id = chain2entity_id[label_chain]
                new_id = str(len(struct_ref["id"]) + 1)
                # struct_ref must be same length as struct_ref_seq, as struct_ref has pdbx_align_begin field,
                # but its ignored by gemmi
                struct_ref["id"].append(new_id)
                struct_ref["entity_id"].append(entity_id)
                struct_ref["db_name"].append("UNP")
                struct_ref["db_code"].append(None)
                struct_ref["pdbx_db_accession"].append(mapping.uniprot_accession)
                struct_ref["pdbx_db_isoform"].append(None)

                new_seq_id = len(struct_ref_seq["align_id"]) + 1
                struct_ref_seq["align_id"].append(new_seq_id)
                struct_ref_seq["ref_id"].append(new_id)
                struct_ref_seq["pdbx_strand_id"].append(auth_chain)
                struct_ref_seq["pdbx_db_accession"].append(mapping.uniprot_accession)
                struct_ref_seq["db_align_beg"].append(chain_range.start)
                struct_ref_seq["db_align_end"].append(chain_range.end)
                for col, values in struct_ref_seq.items():
                    if col not in fillable_struct_ref_seq_cols:
                        values.append(None)

    block.set_mmcif_category("_struct_ref.", struct_ref)
    block.set_mmcif_category("_struct_ref_seq.", struct_ref_seq)
    return gemmi.make_structure_from_block(block)


def _force_auth_system(
    structure: gemmi.Structure, mappings: UniprotChainMappings, chain_system: ChainIdSystem
) -> UniprotChainMappings:
    if chain_system == "label":
        # Translate from label to auth chain ids
        # as all functions expect auth chain ids as input/output.
        label2auth = get_label2auth_chains(structure)
        try:
            return {
                UniprotChainMapping(
                    uniprot_accession=mapping.uniprot_accession,
                    chain_ranges=tuple(
                        UniprotChainRange(
                            chain_ids=tuple(label2auth[label_chain] for label_chain in chain_range.chain_ids),
                            start=chain_range.start,
                            end=chain_range.end,
                        )
                        for chain_range in mapping.chain_ranges
                    ),
                )
                for mapping in mappings
            }
        except KeyError as e:
            raise ChainNotFoundError(e.args[0], structure.name, set(label2auth.keys())) from None
    return mappings


def _mapping_pairs(mappings: UniprotChainMappings) -> set[ChainUniprotPair]:
    return {
        ChainUniprotPair(chain_id, mapping.uniprot_accession)
        for mapping in mappings
        for chain_id in all_chain_ids(mapping.chain_ranges)
    }


def apply_chain_provenance_to_uniprot_mappings(
    mappings: UniprotChainMappings, chain_provenance: ChainExtractionProvenance
) -> UniprotChainMappings:
    """Apply chain extraction provenance to a set of UniprotChainMapping.

    Args:
        mappings: Set of UniprotChainMapping to apply provenance to.
        chain_provenance: ChainExtractionProvenance to apply to mappings.

    Returns:
        Set of UniprotChainMapping with chain ids updated based on provenance.
    """
    renamed_mappings: UniprotChainMappings = set()
    for mapping in mappings:
        renamed_ranges = []
        for chain_range in mapping.chain_ranges:
            renamed_chain_ids = tuple(
                chain_id for chain_id in chain_range.chain_ids if chain_id != chain_provenance.chain2keep
            )
            if len(renamed_chain_ids) != len(chain_range.chain_ids):
                renamed_chain_ids = tuple(sorted((chain_provenance.out_chain, *renamed_chain_ids)))
            renamed_ranges.append(
                UniprotChainRange(
                    chain_ids=renamed_chain_ids,
                    start=chain_range.start,
                    end=chain_range.end,
                )
            )
        renamed_mappings.add(
            UniprotChainMapping(
                uniprot_accession=mapping.uniprot_accession,
                chain_ranges=tuple(renamed_ranges),
            )
        )
    return renamed_mappings


def _rename_chain_based_on_provenance(
    structure: gemmi.Structure, mappings: UniprotChainMappings
) -> UniprotChainMappings:
    prov = retrieve_chain_extraction_provenance(structure)

    if not prov:
        return mappings

    _, chain_provenance = prov
    logger.info(
        "Structure %s has provenance information indicating it was extracted from chain %s to %s. "
        "Using this information to verify/add UniProt accessions.",
        structure.name,
        chain_provenance.chain2keep,
        chain_provenance.out_chain,
    )
    return apply_chain_provenance_to_uniprot_mappings(mappings, chain_provenance)


def _filter_mappings_by_pairs(mappings: UniprotChainMappings, pairs: set[ChainUniprotPair]) -> UniprotChainMappings:
    filtered_mappings: UniprotChainMappings = set()
    for mapping in mappings:
        filtered_ranges = []
        for chain_range in mapping.chain_ranges:
            matching_chain_ids = tuple(
                chain_id for chain_id in chain_range.chain_ids if (chain_id, mapping.uniprot_accession) in pairs
            )
            if matching_chain_ids:
                filtered_ranges.append(
                    UniprotChainRange(
                        chain_ids=matching_chain_ids,
                        start=chain_range.start,
                        end=chain_range.end,
                    )
                )
        if filtered_ranges:
            filtered_mappings.add(
                UniprotChainMapping(
                    uniprot_accession=mapping.uniprot_accession,
                    chain_ranges=tuple(filtered_ranges),
                )
            )
    return filtered_mappings


def add_uniprot_accessions2structure(
    structure: gemmi.Structure,
    pdb2uniprot: Pdb2UniprotChainsMapping | None,
    *,
    chain_system: ChainIdSystem = "auth",
) -> tuple[gemmi.Structure, bool, UniprotChainMappings]:
    """Add UniProt accessions to a structure if they are missing, based on the provided pdb2uniprot mapping.

    If structure has UniProt accession that is not in `pdb2uniprot`, it will be left unchanged.
    If structure has chain extraction provenance, the chain names from pdb2uniprot
    will be renamed to match the output chain name in the provenance.

    Args:
        structure: The gemmi Structure object to add UniProt accessions to.
        pdb2uniprot: Dictionary mapping PDB ID to structured UniProt chain mappings.
            If provided, will be used to inject UniProt accessions into the structure if they are missing.
            If None, the structure is returned unchanged.
        chain_system: System of chain ids in ``pdb2uniprot`` mapping.

    Returns:
        A tuple of (structure, injected, uniprot_mappings) where:
        - structure: gemmi Structure object with UniProt accessions added if they were missing
        - injected: bool indicating whether UniProt accessions were injected
        - uniprot_chain_mappings: set of UniprotChainMapping that were considered for injection (from pdb2uniprot),
          empty set if none
    """
    if not pdb2uniprot:
        return structure, False, set()
    pdb_id = structure.name
    if pdb_id not in pdb2uniprot:
        logger.warning(
            "PDB ID %s not found in pdb2uniprot mapping. Leaving structure unverified and unchanged.", pdb_id
        )
        return structure, False, set()

    expected_mappings = _force_auth_system(structure, pdb2uniprot[pdb_id], chain_system)
    expected_mappings = _rename_chain_based_on_provenance(structure, expected_mappings)
    expected_pairs = _mapping_pairs(expected_mappings)

    sift_mappings = uniprot_chain_mappings_from_sifts(structure)
    struct_ref_mappings = uniprot_chain_mappings_from_struct_ref_seq(structure)
    known_mappings = sift_mappings | struct_ref_mappings
    known_pairs = _mapping_pairs(known_mappings)

    missing = expected_pairs - known_pairs
    if not missing:
        return structure, False, set()

    if known_pairs:
        logger.warning(
            "Structure %s has some UniProt accessions that do not match the provided mapping. "
            "Existing: %s, Expected: %s, Missing: %s. Injecting missing accessions.",
            pdb_id,
            known_pairs,
            expected_pairs,
            missing,
        )
    else:
        logger.info("Injecting UniProt accessions into structure %s: %s", structure.name, missing)

    missing_mappings = _filter_mappings_by_pairs(expected_mappings, missing)
    new_structure = _append_uniprot_to_structure(structure, missing_mappings)
    return new_structure, True, missing_mappings
