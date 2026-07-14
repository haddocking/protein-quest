"""UniProt extraction and injection helpers for structures."""

import logging
from collections import defaultdict
from collections.abc import Generator
from typing import Literal

import gemmi

from protein_quest.structure.chains import (
    ChainIdSystem,
    get_label2auth_chains,
    retrieve_chain_extraction_provenance,
)
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.types import Pdb2UniprotMapping, StructRefSeq

logger = logging.getLogger(__name__)


def struct_ref_seqs_columns_to_records(struct_ref_seqs_columns: dict[str, list[str | int]]) -> Generator[StructRefSeq]:
    """Convert `_struct_ref_seq` columns into collapsed alignment records.

    Args:
        struct_ref_seqs_columns: Struct ref seqs columns.
            As returned by `gemmi.Structure.make_mmcif_block().get_mmcif_category("_struct_ref_seq.")`.

    Yields:
        StructRefSeq records for each unique UniProt accession and chain combination."""
    struct_ref_seqs_rows = [
        dict(zip(struct_ref_seqs_columns.keys(), vals, strict=False))
        for vals in zip(*struct_ref_seqs_columns.values(), strict=False)
    ]

    grouped_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in struct_ref_seqs_rows:
        grouped_rows[(row["pdbx_db_accession"], row["pdbx_strand_id"])].append(row)

    for (accession, chain_id), rows in grouped_rows.items():
        starts = [int(row["db_align_beg"]) for row in rows]
        ends = [int(row["db_align_end"]) for row in rows]
        uniprot_start = min(starts)
        uniprot_end = max(ends)
        aligned_residue_count = sum(end - start + 1 for start, end in zip(starts, ends, strict=False))
        reference_span = uniprot_end - uniprot_start + 1
        sequence_identity = aligned_residue_count / reference_span
        yield StructRefSeq(
            uniprot_accession=accession,
            uniprot_start=uniprot_start,
            uniprot_end=uniprot_end,
            chain_id=chain_id,
            sequence_identity=sequence_identity,
            aligned_residue_count=aligned_residue_count,
        )


def _matching_struct_ref_seqs(structure: gemmi.Structure, accessions: set[str]) -> list[StructRefSeq]:
    block = structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
    struct_ref_seqs_columns = block.get_mmcif_category("_struct_ref_seq.")
    struct_ref_seqs = struct_ref_seqs_columns_to_records(struct_ref_seqs_columns)
    return [struct_ref_seq for struct_ref_seq in struct_ref_seqs if struct_ref_seq.uniprot_accession in accessions]


def _group_struct_ref_seqs_by_chain(struct_ref_seqs: list[StructRefSeq]) -> dict[str, list[StructRefSeq]]:
    grouped_struct_ref_seqs: dict[str, list[StructRefSeq]] = defaultdict(list)
    for struct_ref_seq in struct_ref_seqs:
        grouped_struct_ref_seqs[struct_ref_seq.chain_id].append(struct_ref_seq)
    return grouped_struct_ref_seqs


def _select_best_struct_ref_seq(struct_ref_seqs: list[StructRefSeq]) -> StructRefSeq:
    # Fast path: one UniProt for this chain needs no sorting.
    if len(struct_ref_seqs) == 1:
        return struct_ref_seqs[0]
    # Deterministic ranking: max aligned residues, then alphabetical accession.
    return min(struct_ref_seqs, key=lambda s: (-s.aligned_residue_count, s.uniprot_accession))


def selected_struct_ref_seqs_by_chain(structure: gemmi.Structure, accessions: set[str]) -> dict[str, StructRefSeq]:
    """Select one best ``StructRefSeq`` per chain for matching accessions.

    Selection for each chain is based on highest aligned residue count,
    with ties broken alphabetically by accession for deterministic behavior.

    Args:
        structure: The structure containing ``_struct_ref_seq`` records.
        accessions: UniProt accessions to match against.

    Returns:
        Mapping from chain id to selected ``StructRefSeq``.
        The key is in 'auth' [chain id system][protein_quest.structure.chains.ChainIdSystem].
        StructRefSeq.chain_id is also in 'auth' chain id system.
        Empty when no matching records are present.
    """
    matching_struct_ref_seqs = _matching_struct_ref_seqs(structure, accessions)
    struct_ref_seqs_by_chain = _group_struct_ref_seqs_by_chain(matching_struct_ref_seqs)
    return {
        chain_id: _select_best_struct_ref_seq(chain_struct_ref_seqs)
        for chain_id, chain_struct_ref_seqs in struct_ref_seqs_by_chain.items()
    }


def _extract_chain_uniprots_from_sifts(structure: gemmi.Structure) -> set[tuple[str, str]]:
    chain_uniprots: set[tuple[str, str]] = set()
    label2auth = get_label2auth_chains(structure)
    for entity in structure.entities:
        uniprots = entity.sifts_unp_acc
        if not uniprots:
            continue
        for label_subchain in entity.subchains:
            try:
                auth_subchain = label2auth[label_subchain]
            except KeyError:
                raise ChainNotFoundError(label_subchain, structure.name, set(label2auth.keys())) from None
            for uniprot in uniprots:
                chain_uniprots.add((auth_subchain, uniprot))
    return chain_uniprots


def _extract_chain_uniprots_from_struct_ref_seq(structure: gemmi.Structure) -> set[tuple[str, str]]:
    chain_uniprots: set[tuple[str, str]] = set()
    block = structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
    struct_ref = block.get_mmcif_category("_struct_ref.")
    struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
    if not (struct_ref and struct_ref_seq):
        logger.warning("No struct_ref data category found in structure %s", structure.name)
        return chain_uniprots

    db_name_unp_idx = [i for i, db_name in enumerate(struct_ref["db_name"]) if db_name == "UNP"]
    uniprot_accessions = {struct_ref["pdbx_db_accession"][i] for i in db_name_unp_idx}
    for i, uniprot_accession in enumerate(struct_ref_seq["pdbx_db_accession"]):
        if uniprot_accession in uniprot_accessions:
            auth_chain = struct_ref_seq["pdbx_strand_id"][i]
            chain_uniprots.add((auth_chain, uniprot_accession))

    if not chain_uniprots:
        logger.warning("No UniProt accessions found in structure %s", structure.name)
    return chain_uniprots


UniprotSource = Literal["both", "sifts", "struct_ref_seq", "fallback"]
"""From which source to extract UniProt accessions from a structure."""


def structure_to_uniprot(structure: gemmi.Structure, source: UniprotSource = "fallback") -> Pdb2UniprotMapping:
    """Extract chain-to-UniProt mapping from a structure.

    Args:
        structure: The gemmi Structure object to extract UniProt accessions from.
        source: UniProt source to read from.

            - ``sifts``: Read from entity ``sifts_unp_acc`` values.
            - ``struct_ref_seq``: Read from ``_struct_ref_seq`` filtered by
                ``_struct_ref`` records with ``db_name=UNP``.
            - ``both``: Merge SIFTS and struct_ref_seq results.
            - ``fallback``: Return SIFTS when available, otherwise ``struct_ref_seq``.

    Returns:
        A dictionary mapping the structure name to a set of tuples containing chain and UniProt accession.
            The chain names are in 'auth' [chain id system][protein_quest.structure.chains.ChainIdSystem].

    """
    if source == "fallback":
        sifts_uniprots = _extract_chain_uniprots_from_sifts(structure)
        if sifts_uniprots:
            return {structure.name: sifts_uniprots}
        return {structure.name: _extract_chain_uniprots_from_struct_ref_seq(structure)}

    chain_uniprots: set[tuple[str, str]] = set()
    if source in {"both", "sifts"}:
        chain_uniprots.update(_extract_chain_uniprots_from_sifts(structure))

    if source in {"both", "struct_ref_seq"}:
        chain_uniprots.update(_extract_chain_uniprots_from_struct_ref_seq(structure))

    return {structure.name: chain_uniprots}


def structure2uniprot_accessions(structure: gemmi.Structure) -> set[str]:
    """Extract UniProt accessions from a gemmi Structure object.

    Logs a warning and returns an empty set if no accessions are found in structure.

    Args:
        structure: The gemmi Structure object to extract UniProt accessions from.

    Returns:
        A set of UniProt accessions found in the structure."""
    s2cu = structure_to_uniprot(structure)
    cus = s2cu[structure.name]
    return {uniprot for _chain, uniprot in cus}


def _append_uniprot_to_structure(
    structure: gemmi.Structure, chain_uniprot_pairs: set[tuple[str, str]]
) -> gemmi.Structure:
    if not chain_uniprot_pairs:
        return structure
    block = structure.make_mmcif_block()
    struct_ref = block.get_mmcif_category("_struct_ref.")
    struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
    label2auth = get_label2auth_chains(structure)
    auth2label = {auth: label for label, auth in label2auth.items()}
    chain2entity_id: dict[str, str] = {
        chain: entity.name for entity in structure.entities for chain in entity.subchains
    }
    fillable_struct_ref_seq_cols = {"align_id", "ref_id", "pdbx_strand_id", "pdbx_db_accession"}
    for auth_chain, uniprot_accession in chain_uniprot_pairs:
        # Chain in chain_uniprot_pairs is of auth system, but entity subchains are label system
        # so we need to convert to label system to get entity id.
        try:
            label_chain = auth2label[auth_chain]
        except KeyError:
            raise ChainNotFoundError(auth_chain, structure.name, set(auth2label.keys())) from None
        entity_id = chain2entity_id[label_chain]
        new_id = str(len(struct_ref["id"]) + 1)
        struct_ref["id"].append(new_id)
        struct_ref["entity_id"].append(entity_id)
        struct_ref["db_name"].append("UNP")
        struct_ref["db_code"].append(None)
        struct_ref["pdbx_db_accession"].append(uniprot_accession)
        struct_ref["pdbx_db_isoform"].append(None)

        new_seq_id = len(struct_ref_seq["align_id"]) + 1
        struct_ref_seq["align_id"].append(new_seq_id)
        struct_ref_seq["ref_id"].append(new_id)
        struct_ref_seq["pdbx_strand_id"].append(auth_chain)
        struct_ref_seq["pdbx_db_accession"].append(uniprot_accession)
        for col in struct_ref_seq:
            if col not in fillable_struct_ref_seq_cols:
                struct_ref_seq[col].append(None)

    block.set_mmcif_category("_struct_ref.", struct_ref)
    block.set_mmcif_category("_struct_ref_seq.", struct_ref_seq)
    return gemmi.make_structure_from_block(block)


def _force_auth_system(
    structure: gemmi.Structure, pairs: set[tuple[str, str]], chain_system: ChainIdSystem
) -> set[tuple[str, str]]:
    if chain_system == "label":
        # Translate from label to auth chain ids
        # as all functions expect auth chain ids as input/output.
        label2auth = get_label2auth_chains(structure)
        try:
            return {(label2auth[label_chain], uniprot) for label_chain, uniprot in pairs}
        except KeyError as e:
            raise ChainNotFoundError(e.args[0], structure.name, set(label2auth.keys())) from None
    return pairs


def _rename_chain_based_on_provenance(structure: gemmi.Structure, pairs: set[tuple[str, str]]) -> set[tuple[str, str]]:
    prov = retrieve_chain_extraction_provenance(structure)

    if not prov:
        return pairs

    _, chain_provenance = prov
    logger.info(
        "Structure %s has provenance information indicating it was extracted from chain %s to %s. "
        "Using this information to verify/add UniProt accessions.",
        structure.name,
        chain_provenance.chain2keep,
        chain_provenance.out_chain,
    )
    return {(chain_provenance.out_chain, uniprot) for chain, uniprot in pairs if chain == chain_provenance.chain2keep}


def add_uniprot_accessions2structure(
    structure: gemmi.Structure,
    pdb2uniprot: Pdb2UniprotMapping | None,
    *,
    chain_system: ChainIdSystem = "auth",
) -> gemmi.Structure:
    """Add UniProt accessions to a structure if they are missing, based on the provided pdb2uniprot mapping.

    If structure has UniProt accession that is not in `pdb2uniprot`, it will be left unchanged.
    If structure has chain extraction provenance, the chain names from pdb2uniprot
    will be renamed to match the output chain name in the provenance.

    Args:
        structure: The gemmi Structure object to add UniProt accessions to.
        pdb2uniprot: Dictionary mapping PDB ID to set of tuples containing chain and UniProt accession.
            If provided, will be used to inject UniProt accessions into the structure if they are missing.
            If None, the structure is returned unchanged.
        chain_system: System of chain ids in ``pdb2uniprot`` mapping.

    Returns:
        A gemmi Structure object with UniProt accessions added if they were missing
        or the structure unchanged if all accessions were already present."""
    if not pdb2uniprot:
        return structure
    pdb_id = structure.name
    if pdb_id not in pdb2uniprot:
        logger.warning(
            "PDB ID %s not found in pdb2uniprot mapping. Leaving structure unverified and unchanged.", pdb_id
        )
        return structure

    expected_pairs = _force_auth_system(structure, pdb2uniprot[pdb_id], chain_system)
    expected_pairs = _rename_chain_based_on_provenance(structure, expected_pairs)

    known = structure_to_uniprot(structure)
    missing = expected_pairs - known[pdb_id]
    if not missing:
        return structure

    if known[pdb_id]:
        logger.warning(
            "Structure %s has some UniProt accessions that do not match the provided mapping. "
            "Existing: %s, Expected: %s, Missing: %s. Injecting missing accessions.",
            pdb_id,
            known[pdb_id],
            expected_pairs,
            missing,
        )
    else:
        logger.info("Injecting UniProt accessions into structure %s: %s", structure.name, missing)
    return _append_uniprot_to_structure(structure, missing)


def selected_struct_ref_seqs_from_sifts_by_chain(
    structure: gemmi.Structure, accessions: set[str]
) -> dict[str, StructRefSeq]:
    """Construct struct ref seq blocks from SIFTS data.

    With approximate sequence identity based on highest Uniprot position in structure.

    Args:
        structure: The structure to extract SIFTS data from.
        accessions: UniProt accessions to match against.

    Returns:
        Mapping from chain id to selected ``StructRefSeq``.
    """
    sc2ua = {}
    for entity in structure.entities:
        entity_uniprots = entity.sifts_unp_acc
        subchains = entity.subchains
        for subchain in subchains:
            sc2ua[subchain] = entity_uniprots

    # (chain, uniprot) -> list of uniprot positions in structure
    chain_ua2up: dict[tuple[str, str], list[int]] = {}
    for model in structure:
        for chain in model:
            polymer = chain.get_polymer()
            subchain_id = polymer.subchain_id()
            entity_uniprots = sc2ua[subchain_id]
            for residue in polymer:
                # valid sifts_unp looks like ('D', 0, 2865)
                sifts_unp: tuple[str, int, int] = residue.sifts_unp
                uniprot_pos = sifts_unp[2]
                if uniprot_pos == 0:
                    continue
                uniprot_accession = entity_uniprots[sifts_unp[1]]
                if uniprot_accession not in accessions:
                    continue
                key = (chain.name, uniprot_accession)
                if key not in chain_ua2up:
                    chain_ua2up[key] = []
                chain_ua2up[key].append(uniprot_pos)

    struct_ref_seqs_by_chain: dict[str, StructRefSeq] = {}
    for (chain_id, uniprot_accession), uniprot_positions in chain_ua2up.items():
        uniprot_start = min(uniprot_positions)
        uniprot_end = max(uniprot_positions)
        aligned_residue_count = len(uniprot_positions)
        # Do not know how long uniprot is so use highest mapped position as uniprot length
        reference_span = uniprot_end
        sequence_identity = aligned_residue_count / reference_span
        struct_ref_seqs_by_chain[chain_id] = StructRefSeq(
            uniprot_accession=uniprot_accession,
            uniprot_start=uniprot_start,
            uniprot_end=uniprot_end,
            chain_id=chain_id,
            sequence_identity=round(sequence_identity, 4),
            aligned_residue_count=aligned_residue_count,
        )

    return struct_ref_seqs_by_chain
