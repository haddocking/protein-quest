"""UniProt injection and provenance helpers for structures."""

import logging
from pathlib import Path

import gemmi

from protein_quest.structure.chains import (
    ChainExtractionProvenance,
    ChainIdSystem,
    _label_mappings_to_auth_system,
    get_label2auth_chains,
    retrieve_chain_extraction_provenance,
)
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.formats import read_structure, read_structure_as_cif_block
from protein_quest.structure.sifts import uniprot_chain_mappings_from_cif
from protein_quest.structure.uniprot_extraction import (
    ChainUniprotPair,
    uniprot_chain_mappings_from_struct_ref_seq,
)
from protein_quest.uniprot_chains import (
    Pdb2UniprotChainsMapping,
    UniprotChainMapping,
    UniprotChainMappings,
    UniprotChainRange,
    all_chain_ids,
)

logger = logging.getLogger(__name__)


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
    structure_file: Path,
    pdb2uniprot: Pdb2UniprotChainsMapping | None = None,
    *,
    structure: gemmi.Structure | None = None,
    chain_system: ChainIdSystem = "auth",
) -> tuple[gemmi.Structure, bool, UniprotChainMappings]:
    """Add UniProt accessions to a structure if they are missing, based on the provided pdb2uniprot mapping.

    If structure has UniProt accession that is not in `pdb2uniprot`, it will be left unchanged.
    If structure has chain extraction provenance, the chain names from pdb2uniprot
    will be renamed to match the output chain name in the provenance.

    Args:
        structure_file: Source mmCIF path to preserve raw SIFTS segment data
            when checking which UniProt mappings are already present.
        pdb2uniprot: Dictionary mapping PDB ID to structured UniProt chain mappings.
            If provided, will be used to inject UniProt accessions into the structure if they are missing.
            If None, the structure is returned unchanged.
        structure: The gemmi Structure object to add UniProt accessions to.
            Can be passed if caller already read structure.
            If not passed will use [protein_quest.structure.formats.read_structure][] to read `structure_file`.
        chain_system: System of chain ids in ``pdb2uniprot`` mapping.

    Returns:
        A tuple of (structure, injected, uniprot_mappings) where:
        - structure: gemmi Structure object with UniProt accessions added if they were missing
        - injected: bool indicating whether UniProt accessions were injected
        - uniprot_chain_mappings: set of UniprotChainMapping that were considered for injection (from pdb2uniprot),
            empty set if none
    """
    if structure is None:
        structure = read_structure(structure_file)
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

    block = read_structure_as_cif_block(structure_file)
    raw_sift_mappings = uniprot_chain_mappings_from_cif(block) if block is not None else set()
    sift_mappings = _label_mappings_to_auth_system(structure, raw_sift_mappings)
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
