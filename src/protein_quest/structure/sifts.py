"""SIFTS mmCIF helpers.

These helpers operate on ``_pdbx_sifts_unp_segments`` category in mmCIF files,
which are used to map PDB chains to UniProt accessions and residue ranges.
"""

import logging
from operator import attrgetter

import gemmi

from protein_quest.uniprot_chains import UniprotChainMapping, UniprotChainMappings, UniprotChainRange

logger = logging.getLogger(__name__)


def uniprot_chain_mappings_from_cif(block: gemmi.cif.Block, best_only: bool = True) -> UniprotChainMappings:
    """Extract UniProt chain mappings from raw ``_pdbx_sifts_unp_segments`` rows.

    Args:
        block: mmCIF block containing ``_pdbx_sifts_unp_segments`` rows.
        best_only: If True, only return rows where ``best_mapping`` is ``y``.

    Returns:
        Set of UniProt chain mappings with ranges per chain. Empty if no SIFTS
        segment data is found. Returned chains use label system.
    """
    sifts_segments = block.get_mmcif_category("_pdbx_sifts_unp_segments.")
    if not sifts_segments:
        return set()

    acc_to_ranges: dict[str, list[UniprotChainRange]] = {}
    for i, acc in enumerate(sifts_segments["unp_acc"]):
        if best_only and sifts_segments["best_mapping"][i] != "y":
            continue
        chain_id = sifts_segments["asym_id"][i]
        try:
            start = int(sifts_segments["unp_start"][i])
            end = int(sifts_segments["unp_end"][i])
        except (ValueError, TypeError):
            logger.info(
                "Skipping pdbx_sifts_unp_segments row with accession %s in block %s "
                "due to non-numeric unp_start/unp_end",
                acc,
                block.name,
            )
            continue
        acc_to_ranges.setdefault(acc, []).append(UniprotChainRange(chain_ids=(chain_id,), start=start, end=end))

    return {
        UniprotChainMapping(uniprot_accession=acc, chain_ranges=tuple(ranges)) for acc, ranges in acc_to_ranges.items()
    }


def uniprot_chain_mappings_to_cif(mappings: UniprotChainMappings, block: gemmi.cif.Block) -> None:
    """Write minimal ``_pdbx_sifts_unp_segments`` rows to an mmCIF block.

    The emitted category only contains the subset consumed by
    ``uniprot_chain_mappings_from_cif()``.

    Args:
        mappings: UniProt chain mappings in label chain-id system.
        block: mmCIF block to update.
    """
    rows: dict[str, list[str]] = {
        "asym_id": [],
        "unp_acc": [],
        "unp_start": [],
        "unp_end": [],
        "best_mapping": [],
    }
    for mapping in sorted(mappings, key=attrgetter("uniprot_accession")):
        for chain_range in mapping.chain_ranges:
            for chain_id in chain_range.chain_ids:
                rows["asym_id"].append(chain_id)
                rows["unp_acc"].append(mapping.uniprot_accession)
                rows["unp_start"].append(str(chain_range.start))
                rows["unp_end"].append(str(chain_range.end))
                rows["best_mapping"].append("y")

    if rows["asym_id"]:
        block.set_mmcif_category("_pdbx_sifts_unp_segments.", rows)
