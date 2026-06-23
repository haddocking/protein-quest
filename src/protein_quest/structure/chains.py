"""Chain-level structure helpers and transformations."""

import logging
from datetime import UTC, datetime
from pathlib import Path

import gemmi
from gemmi import Structure

from protein_quest.__version__ import __version__
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.files import split_name_and_extension
from protein_quest.structure.formats import read_structure, write_structure
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


def find_chain_in_model(model: gemmi.Model, wanted_chain: str) -> gemmi.Chain | None:
    """Find a chain in a model.

    Args:
        model: The gemmi model to search in.
        wanted_chain: The chain identifier to search for.

    Returns:
        The found chain or None if not found."""
    chain = model.find_chain(wanted_chain)
    if chain is None:
        mchains = [c for c in model if c.name.endswith(wanted_chain)]
        if mchains:
            return mchains[0]
    return chain


def find_chain_in_structure(structure: gemmi.Structure, wanted_chain: str) -> gemmi.Chain | None:
    """Find a chain in a structure.

    Args:
        structure: The gemmi structure to search in.
        wanted_chain: The chain identifier to search for.

    Returns:
        The found chain or None if not found."""
    for model in structure:
        chain = find_chain_in_model(model, wanted_chain)
        if chain is not None:
            return chain
    return None


def nr_residues_in_chain(file: Path, chain: str = "A") -> int:
    """Returns the number of residues in a specific chain from a structure file.

    Args:
        file: Path to the input structure file.
        chain: Chain to count residues of.

    Returns:
        The number of residues in the specified chain."""
    structure = read_structure(file)
    gchain = find_chain_in_structure(structure, chain)
    if gchain is None:
        logger.warning("Chain %s not found in %s. Returning 0.", chain, file)
        return 0
    return len(gchain)


def _dedup_helices(structure: gemmi.Structure):
    helix_starts: set[str] = set()
    duplicate_helix_indexes: list[int] = []
    for hindex, helix in enumerate(structure.helices):
        if str(helix.start) in helix_starts:
            logger.debug("Duplicate start helix found: %s %s, removing", hindex, helix.start)
            duplicate_helix_indexes.append(hindex)
        else:
            helix_starts.add(str(helix.start))
    for helix_index in reversed(duplicate_helix_indexes):
        structure.helices.pop(helix_index)


def _dedup_sheets(structure: gemmi.Structure, chain2keep: str):
    duplicate_sheet_indexes: list[int] = []
    for sindex, sheet in enumerate(structure.sheets):
        if sheet.name != chain2keep:
            duplicate_sheet_indexes.append(sindex)
    for sheet_index in reversed(duplicate_sheet_indexes):
        structure.sheets.pop(sheet_index)


def _add_provenance_info(structure: gemmi.Structure, chain2keep: str, out_chain: str):
    old_id = structure.name
    new_id = structure.name + f"{chain2keep}2{out_chain}"
    structure.name = new_id
    structure.info["_entry.id"] = new_id
    new_title = f"From {old_id} chain {chain2keep} to {out_chain}"
    structure.info["_struct.title"] = new_title
    structure.info["_struct_keywords.pdbx_keywords"] = new_title.upper()
    new_si = gemmi.SoftwareItem()
    new_si.classification = gemmi.SoftwareItem.Classification.DataExtraction
    new_si.name = "protein-quest.pdbe.io.write_single_chain_pdb_file"
    new_si.version = __version__
    new_si.date = str(datetime.now(tz=UTC).date())
    structure.meta.software = [*structure.meta.software, new_si]


def chains_in_structure(structure: gemmi.Structure) -> set[gemmi.Chain]:
    """Get a list of chains in a structure.

    Args:
        structure: The gemmi structure to get chains from.

    Returns:
        A set of chains in the structure."""
    return {c for model in structure for c in model}


def write_single_chain_structure_file(
    input_file: Path,
    chain2keep: str,
    output_dir: Path,
    out_chain: str = "A",
    copy_method: CopyMethod = "copy",
) -> Path:
    """Write a single chain from a structure file to a new structure file.

    Also

    - removes ligands and waters
    - renumbers atoms ids
    - removes chem_comp section from cif files
    - adds provenance information to the header like software and input file+chain

    This function is equivalent to the following gemmi commands:

    ```shell
    gemmi convert --remove-lig-wat --select=B --to=cif chain-in/3JRS.cif - | \
    gemmi convert --from=cif --rename-chain=B:A - chain-out/3JRS_B2A.gemmi.cif
    ```

    Args:
        input_file: Path to the input structure file.
        chain2keep: The chain to keep.
        output_dir: Directory to save the output file.
        out_chain: The chain identifier for the output file.
        copy_method: How to copy when no changes are needed to output file.

    Returns:
        Path to the output structure file

    Raises:
        FileNotFoundError: If the input file does not exist.
        ChainNotFoundError: If the specified chain is not found in the input file."""
    logger.debug("chain2keep: %s, out_chain: %s", chain2keep, out_chain)
    structure = read_structure(input_file)
    structure.setup_entities()

    chain = find_chain_in_structure(structure, chain2keep)
    chainnames_in_structure = {c.name for c in chains_in_structure(structure)}
    if chain is None:
        raise ChainNotFoundError(chain2keep, input_file, chainnames_in_structure)
    chain_name = chain.name
    name, extension = split_name_and_extension(input_file.name)
    output_file = output_dir / f"{name}_{chain_name}2{out_chain}{extension}"

    if output_file.exists():
        logger.info("Output file %s already exists for input file %s. Skipping.", output_file, input_file)
        return output_file

    if chain_name == out_chain and len(chainnames_in_structure) == 1:
        logger.info(
            "%s only has chain %s and out_chain is also %s. Copying file to %s.",
            input_file,
            chain_name,
            out_chain,
            output_file,
        )
        copyfile(input_file, output_file, copy_method)
        return output_file

    gemmi.Selection(f"/1/{chain_name}").remove_not_selected(structure)
    for model in structure:
        model.remove_ligands_and_waters()
    structure.setup_entities()
    structure.rename_chain(chain_name, out_chain)
    _dedup_helices(structure)
    _dedup_sheets(structure, out_chain)
    _add_provenance_info(structure, chain_name, out_chain)

    if not (len(structure) == 1 and len(structure[0]) == 1 and len(structure[0][out_chain]) > 0):
        msg = (
            f"After processing, structure does not have exactly one model ({len(structure)}) "
            f"with one chain (found {len(structure[0])}) called {out_chain} "
            f"with some residues ({len(structure[0][out_chain])})."
        )
        raise ValueError(msg)

    write_structure(structure, output_file)

    return output_file


def nr_of_residues_in_total(structure: Structure) -> int:
    """Count the total number of residues in the structure.

    Args:
        structure: The gemmi Structure object to analyze.

    Returns:
        The total number of residues in the structure."""
    count = 0
    for model in structure:
        for chain in model:
            count += len(chain)
    return count
