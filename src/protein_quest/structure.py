"""Module for querying and modifying [gemmi structures][gemmi.Structure]."""

import logging
from datetime import UTC, datetime
from pathlib import Path

import gemmi

from protein_quest.__version__ import __version__
from protein_quest.io import read_structure, split_name_and_extension, write_structure
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)


def find_chain_in_model(model: gemmi.Model, wanted_chain: str) -> gemmi.Chain | None:
    """Find a chain in a model.

    Args:
        model: The gemmi model to search in.
        wanted_chain: The chain identifier to search for.

    Returns:
        The found chain or None if not found.
    """
    chain = model.find_chain(wanted_chain)
    if chain is None:
        # For chain A in 4v92 the find_chain method returns None,
        # however it is prefixed with 'B',
        # so we try again as last char of chain name
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
        The found chain or None if not found.
    """
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
        The number of residues in the specified chain.
    """
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
            logger.debug(f"Duplicate start helix found: {hindex} {helix.start}, removing")
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
    new_si.version = str(__version__)
    new_si.date = str(datetime.now(tz=UTC).date())
    structure.meta.software = [*structure.meta.software, new_si]


def chains_in_structure(structure: gemmi.Structure) -> set[gemmi.Chain]:
    """Get a list of chains in a structure.

    Args:
        structure: The gemmi structure to get chains from.

    Returns:
        A set of chains in the structure.
    """
    return {c for model in structure for c in model}


class ChainNotFoundError(IndexError):
    """Exception raised when a chain is not found in a structure."""

    def __init__(self, chain_id: str, file: Path | str, available_chains: set[str]):
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
    gemmi convert --remove-lig-wat --select=B --to=cif chain-in/3JRS.cif - | \\
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
        ChainNotFoundError: If the specified chain is not found in the input file.
    """

    logger.debug(f"chain2keep: {chain2keep}, out_chain: {out_chain}")
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
    for m in structure:
        m.remove_ligands_and_waters()
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


def structure2uniprot_accessions(structure: gemmi.Structure) -> set[str]:
    """Extract UniProt accessions from a gemmi Structure object.

    Logs a warning and returns an empty set if no accessions are found in structure.

    Args:
        structure: The gemmi Structure object to extract UniProt accessions from.

    Returns:
        A set of UniProt accessions found in the structure.
    """
    block = structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
    struct_ref = block.get_mmcif_category("_struct_ref.")
    uniprot_accessions: set[str] = set()
    for i, db_name in enumerate(struct_ref["db_name"]):
        if db_name != "UNP":
            continue
        pdbx_db_accession = struct_ref["pdbx_db_accession"][i]
        uniprot_accessions.add(pdbx_db_accession)
    if not uniprot_accessions:
        logger.warning("No UniProt accessions found in structure %s", structure.name)
    return uniprot_accessions
