"""Chain-level structure helpers and transformations."""

import logging
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import gemmi
from gemmi import SoftwareItem, Structure

from protein_quest.__version__ import __version__
from protein_quest.converter import converter
from protein_quest.structure.errors import ChainNotFoundError
from protein_quest.structure.files import split_name_and_extension
from protein_quest.structure.formats import read_structure, write_structure
from protein_quest.utils import CopyMethod, copyfile

logger = logging.getLogger(__name__)

CHAIN_PROVENANCE_SOFTWARE_NAME = "protein-quest.structure.chains.write_single_chain_structure_file"
"""Name stored in structure metadata to record chain extraction provenance."""
ChainIdSystem = Literal["auth", "label"]
"""Which chain identifier system is used.

* ``label``: PDB-assigned chain id (``label_asym_id`` in mmcif).
* ``auth``: author-reported chain id (``auth_asym_id`` in mmcif).

If they differ, chain ids are shown as
``label_asym_id [auth auth_asym_id]`` on [https://www.rcsb.org/](https://www.rcsb.org/).
"""


def find_chain_in_model(model: gemmi.Model, wanted_chain: str) -> gemmi.Chain | None:
    """Find a chain in a model.

    Args:
        model: The gemmi model to search in.
        wanted_chain: The chain identifier to search for.
            Interpreted in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].

    Returns:
        The found chain or None if not found.
            Returned chain object is in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem]."""
    chain = model.find_chain(wanted_chain)
    if chain is None:
        mchains = [c for c in model if c.name.endswith(wanted_chain)]
        if mchains:
            return mchains[0]
    return chain


def get_label2auth_chains(structure: gemmi.Structure) -> dict[str, str]:
    """Build a label-to-author chain mapping from a structure.

    This function primarily reads mmCIF ``_atom_site.label_asym_id`` and
    ``_atom_site.auth_asym_id`` columns from ``group_PDB == 'ATOM'`` rows to derive
    ``label_asym_id -> auth_asym_id``.

    Args:
        structure: The structure to inspect.

    Returns:
        A dictionary mapping label chain ids to author chain ids.
            Keys are in 'label'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].
            Values are in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].
            If the same label appears multiple times, the first observed mapping
            is kept to ensure deterministic output.
    """
    # as atoms site is largest block we do not filter with MmcifOutputGroups
    block = structure.make_mmcif_block()
    atom_site = block.get_mmcif_category("_atom_site.")
    label_asym_ids = atom_site.get("label_asym_id", [])
    auth_asym_ids = atom_site.get("auth_asym_id", [])
    group_pdb_values = atom_site.get("group_PDB", [])

    # Empirical note: in a sample of ~14k structures, label/auth chain pairs were always 1:1.
    # so no many to one (B2A + C2A) or one to many (A2B + A2C) mappings were observed.
    # We therefore keep first-seen label->auth pairs in a dict and can safely invert when needed.
    label2auth: dict[str, str] = {}

    for label_asym_id, auth_asym_id, group_pdb in zip(label_asym_ids, auth_asym_ids, group_pdb_values, strict=False):
        if group_pdb != "ATOM":
            # Skip HETATM
            continue
        if label_asym_id not in label2auth:
            label2auth[label_asym_id] = auth_asym_id
    return label2auth


def find_chain_in_structure(
    structure: gemmi.Structure, wanted_chain: str, chain_system: ChainIdSystem = "auth"
) -> gemmi.Chain | None:
    """Find a chain in a structure.

    Args:
        structure: The gemmi structure to search in.
        wanted_chain: The chain identifier to search for.
        chain_system: System of ``wanted_chain``.

    Returns:
        The found chain or None if not found.
            Returned chain object is in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].
    """
    if chain_system == "label":
        label2auth = get_label2auth_chains(structure)
        try:
            wanted_chain = label2auth[wanted_chain]
        except KeyError:
            logger.warning("Label chain %s not found in structure. Returning None.", wanted_chain)
            return None
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
            Interpreted in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].

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


@dataclass(frozen=True)
class ChainExtractionProvenance:
    """Provenance information for chain extraction.

    Attributes:
        chain2keep: The chain identifier that was kept from the input structure.
            Stored in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].
        out_chain: The chain identifier that is used in this output structure.
            Stored in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem]."""

    chain2keep: str
    out_chain: str


def _add_provenance_info(structure: gemmi.Structure, chain2keep: str, out_chain: str):
    new_si = gemmi.SoftwareItem()
    new_si.classification = gemmi.SoftwareItem.Classification.DataExtraction
    new_si.name = CHAIN_PROVENANCE_SOFTWARE_NAME
    new_si.version = __version__
    new_si.date = str(datetime.now(tz=UTC).date())
    chain_provenance = converter.dumps(
        ChainExtractionProvenance(
            chain2keep=chain2keep,
            out_chain=out_chain,
        )
    ).decode("utf-8")
    new_si.contact_author = chain_provenance
    structure.meta.software = [*structure.meta.software, new_si]


def retrieve_chain_extraction_provenance(
    structure: gemmi.Structure,
) -> tuple[SoftwareItem, ChainExtractionProvenance] | None:
    """Extract the provenance information from a structure.

    Gives back what chain renamed as a result of call to `protein-quest filter chain` command or
    [write_single_chain_structure_file][protein_quest.structure.chains.write_single_chain_structure_file]
    function.

    Args:
        structure: The gemmi structure to extract provenance from.

    Returns:
        A tuple of the software item and the provenance information, or None if not found.

    Raises:
        JSONDecodeError: If the contact_author field is not valid JSON.
        ClassValidationError: If the contact_author field is valid JSON is incorrect shape.
    """
    for software_item in reversed(structure.meta.software):
        if software_item.name != CHAIN_PROVENANCE_SOFTWARE_NAME or not software_item.contact_author:
            continue
        contact_author = converter.loads(software_item.contact_author, ChainExtractionProvenance)
        return software_item, contact_author
    return None


def chains_in_structure(structure: gemmi.Structure) -> set[gemmi.Chain]:
    """Get a list of chains in a structure.

    Args:
        structure: The gemmi structure to get chains from.

    Returns:
        A set of chains in the structure.
            Returned chain objects are in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem]."""
    return {c for model in structure for c in model}


def _normalize_single_chain_entities(structure: gemmi.Structure, source_entity_id: str, out_chain: str):
    kept_entity_index = next(
        (index for index, entity in enumerate(structure.entities) if entity.name == source_entity_id),
        None,
    )
    if kept_entity_index is None:
        msg = f"Could not find entity {source_entity_id}."
        raise ValueError(msg)

    kept_entity = structure.entities[kept_entity_index]

    kept_entity.name = "1"
    kept_entity.subchains = [out_chain]
    for model in structure:
        for chain in model:
            for residue in chain:
                residue.subchain = chain.name
                residue.entity_id = "1"

    for entity_index in range(len(structure.entities) - 1, -1, -1):
        if entity_index != kept_entity_index:
            structure.entities.pop(entity_index)


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
    - stores chain2keep and out_chain as JSON-ified
        [ChainExtractionProvenance][protein_quest.structure.chains.ChainExtractionProvenance]
        object in the `contact_author` field of a new software item.
        The software item also contains this function name, version and current date.

    This function is equivalent to the following gemmi commands:

    ```shell
    gemmi convert --remove-lig-wat --select=B --to=cif chain-in/3JRS.cif - | \
    gemmi convert --from=cif --rename-chain=B:A - chain-out/3JRS_B2A.gemmi.cif
    ```

    Args:
        input_file: Path to the input structure file.
        chain2keep: The chain to keep.
            Interpreted in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].
        output_dir: Directory to save the output file.
        out_chain: The chain identifier for the output file.
            Written in 'auth'
            [chain id system][protein_quest.structure.chains.ChainIdSystem].
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
    source_entity_id = _extract_source_entity_id(structure, chain_name)
    structure.setup_entities()
    structure.rename_chain(chain_name, out_chain)
    _normalize_single_chain_entities(structure, source_entity_id, out_chain)
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


def _extract_source_entity_id(structure: Structure, chain_name: str) -> str:
    remaining_entity_ids = {
        residue.entity_id for model in structure for chain in model for residue in chain if residue.entity_id
    }
    if len(remaining_entity_ids) != 1:
        msg = f"Could not determine a unique entity for source chain {chain_name}."
        raise ValueError(msg)
    return next(iter(remaining_entity_ids))


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
