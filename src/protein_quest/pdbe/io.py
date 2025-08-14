import logging
from pathlib import Path

import gemmi

from protein_quest import __version__

logger = logging.getLogger(__name__)


def is_chain_in_residues_range(
    file: Path | str,
    min_residues: int,
    max_residues: int,
    chain: str = "A",
) -> bool:
    """Checks if a specific chain in a mmCIF/pdb file has a number of residues within a specified range.

    Args:
        file: Path to the input mmCIF/pdb file.
        min_residues: Minimum number of residues.
        max_residues: Maximum number of residues.
        chain: Chain to check.

    Returns:
        True if the chain has a number of residues within the specified range, False otherwise.
    """
    nr_residues = nr_residues_in_chain(file, chain)
    return min_residues <= nr_residues <= max_residues


def nr_residues_in_chain(file: Path | str, chain: str = "A") -> int:
    """Returns the number of residues in a specific chain from a mmCIF/pdb file.

    Args:
        file: Path to the input mmCIF/pdb file.
        chain: Chain to count residues of.

    Returns:
        The number of residues in the specified chain.
    """
    structure = gemmi.read_structure(str(file))
    model = structure[0]
    gchain = find_chain_in_model(model, chain)
    if gchain is None:
        return 0
    return len(gchain)


def find_chain_in_model(model: gemmi.Model, wanted_chain: str) -> gemmi.Chain | None:
    chain = model.find_chain(wanted_chain)
    if chain is None:
        # For chain A in 4v92 the find_chain method returns None,
        # however it is prefixed with 'B',
        # so we try again as last char of chain name
        mchains = [c for c in model if c.name.endswith(wanted_chain)]
        if mchains:
            return mchains[0]
    return chain


def write_single_chain_pdb_file(
    input_file: Path, chain2keep: str, output_dir: Path, out_chain: str = "A"
) -> Path | None:
    """Write a single chain PDB file from a mmCIF file.

    Args:
        input_file: Path to the input mmCIF file.
        chain2keep: The chain to keep.
        output_dir: Directory to save the output PDB file.
        out_chain: The chain identifier for the output PDB file.

    Returns:
        Path to the output PDB file or None if not created.
    """

    structure = gemmi.read_structure(str(input_file))
    model = structure[0]

    # Only count residues of polymer
    model.remove_ligands_and_waters()

    chain = find_chain_in_model(model, chain2keep)
    if chain is None:
        logger.warning(
            "Chain %s not found in %s. Skipping.",
            chain2keep,
            input_file,
        )
        return None
    stemmed_input_file = input_file.stem.replace(".gz", "").replace(".cif", "").replace(".pdb", "")
    output_file = output_dir / f"{stemmed_input_file}_{chain.name}2{out_chain}.pdb"

    new_structure = gemmi.Structure()
    new_structure.resolution = structure.resolution
    new_id = structure.name + f"{chain2keep}2{out_chain}"
    new_structure.name = new_id
    new_structure.info["_entry.id"] = new_id
    new_title = f"From {structure.info['_entry.id']} chain {chain2keep} to {out_chain}"
    new_structure.info["_struct.title"] = new_title
    new_structure.info["_struct_keywords.pdbx_keywords"] = new_title.upper()
    new_si = gemmi.SoftwareItem()
    new_si.classification = gemmi.SoftwareItem.Classification.DataExtraction
    new_si.name = "protein-quest"
    new_si.version = str(__version__)
    new_structure.meta.software.append(new_si)
    new_model = gemmi.Model(1)
    chain.name = out_chain
    new_model.add_chain(chain)
    new_structure.add_model(new_model)
    new_structure.write_pdb(str(output_file))

    return output_file
