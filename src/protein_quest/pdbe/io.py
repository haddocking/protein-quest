import gzip
import logging
from collections.abc import Generator
from pathlib import Path

import gemmi

from protein_quest import __version__

logger = logging.getLogger(__name__)


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
        logger.warning("Chain %s not found in %s. Returning 0.", chain, file)
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


def write_structure(structure: gemmi.Structure, path: Path):
    """Write a gemmi structure to a file.

    Args:
        structure: The gemmi structure to write.
        path: The file path to write the structure to.
            The format depends on the file extension.
            Supported extensions are .pdb, .pdb.gz, .cif, .cif.gz.

    Raises:
        ValueError: If the file extension is not supported.
    """
    if path.name.endswith(".pdb"):
        body: str = structure.make_pdb_string()
        path.write_text(body)
    elif path.name.endswith(".pdb.gz"):
        body: str = structure.make_pdb_string()
        with gzip.open(path, "wt") as f:
            f.write(body)
    elif path.name.endswith(".cif"):
        doc = structure.make_mmcif_document()
        doc.write_file(str(path))
    elif path.name.endswith(".cif.gz"):
        doc = structure.make_mmcif_document()
        cif_str = doc.as_string()
        with gzip.open(path, "wt") as f:
            f.write(cif_str)
    else:
        msg = f"Unsupported file extension in {path.name}. Supported extensions are .pdb, .pdb.gz, .cif, .cif.gz"
        raise ValueError(msg)


def _split_name_and_extension(name: str) -> tuple[str, str]:
    # 1234.pdb -> (1234, .pdb)
    # 1234.pdb.gz -> (1234, .pdb.gz)
    # 1234.cif -> (1234, .cif)
    # 1234.cif.gz -> (1234, .cif.gz)
    if name.endswith(".pdb.gz"):
        return name.replace(".pdb.gz", ""), ".pdb.gz"
    if name.endswith(".cif.gz"):
        return name.replace(".cif.gz", ""), ".cif.gz"
    if name.endswith(".pdb"):
        return name.replace(".pdb", ""), ".pdb"
    if name.endswith(".cif"):
        return name.replace(".cif", ""), ".cif"

    msg = f"Unknown file extension in {name}. Supported extensions are .pdb, .pdb.gz, .cif, .cif.gz"
    raise ValueError(msg)


def locate_structure_file(root: Path, pdb_id: str) -> Path:
    """Locate a structure file for a given PDB ID in the specified directory.

    Args:
        root: The root directory to search in.
        pdb_id: The PDB ID to locate.

    Returns:
        The path to the located structure file.

    Raises:
        FileNotFoundError: If no structure file is found for the given PDB ID.
    """
    exts = [".cif.gz", ".cif", ".pdb.gz", ".pdb"]
    # files downloaded from https://www.ebi.ac.uk/pdbe/ website
    # have file names like pdb6t5y.ent or pdb6t5y.ent.gz for a PDB formatted file.
    # TODO support pdb6t5y.ent or pdb6t5y.ent.gz file names
    for ext in exts:
        candidate = root / f"{pdb_id.lower()}{ext}"
        if candidate.exists():
            return candidate
    msg = f"No structure file found for {pdb_id} in {root}"
    raise FileNotFoundError(msg)


def glob_structure_files(input_dir: Path) -> Generator[Path]:
    """Glob for structure files in a directory.

    Args:
        input_dir: The input directory to search for structure files.

    Yields:
        Paths to the found structure files.
    """
    for ext in [".cif.gz", ".cif", ".pdb.gz", ".pdb"]:
        yield from input_dir.glob(f"*{ext}")


def write_single_chain_pdb_file(
    input_file: Path, chain2keep: str, output_dir: Path, out_chain: str = "A"
) -> Path | None:
    """Write a single chain from a mmCIF/pdb file to a new mmCIF/pdb file.

    Args:
        input_file: Path to the input mmCIF/pdb file.
        chain2keep: The chain to keep.
        output_dir: Directory to save the output file.
        out_chain: The chain identifier for the output file.

    Returns:
        Path to the output mmCIF/pdb file or None if not created.
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
    name, extension = _split_name_and_extension(input_file.name)
    output_file = output_dir / f"{name}_{chain.name}2{out_chain}{extension}"

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
    write_structure(new_structure, output_file)

    return output_file
