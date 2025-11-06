"""MCP server for protein-quest.

Can be run with:

```shell
# for development
fastmcp dev src/protein_quest/mcp_server.py
# or from inspector
npx @modelcontextprotocol/inspector
# tranport type: stdio
# comand: protein-quest
# arguments: mcp

# or with server and inspector
protein-quest mcp --transport streamable-http
# in another shell
npx @modelcontextprotocol/inspector
# transport type: streamable http
# URL: http://127.0.0.1:8000/mcp

# or with copilot in VS code
# ctrl + shift + p
# mcp: add server...
# Choose STDIO
# command: uv run protein-quest mcp
# id: protein-quest
```

Examples:

   - What are the PDBe structures for `A8MT69` uniprot accession?

"""

from collections.abc import Mapping
from pathlib import Path
from textwrap import dedent
from typing import Annotated

from fastmcp import FastMCP
from pydantic import Field

from protein_quest.alphafold.confidence import ConfidenceFilterQuery, ConfidenceFilterResult, filter_file_on_confidence
from protein_quest.alphafold.fetch import AlphaFoldEntry, DownloadableFormat
from protein_quest.alphafold.fetch import fetch_many as alphafold_fetch
from protein_quest.emdb import fetch as emdb_fetch
from protein_quest.go import search_gene_ontology_term
from protein_quest.io import convert_to_cif_file, glob_structure_files, read_structure
from protein_quest.pdbe.fetch import fetch as pdbe_fetch
from protein_quest.ss import filter_file_on_secondary_structure
from protein_quest.structure import (
    nr_residues_in_chain,
    structure2uniprot_accessions,
    write_single_chain_structure_file,
)
from protein_quest.taxonomy import search_taxon
from protein_quest.uniprot import (
    PdbResult,
    Query,
    search4af,
    search4emdb,
    search4macromolecular_complexes,
    search4pdb,
    search4uniprot,
)

mcp = FastMCP("protein-quest")

# do not want to make dataclasses in non-mcp code into Pydantic models,
# so we use Annotated here to add description on roots.


@mcp.tool
def search_uniprot(
    uniprot_query: Annotated[Query, Field(description=Query.__doc__)],
    limit: Annotated[int, Field(gt=0, description="Limit the number of uniprot accessions returned")] = 100,
) -> set[str]:
    """Search UniProt for proteins matching the given query."""
    return search4uniprot(uniprot_query, limit=limit)


@mcp.tool
def search_pdb(
    uniprot_accs: set[str],
    limit: Annotated[int, Field(gt=0, description="Limit the number of entries returned")] = 100,
) -> Annotated[
    dict[str, set[PdbResult]],
    Field(
        description=dedent(f"""\
            Dictionary with protein IDs as keys and sets of PDB results as values.
            A PDB result is {PdbResult.__doc__}""")
    ),
]:
    """Search PDBe structures for given uniprot accessions."""
    return search4pdb(uniprot_accs, limit=limit)


@mcp.tool
async def fetch_pdbe_structures(pdb_ids: set[str], save_dir: Path) -> Mapping[str, Path]:
    """Fetch the PDBe structures for given PDB IDs.

    Args:
        pdb_ids: A set of PDB IDs.
        save_dir: The directory to save the fetched files.

    Returns:
        A mapping of PDB ID to the path of the fetched structure file.
    """
    return await pdbe_fetch(pdb_ids, save_dir)


@mcp.tool
def extract_single_chain_from_structure(
    input_file: Path,
    chain2keep: str,
    output_dir: Path,
    out_chain: str = "A",
) -> Path:
    """
    Extract a single chain from a structure (mmCIF or pdb) file and write to a new file.

    Args:
        input_file: Path to the input structure (mmCIF or pdb) file.
        chain2keep: The chain to keep.
        output_dir: Directory to save the output file.
        out_chain: The chain identifier for the output file.

    Returns:
        Path to the output structure (mmCIF or pdb) file
    """
    return write_single_chain_structure_file(input_file, chain2keep, output_dir, out_chain)


@mcp.tool
def list_structure_files(path: Path) -> list[Path]:
    """List structure files (.pdb, .pdb.gz, .cif, .cif.gz, .bcif) in the specified directory."""
    return list(glob_structure_files(path))


# TODO replace remaining decorators with wrapper if tool does single function call
# so we do not have to replicate docstring,
# minor con is that it does not show up in api docs
mcp.tool(nr_residues_in_chain)
mcp.tool(search_taxon)
mcp.tool(search_gene_ontology_term)


@mcp.tool
def search_alphafolds(
    uniprot_accs: set[str],
    limit: Annotated[int, Field(gt=0, description="Limit the number of entries returned")] = 100,
) -> Annotated[
    set[str],
    Field(description="Set of uniprot accessions which have an AlphaFold entry"),
]:
    """Search for AlphaFold entries in UniProtKB accessions."""
    # each uniprot accession can have one or more AlphaFold IDs
    # an AlphaFold ID is the same as the uniprot accession
    # so we return a subset of uniprot_accs
    results = search4af(uniprot_accs, limit)
    return {k for k, v in results.items() if v}


mcp.tool(search4emdb, name="search_emdb")
mcp.tool(search4macromolecular_complexes, name="search_macromolecular_complexes")


@mcp.tool
def fetch_alphafold_structures(uniprot_accs: set[str], save_dir: Path) -> list[AlphaFoldEntry]:
    """Fetch the AlphaFold mmCIF file for given UniProt accessions.

    Args:
        uniprot_accs: A set of UniProt accessions.
        save_dir: The directory to save the fetched files.

    Returns:
        A list of AlphaFold entries.
    """
    formats: set[DownloadableFormat] = {"cif"}
    return alphafold_fetch(uniprot_accs, save_dir, formats)


@mcp.tool
async def fetch_emdb_volumes(emdb_ids: set[str], save_dir: Path) -> Mapping[str, Path]:
    """Fetch EMDB volumes for given EMDB IDs.

    Args:
        emdb_ids: A set of EMDB IDs.
        save_dir: The directory to save the fetched files.
    Returns:
        A mapping of EMDB ID to the path of the fetched volume file.
    """
    return await emdb_fetch(emdb_ids=emdb_ids, save_dir=save_dir)


@mcp.tool
def alphafold_confidence_filter(file: Path, query: ConfidenceFilterQuery, filtered_dir: Path) -> ConfidenceFilterResult:
    """Take a mmcif/PDB file and filter it based on confidence (plDDT) scores.

    If passes filter writes file to filtered_dir with residues above confidence threshold.
    """
    return filter_file_on_confidence(file, query, filtered_dir)


mcp.tool(filter_file_on_secondary_structure)

mcp.tool(convert_to_cif_file)


@mcp.tool
def uniprot_accessions_of_structure_file(file: Path) -> set[str]:
    """Extract UniProt accessions from structure file."""
    structure = read_structure(file)
    return structure2uniprot_accessions(structure)


@mcp.prompt
def candidate_structures(
    species: str = "Human",
    cellular_location: str = "nucleus",
    confidence: int = 90,
    min_residues: int = 100,
    max_residues: int = 200,
) -> str:
    """Prompt to find candidate structures.

    Args:
        species: The species to search for (default: "Human").
        cellular_location: The cellular location to search for (default: "nucleus").
        confidence: The confidence threshold for AlphaFold structures (default: 90).
        min_residues: Minimum number of high confidence residues (default: 100).
        max_residues: Maximum number of high confidence residues (default: 200).

    Returns:
        A prompt string to find candidate structures.
    """
    return dedent(f"""\
        Given the species '{species}' and cellular location '{cellular_location}' find the candidate structures.
        Download structures from 2 sources namely PDB and Alphafold.
        For alphafold I only want to use high confidence scores of over {confidence}.
        and only keep structures with number of high confidence residues between {min_residues} and {max_residues}.

        1. Search uniprot for proteins related to {species} and {cellular_location}.
            1. For the species find the NCBI taxonomy id.
            2. For cellular location find the associated GO term.
            3. Find uniprot accessions based on NCBI taxonomy id and cellular location GO term.
        2. For PDB
            1. Search for structures related to the identified proteins.
            2. Download each PDB entry from PDBe
            3. Extract chain for the protein of interest.
        3. For Alphafold
            1. Search for AlphaFold entries related to the identified proteins.
            2. Download each AlphaFold entry.
            3. Filter the structures based on {confidence} as confidence
               and nr residues between {min_residues} and {max_residues}.
    """)
