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
# Prompt: What are the PDBe structures for `A8MT69` uniprot accession?
```

Examples:

    For search pdb use `A8MT69` as input.

"""

from collections.abc import Mapping
from pathlib import Path
from textwrap import dedent
from typing import Annotated

from fastmcp import FastMCP
from pydantic import Field

from protein_quest.alphafold.confidence import ConfidenceFilterQuery, ConfidenceFilterResult, filter_file_on_residues
from protein_quest.alphafold.fetch import AlphaFoldEntry, DownloadableFormat, fetch_many
from protein_quest.emdb import fetch as emdb_fetch
from protein_quest.go import search_gene_ontology_term
from protein_quest.pdbe.fetch import fetch as pdbe_fetch
from protein_quest.pdbe.io import glob_structure_files, nr_residues_in_chain, write_single_chain_pdb_file
from protein_quest.taxonomy import search_taxon
from protein_quest.uniprot import PdbResult, Query, search4af, search4emdb, search4pdb, search4uniprot

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
def fetch_structures_from_pdbe(
    ids: set[str], save_dir: Path
) -> Annotated[Mapping[str, Path], Field(description="Mapping of PDB IDs to their file paths.")]:
    """Fetch PDB structures as mmCIF files from PDBe and save them to the specified directory."""
    return pdbe_fetch(ids, save_dir)


@mcp.tool
def extract_single_chain_from_structure(
    input_file: Path,
    chain2keep: str,
    output_dir: Path,
    out_chain: str = "A",
) -> Path | None:
    """
    Extract a single chain from a mmCIF/pdb file and write to a new file.

    Args:
        input_file: Path to the input mmCIF/pdb file.
        chain2keep: The chain to keep.
        output_dir: Directory to save the output file.
        out_chain: The chain identifier for the output file.

    Returns:
        Path to the output mmCIF/pdb file or None if not created.
    """
    return write_single_chain_pdb_file(input_file, chain2keep, output_dir, out_chain)


@mcp.tool
def list_structure_files(path: Path) -> list[Path]:
    """List structure files (.pdb, .pdb.gz, .cif, .cif.gz) in the specified directory."""
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
    # each uniprot accesion can have one or more AlphaFold IDs
    # an AlphaFold ID is the same as the uniprot accession
    # so we return a subset of uniprot_accs
    results = search4af(uniprot_accs, limit)
    return {k for k, v in results.items() if v}


mcp.tool(search4emdb, name="search_emdb")


@mcp.tool
def fetch_alphafold_structures(uniprot_accs: set[str], save_dir: Path) -> list[AlphaFoldEntry]:
    """Fetch the AlphaFold summary and mmcif file for given UniProt accessions.

    Args:
        uniprot_accs: A set of UniProt accessions.
        save_dir: The directory to save the fetched files.

    Returns:
        A list of AlphaFold entries.
    """
    what: set[DownloadableFormat] = {"cif"}
    return fetch_many(uniprot_accs, save_dir, what)


mcp.tool(emdb_fetch, name="fetch_emdb_volumes")


@mcp.tool
def alphafold_confidence_filter(file: Path, query: ConfidenceFilterQuery, filtered_dir: Path) -> ConfidenceFilterResult:
    """Take a mmcif/PDB file and filter it based on confidence (plDDT) scores.

    If passes filter writes file to filtered_dir with residues above confidence threshold.
    """
    return filter_file_on_residues(file, query, filtered_dir)


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
