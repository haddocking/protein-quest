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

from collections.abc import Collection, Generator, Mapping
from pathlib import Path
from textwrap import dedent
from typing import Annotated

from fastmcp import FastMCP
from pydantic import Field

from protein_quest.go import Aspect, search_go_term
from protein_quest.pdbe.fetch import fetch as pdbe_fetch
from protein_quest.pdbe.io import glob_structure_files, nr_residues_in_chain, write_single_chain_pdb_file
from protein_quest.uniprot import PdbResult, Query, Taxon, search4pdb, search4taxon, search4uniprot

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
    limit: Annotated[int, Field(gt=0, description="Limit the number of pdb entries returned")] = 100,
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
def fetch_structure_from_pdbe(
    ids: set[str], save_dir: Path
) -> Annotated[Mapping[str, Path], Field(description="Mapping of PDB IDs to their file paths.")]:
    """Fetch PDB structures as mmCIF files from PDBe and save them to the specified directory."""
    return pdbe_fetch(ids, save_dir)


# MCP tool for extracting a single chain from a structure file
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
def list_structure_files(path: Path) -> Generator[Path]:
    """List structure files (.pdb, .pdb.gz, .cif, .cif.gz) in the specified directory."""
    yield from glob_structure_files(path)


@mcp.tool
def count_residues_in_chain(file: Path, chain: str = "A") -> int:
    """Count the number of residues in a specific chain of a structure file."""
    return nr_residues_in_chain(file, chain)


@mcp.tool
def search_taxon_by_name(term: str) -> Collection[Taxon]:
    """Search NCBI Taxonomy by common or scientific name."""
    return search4taxon(term)


@mcp.tool
def search_gene_ontology_term(term: str, aspect: Aspect | None = None):
    """Search Gene Ontology (GO) terms by name and aspect.

    If aspect is not provided, all aspects are included.
    """
    return search_go_term(term, aspect)

# TODO add all cli subcommands as mcp tools
# - Alphafold fetch and filter
# - use @mcp.resource and @mcp.prompt
