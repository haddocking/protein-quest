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

from textwrap import dedent
from typing import Annotated

from fastmcp import FastMCP
from pydantic import Field

from protein_quest.uniprot import PdbResult, Query, search4pdb, search4uniprot

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


# TODO add all cli subcommands as mcp tools
