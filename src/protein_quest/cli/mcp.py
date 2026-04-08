"""MCP subcommand for protein-quest."""

import sys
from typing import Literal

from cyclopts import App

from protein_quest.cli.common import Common, console

# Copy from from fastmcp.server.server import Transport
# Duplicated here so cli works without fastmcpp dep installed
Transport = Literal["stdio", "http", "sse", "streamable-http"]

rprint = console.print
mcp_app = App(name="mcp", help="Run Model Context Protocol (MCP) server")


@mcp_app.default
def mcp(
    *,
    transport: Transport | None = "stdio",
    host: str = "127.0.0.1",
    port: int = 8000,
    _: Common | None = None,
) -> None:
    """Run Model Context Protocol (MCP) server.

    Run Model Context Protocol (MCP) server.
    Can be used by agentic LLMs like Claude Sonnet 4 as a set of tools.

    Args:
        transport: Transport protocol to use.
        host: Host to bind the server to.
        port: Port to bind the server to.
        _: Common CLI options.
    """
    try:
        from protein_quest.mcp_server import mcp as mcp_server  # noqa: PLC0415

        if transport == "stdio":
            mcp_server.run(transport=transport)
        else:
            mcp_server.run(transport=transport, host=host, port=port)
    except ImportError:
        rprint(
            "[red]Error:[/red] The `protein-quest mcp` command requires the 'fastmcp' Python package to be installed."
        )
        sys.exit(1)
