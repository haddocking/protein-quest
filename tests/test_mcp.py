from pathlib import Path

import pytest
from fastmcp import Client

from protein_quest.mcp_server import mcp


@pytest.mark.asyncio
async def test_nr_residues_in_chain(sample_cif: Path):
    async with Client(mcp) as client:
        result = await client.call_tool("nr_residues_in_chain", {"file": sample_cif})
        assert result.data == 173
