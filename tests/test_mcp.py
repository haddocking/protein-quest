from pathlib import Path

import pytest
from fastmcp import Client

from protein_quest.mcp_server import mcp


@pytest.fixture
def cif_path() -> Path:
    return Path(__file__).parent / "pdbe" / "fixtures" / "2y29.cif"


@pytest.mark.asyncio
async def test_nr_residues_in_chain(cif_path: Path):
    async with Client(mcp) as client:
        result = await client.call_tool("nr_residues_in_chain", {"file": cif_path})
        assert result.data == 8
