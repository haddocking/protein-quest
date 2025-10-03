from pathlib import Path

import pytest


@pytest.fixture
def sample_cif() -> Path:
    """Downloaded from https://www.rcsb.org/structure/3JRS
    and filtered with
    `write_single_chain_pdb_file(Path('tests/fixtures/3JRS.cif.gz'), 'B', Path('tests/fixtures/'))`
    """
    return Path(__file__).parent / "fixtures" / "3JRS_B2A.cif.gz"


@pytest.fixture
def sample2_cif() -> Path:
    """Downloaded from https://www.rcsb.org/structure/2Y29"""
    return Path(__file__).parent / "fixtures" / "2Y29.cif.gz"
