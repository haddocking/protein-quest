from pathlib import Path

import pytest


@pytest.fixture
def sample_cif() -> Path:
    """Downloaded from https://www.rcsb.org/structure/3JRS
    and filtered with
    `write_single_chain_structure_file(Path('tests/fixtures/3JRS.cif.gz'), 'B', Path('tests/fixtures/'))`
    """
    return Path(__file__).parent / "fixtures" / "3JRS_B2A.cif.gz"


@pytest.fixture
def sample2_cif() -> Path:
    """Downloaded from https://www.rcsb.org/structure/2Y29"""
    return Path(__file__).parent / "fixtures" / "2Y29.cif.gz"


@pytest.fixture
def af_cif() -> Path:
    """AlphaFold structure for A0A0C5B5G6."""
    return Path(__file__).parent / "fixtures" / "AF-A0A0C5B5G6-F1-model_v6.cif.gz"


@pytest.fixture
def nmr_cif() -> Path:
    """NMR structure for P05067 (1AMB, no resolution)."""
    return Path(__file__).parent / "fixtures" / "1amb_updated.cif.gz"
