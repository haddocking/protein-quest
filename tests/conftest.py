from pathlib import Path

import gemmi
import pytest

from protein_quest.structure.formats import read_structure, write_structure


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


@pytest.fixture
def em_cif() -> Path:
    """EM structure for P0ABE7 (8w77, resolution in weird place)."""
    return Path(__file__).parent / "fixtures" / "8w77_updated.cif.gz"


@pytest.fixture
def sample_multispan_cif() -> Path:
    """6O5I structure with multiple `_struct_ref_seq` spans for one chain."""
    return Path(__file__).parent / "fixtures" / "6O5I.cif.gz"


@pytest.fixture
def multi_accession_cif() -> Path:
    """1A02 structure with multiple UniProt accessions in separate chains."""
    return Path(__file__).parent / "fixtures" / "1a02.cif.gz"


@pytest.fixture
def multi_accession_chain_cif() -> Path:
    """1UN5 structure with multiple UniProt accessions in the same chain."""
    return Path(__file__).parent / "fixtures" / "1un5.cif.gz"


@pytest.fixture
def no_uniprot_cif(sample2_cif: Path, tmp_path: Path) -> Path:
    """2Y29 structure with no UniProt accession."""
    structure_with_unp = read_structure(sample2_cif)
    block_without_struct_ref = structure_with_unp.make_mmcif_block(
        gemmi.MmcifOutputGroups(True, chem_comp=False, struct_ref=False)
    )
    structure = gemmi.make_structure_from_block(block_without_struct_ref)
    write_structure(structure, tmp_path / "no_uniprot.cif")
    return tmp_path / "no_uniprot.cif"


@pytest.fixture
def all_cifs(
    sample_cif: Path,
    sample2_cif: Path,
    sample_multispan_cif: Path,
    multi_accession_chain_cif: Path,
    af_cif: Path,
    nmr_cif: Path,
    em_cif: Path,
) -> list[Path]:
    """List of all CIF fixtures except multi_accession_cif and no_uniprot_cif as they raise an error."""
    return [
        sample_cif,
        sample2_cif,
        sample_multispan_cif,
        multi_accession_chain_cif,
        af_cif,
        nmr_cif,
        em_cif,
    ]
