from pathlib import Path

import gemmi
import pooch
import pytest
from platformdirs import user_cache_dir

from protein_quest.structure.formats import read_structure, write_structure


@pytest.fixture
def sample_cif() -> Path:
    """Downloaded from https://www.ebi.ac.uk/pdbe/entry-files/download/3jrs_updated.cif.gz
    and filtered with
    `write_single_chain_structure_file(Path('tests/fixtures/3jrs_updated.cif.gz'), 'B', Path('tests/fixtures/'))`
    """
    return Path(__file__).parent / "fixtures" / "3jrs_updated_B2A.cif.gz"


@pytest.fixture
def sample2_cif() -> Path:
    """X-ray structure downloaded from https://www.rcsb.org/structure/2Y29"""
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
def unreadable_cif(tmp_path: Path) -> Path:
    """Create a CIF file that cannot be read by gemmi."""
    unreadable_cif_path = tmp_path / "unreadable.cif"
    unreadable_cif_path.write_text("This is not a valid CIF file.")
    return unreadable_cif_path


@pytest.fixture
def atomless_cif(tmp_path: Path) -> Path:
    structure = gemmi.Structure()
    structure.name = "1ABC"
    structure.info["_entry.id"] = "1ABC"
    fn = tmp_path / "atomless.cif"
    write_structure(structure, fn)
    # Gemmi has shortcut for atomless structures, which causes name to be  " " when read.
    return fn


@pytest.fixture
def download_cache_dir() -> Path:
    cache_dir = Path(user_cache_dir("protein-quest-tests"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def fetch_cif(filename: str, sha256: str) -> Path:
    base_url = "https://www.ebi.ac.uk/pdbe/entry-files/download/"
    path = pooch.retrieve(
        url=f"{base_url}{filename}",
        known_hash=f"sha256:{sha256}",
        fname=filename,
        # keep path same as download_cache_dir fixture, for ci caching
        path=pooch.os_cache("protein-quest-tests"),
    )
    return Path(path)


@pytest.fixture
def multi_entity_cif() -> Path:
    """1f66 structure with multiple entities."""
    return fetch_cif(
        "1f66_updated.cif.gz",
        "208386eac478c4deb9767f66f6ade97dc869a8d5395cdadfb385c9597442a56e",
    )


@pytest.fixture
def cif_8rw8() -> Path:
    """8rw8 structure with a mismatch between label (A) and auth (B) chain."""
    return fetch_cif(
        "8rw8_updated.cif.gz",
        "882e27541671578935db963787a720c3c942199049f9b4cfa99edbac5239a946",
    )


@pytest.fixture
def cif_3jrs() -> Path:
    """3jrs structure which is used to create sample_cif fixture."""
    return fetch_cif(
        "3jrs_updated.cif.gz",
        "6117a1ef3d5d655491367588c56747fe9a4c5132dd240401f03f6ad3645d7603",
    )


@pytest.fixture
def cif_2y2a() -> Path:
    """2y2a x-ray structure with same uniprot accession as sample2_cif fixture aka 2Y29."""
    return fetch_cif(
        "2y2a_updated.cif.gz",
        "eaec5ba1d0b744fc7561e18ecbe1951e875913455450c102f3707519a616d092",
    )


@pytest.fixture
def cif_2fui():
    """2fui x-ray structure with uniprot from just sift."""
    return fetch_cif(
        "2fui_updated.cif.gz",
        "f70adacd1ec1b8bde59c3754a6774cf190a9d3dc0710b6041bab9ad96c7118f5",
    )


@pytest.fixture
def xray_p05067_cifs(sample2_cif: Path, cif_2y2a: Path) -> list[Path]:
    return [sample2_cif, cif_2y2a]


@pytest.fixture
def fake_archive_em_structure() -> gemmi.Structure:
    structure = gemmi.Structure()
    structure.name = "1ABC"
    # Archived structures use UPPERCASE as exp method value
    structure.info["_exptl.method"] = "ELECTRON MICROSCOPY"
    structure.resolution = 4.2
    atom = gemmi.Atom()
    atom.name = "CA"
    atom.element = gemmi.Element("C")
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.add_atom(atom)
    residue.entity_type = gemmi.EntityType.Polymer
    chain = gemmi.Chain("A")
    chain.add_residue(residue)
    model = gemmi.Model(1)
    model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    structure.assign_subchains()
    return structure


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


@pytest.fixture
def bad_cifs(
    no_uniprot_cif: Path,
    unreadable_cif: Path,
    atomless_cif: Path,
) -> list[Path]:
    """List of all CIF fixtures that raise an error when read or when metadata is extracted."""
    return [
        no_uniprot_cif,
        unreadable_cif,
        atomless_cif,
    ]
