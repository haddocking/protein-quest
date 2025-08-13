from pathlib import Path

import gemmi
import pytest

from protein_quest.pdbe.io import first_chain_from_uniprot_chains, write_single_chain_pdb_file


@pytest.mark.parametrize(
    "query,expected",
    [
        ("O=1-300", "O"),  #  uniprot:A8MT69 pdb:7R5S
        ("B/D=1-81", "B"),  # uniprot:A8MT69 pdb:4E44
        (
            "B/D/H/L/M/N/U/V/W/X/Z/b/d/h/i/j/o/p/q/r=8-81",  # uniprot:A8MT69 pdb:4NE1
            "B",
        ),
        ("A/B=2-459,A/B=520-610", "A"),  # uniprot/O00255 pdb/3U84
        ("DD/Dd=1-1085", "DD"),  # uniprot/O00268 pdb/7ENA
        ("A=398-459,A=74-386,A=520-584,A=1-53", "A"),  # uniprot/O00255 pdb/7O9T
    ],
)
def test_first_chain_from_uniprot_chains(query, expected):
    result = first_chain_from_uniprot_chains(query)

    assert result == expected


@pytest.fixture
def cif_path() -> Path:
    return Path(__file__).parent / "fixtures" / "2y29.cif"


def test_write_single_chain_pdb_file_happypath(cif_path: Path, tmp_path: Path):
    chain2keep = "A=1-6"

    output_file = write_single_chain_pdb_file(
        input_file=cif_path,
        uniprot_chain=chain2keep,
        uniprot_acc="SOMEACC",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is not None
    assert output_file.exists()
    structure = gemmi.read_structure(str(output_file))
    assert len(structure) == 1  # One model
    model = structure[0]
    assert len(model) == 1  # One chain
    chain = model[0]
    assert chain.name == "Z"
    assert len(chain) == 6  # 6 residues in chain Z


def test_write_single_chain_pdb_file_unknown_chain(cif_path: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture):
    chain2keep = "B=1-20"

    output_file = write_single_chain_pdb_file(
        cif_path,
        uniprot_chain=chain2keep,
        uniprot_acc="SOMEACC",
        output_dir=tmp_path,
        out_chain="Z",
    )

    assert output_file is None
    assert "Chain B not found in" in caplog.text
