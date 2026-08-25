import gzip
from pathlib import Path

import gemmi
import pytest
from gemmi import cif

from protein_quest.structure.formats import (
    gunzip_file,
    read_structure,
    read_structure_as_cif_block,
    structure2cifgz,
    write_structure,
)
from protein_quest.structure.sifts import uniprot_chain_mappings_from_cif
from protein_quest.structure.types import valid_structure_file_extensions
from protein_quest.structure.uniprot_extraction import structure_to_uniprot


def test_invalid_extension(sample2_cif: Path, tmp_path: Path):
    structure = read_structure(sample2_cif)
    output_file = tmp_path / "bla.txt"

    with pytest.raises(ValueError, match="Unsupported file extension"):
        write_structure(structure, output_file)


@pytest.mark.parametrize("extension", valid_structure_file_extensions)
def test_read_structure(sample2_cif: Path, tmp_path: Path, extension: str):
    # We only have cif as fixture, so convert to other formats first
    structure_from_cif = read_structure(sample2_cif)
    thefile = tmp_path / f"foo{extension}"
    write_structure(structure_from_cif, thefile)

    structure_from_extension = read_structure(thefile)

    assert structure_from_extension.make_pdb_string() == structure_from_cif.make_pdb_string()


def test_read_structure_as_cif_block(sample2_cif: Path):
    block = read_structure_as_cif_block(sample2_cif)

    assert block is not None
    assert block.name == "2Y29"


def test_read_structure_as_cif_block_no_blocks(comments_only_cif: Path):
    block = read_structure_as_cif_block(comments_only_cif)

    assert block is None


@pytest.mark.parametrize("extension", [".bcif", ".bcif.gz"])
def test_read_structure_as_cif_block_binary_formats(sample2_cif: Path, tmp_path: Path, extension: str):
    structure = read_structure(sample2_cif)
    thefile = tmp_path / f"foo{extension}"
    write_structure(structure, thefile)

    block = read_structure_as_cif_block(thefile)

    assert block is not None
    assert block.name == "2Y29"


def test_read_structure_as_cif_block_logs_read_error(af_pdb: Path, caplog: pytest.LogCaptureFixture):
    caplog.set_level("INFO")

    block = read_structure_as_cif_block(af_pdb)

    assert block is None
    assert "Unable to read mmCIF block from" in caplog.text
    assert str(af_pdb) in caplog.text
    assert "expected block header (data_)" in caplog.text


@pytest.mark.parametrize("extension", [".bcif", ".bcif.gz"])
def test_read_structure_as_cif_block_logs_binary_read_error(
    tmp_path: Path, extension: str, caplog: pytest.LogCaptureFixture
):
    caplog.set_level("INFO")
    bad_file = tmp_path / f"bad{extension}"
    if extension == ".bcif":
        bad_file.write_bytes(b"not a bcif")
    else:
        bad_file.write_bytes(gzip.compress(b"not a bcif"))

    block = read_structure_as_cif_block(bad_file)

    assert block is None
    assert "Unable to read mmCIF block from" in caplog.text
    assert str(bad_file) in caplog.text


def test_structure2cifgz(sample2_cif: Path):
    structure = read_structure(sample2_cif)

    cif_gz_bytes = structure2cifgz(structure)

    assert isinstance(cif_gz_bytes, bytes)
    assert len(cif_gz_bytes) > 0
    observed_cif_text = gzip.decompress(cif_gz_bytes).decode("utf-8")
    assert observed_cif_text.startswith("data_")


def test_gunzip_of_non_gz_file(tmp_path: Path):
    input_file = tmp_path / "bla.zip"
    input_file.write_text("some text")

    with pytest.raises(ValueError, match="must end with \\.gz"):
        gunzip_file(input_file)


class TestWriteStructure:
    @pytest.mark.parametrize("extension", valid_structure_file_extensions)
    def test_extensions(self, sample2_cif: Path, tmp_path: Path, extension: str):
        structure = read_structure(sample2_cif)
        output_file = tmp_path / f"bla{extension}"

        write_structure(structure, output_file)

        assert output_file.exists()
        assert output_file.stat().st_size > 0

    @pytest.mark.parametrize(
        ["output_fn"],
        [
            ("em_structure.cif",),
            ("em_structure.cif.gz",),
        ],
    )
    def test_em_structure_retains_resolution(self, em_cif: Path, tmp_path: Path, output_fn: str):
        structure = read_structure(em_cif)
        assert structure.resolution == 3.61
        output_file = tmp_path / output_fn

        write_structure(structure, output_file)

        written_structure = read_structure(output_file)

        assert written_structure.resolution == 3.61

    def test_em_archived_structure_retains_resolution(self, fake_archive_em_structure: gemmi.Structure, tmp_path: Path):
        structure = fake_archive_em_structure
        assert structure.resolution == 4.2
        output_file = tmp_path / "em_structure.cif"

        write_structure(structure, output_file)

        written_structure = read_structure(output_file)

        assert written_structure.resolution == 4.2

    def test_keeps_sifts_roundtrip(self, cif_3plz: Path, tmp_path: Path):
        structure = read_structure(cif_3plz)
        expected = structure_to_uniprot(cif_3plz, structure=structure, source="sifts", one_uniprot_per_chain=False)
        assert len(expected) == 4
        output_file = tmp_path / "3plz_roundtrip.cif.gz"
        block = read_structure_as_cif_block(cif_3plz)

        assert block is not None

        write_structure(structure, output_file, uniprot_chain_mappings=uniprot_chain_mappings_from_cif(block))

        written_structure = read_structure(output_file)
        result = structure_to_uniprot(
            output_file, structure=written_structure, source="sifts", one_uniprot_per_chain=False
        )
        assert len(result) == 4

        assert result == expected

    def test_keeps_sifts_with_empty_polymer_chains(self, cif_1l5w: Path, tmp_path: Path):
        structure = read_structure(cif_1l5w)
        expected = structure_to_uniprot(cif_1l5w, structure=structure, source="sifts", one_uniprot_per_chain=False)
        assert len(expected) == 2
        output_file = tmp_path / "1l5w_roundtrip.cif.gz"
        block = read_structure_as_cif_block(cif_1l5w)

        assert block is not None

        write_structure(structure, output_file, uniprot_chain_mappings=uniprot_chain_mappings_from_cif(block))

        written_structure = read_structure(output_file)
        result = structure_to_uniprot(
            output_file, structure=written_structure, source="sifts", one_uniprot_per_chain=False
        )

        assert result == expected

    def test_omits_sifts_block_without_sifts(self, sample2_cif: Path, tmp_path: Path):
        structure = read_structure(sample2_cif)
        assert (
            structure_to_uniprot(sample2_cif, structure=structure, source="sifts", one_uniprot_per_chain=False) == set()
        )
        output_file = tmp_path / "2y29_roundtrip.cif.gz"

        write_structure(structure, output_file)

        observed_cif_text = gzip.decompress(output_file.read_bytes()).decode("utf-8")

        assert "_pdbx_sifts_xref_db." not in observed_cif_text

    def test_retains_chem_comp_id_not_type(self, cif_1l5w: Path, tmp_path: Path):
        original_block = cif.read_file(str(cif_1l5w)).sole_block()
        structure = read_structure(cif_1l5w)
        output_file = tmp_path / "1l5w_roundtrip.cif.gz"

        write_structure(structure, output_file)

        written_block = cif.read_file(str(output_file)).sole_block()
        expected_block: dict[str, list[str] | list[bool]] = {
            "id": list(original_block.find_values("_chem_comp.id")),
        }
        expected_block["type"] = [False] * len(expected_block["id"])
        assert written_block.get_mmcif_category("_chem_comp") == expected_block
