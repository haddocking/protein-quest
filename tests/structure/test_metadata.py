from pathlib import Path

import gemmi
import pytest

from protein_quest.structure.formats import read_structure
from protein_quest.structure.metadata import StructureMetadata, structure_metadata


class TestStructureMetadata:
    @pytest.mark.parametrize(
        "cif_fixture, expected",
        [
            pytest.param(
                "sample_cif",
                StructureMetadata(
                    id="3JRS",
                    uniprot_accession="Q8VZS8",
                    resolution=2.05,
                    total_residue_count=173,
                    is_alphafold=False,
                    uniprot_start=8,
                    uniprot_end=211,
                    sequence_identity=1.0,
                    chain_length=173,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                id="3JRS_B2A",
            ),
            pytest.param(
                "sample2_cif",
                StructureMetadata(
                    id="2Y29",
                    uniprot_accession="P05067",
                    resolution=2.3,
                    total_residue_count=8,
                    is_alphafold=False,
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=8,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                id="2Y29",
            ),
            pytest.param(
                "af_cif",
                StructureMetadata(
                    id="AF-A0A0C5B5G6-F1",
                    uniprot_accession="A0A0C5B5G6",
                    resolution=0.0,
                    total_residue_count=16,
                    is_alphafold=True,
                    uniprot_start=1,
                    uniprot_end=16,
                    sequence_identity=1.0,
                    chain_length=16,
                    auth_chain="A",
                    label_chain="A",
                    method="Predicted",
                ),
                id="AF-A0A0C5B5G6-F1",
            ),
            pytest.param(
                "nmr_cif",
                StructureMetadata(
                    id="1AMB",
                    uniprot_accession="P05067",
                    resolution=0.0,
                    total_residue_count=28,
                    is_alphafold=False,
                    uniprot_start=672,
                    uniprot_end=699,
                    sequence_identity=1.0,
                    chain_length=28,
                    auth_chain="A",
                    label_chain="A",
                    method="NMR",
                ),
                id="1AMB",
            ),
            pytest.param(
                "em_cif",
                StructureMetadata(
                    id="8W77",
                    uniprot_accession="P0ABE7",
                    resolution=3.61,
                    total_residue_count=260,
                    is_alphafold=False,
                    uniprot_start=23,
                    uniprot_end=127,
                    sequence_identity=1.0,
                    chain_length=260,
                    auth_chain="A",
                    label_chain="A",
                    method="EM",
                ),
                id="8W77",
            ),
            pytest.param(
                "sample_multispan_cif",
                StructureMetadata(
                    id="6O5I",
                    uniprot_accession="O00255",
                    resolution=1.240,
                    total_residue_count=1346,
                    is_alphafold=False,
                    uniprot_start=1,
                    uniprot_end=593,
                    sequence_identity=0.816,
                    chain_length=1346,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                id="6O5I_multispan",
            ),
            pytest.param(
                "multi_accession_cif",
                StructureMetadata(
                    id="1A02",
                    uniprot_accession=None,
                    resolution=2.7,
                    total_residue_count=513,
                    is_alphafold=False,
                    uniprot_start=0,
                    uniprot_end=0,
                    sequence_identity=0.0,
                    chain_length=40,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                id="1A02_multi_accessions-separate_chains",
            ),
            pytest.param(
                "multi_accession_chain_cif",
                StructureMetadata(
                    id="1UN5",
                    uniprot_accession="P03950",
                    resolution=2.600,
                    total_residue_count=131,
                    is_alphafold=False,
                    uniprot_start=25,
                    uniprot_end=147,
                    sequence_identity=119 / 123,
                    chain_length=131,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                id="1UN5_multi_accession_chain",
            ),
            pytest.param(
                "cif_8rw8",
                StructureMetadata(
                    id="8RW8",
                    uniprot_accession="O00327",
                    resolution=2.16,
                    total_residue_count=116,
                    is_alphafold=False,
                    uniprot_start=337,
                    uniprot_end=449,
                    sequence_identity=1.0,
                    chain_length=116,
                    auth_chain="B",
                    label_chain="A",
                    method="X-ray",
                ),
                id="label-auth-chain-mismatch",
            ),
            pytest.param(
                "no_uniprot_cif",
                StructureMetadata(
                    id="2Y29",
                    uniprot_accession=None,
                    resolution=2.3,
                    total_residue_count=8,
                    is_alphafold=False,
                    uniprot_start=0,
                    uniprot_end=0,
                    sequence_identity=0.0,
                    chain_length=8,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                id="no_uniprot_accession",
            ),
            pytest.param(
                "cif_2fui",
                StructureMetadata(
                    id="2FUI",
                    uniprot_accession="Q12830",
                    resolution=0.0,
                    total_residue_count=1240,
                    is_alphafold=False,
                    uniprot_start=2865,
                    uniprot_end=2921,
                    sequence_identity=1.0,
                    chain_length=62,
                    auth_chain="A",
                    label_chain="A",
                    method="NMR",
                ),
                id="sift-wins",
            ),
        ],
    )
    def test_cif_fixtures(self, cif_fixture: str, expected: StructureMetadata, request: pytest.FixtureRequest):
        path = request.getfixturevalue(cif_fixture)
        result = structure_metadata(read_structure(path), path=path)

        assert result == expected

    def test_multiple_accessions_warns(self, multi_accession_cif: Path, caplog: pytest.LogCaptureFixture):
        structure_metadata(read_structure(multi_accession_cif), path=multi_accession_cif)

        message = caplog.text
        assert "Multiple UniProt mappings found in structure" in message
        assert "Q13469" in message
        assert "P01100" in message
        assert "P05412" in message

    def test_empty_structure(self):
        structure = gemmi.Structure()
        structure.name = "EMPTY"

        with pytest.raises(ValueError, match="No chains found in structure EMPTY"):
            structure_metadata(structure)

    def test_archived_em_structure(self, fake_archive_em_structure: gemmi.Structure):
        result = structure_metadata(fake_archive_em_structure)

        expected = StructureMetadata(
            id="1ABC",
            uniprot_accession=None,
            resolution=4.2,
            total_residue_count=1,
            is_alphafold=False,
            uniprot_start=0,
            uniprot_end=0,
            sequence_identity=0.0,
            chain_length=1,
            auth_chain="A",
            label_chain="Axp",
            method="EM",
        )
        assert result == expected
