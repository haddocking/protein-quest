from pathlib import Path

import gemmi
import pytest

from protein_quest.structure.formats import read_structure
from protein_quest.structure.types import Pdb2UniprotMapping, StructRefSeq
from protein_quest.structure.uniprot import (
    add_uniprot_accessions2structure,
    struct_ref_seqs_columns_to_records,
    structure2uniprot_accessions,
    structure_to_uniprot,
)


def test_structure2uniprot_accessions_present(sample2_cif: Path):
    structure = read_structure(sample2_cif)
    accessions = structure2uniprot_accessions(structure)

    assert accessions == {"P05067"}


def test_structure2uniprot_accessions_multiple(multi_accession_cif: Path):
    structure = read_structure(multi_accession_cif)
    accessions = structure2uniprot_accessions(structure)

    assert accessions == {"Q13469", "P01100", "P05412"}


def test_structure2uniprot_accessions_missing(sample_cif: Path, caplog: pytest.LogCaptureFixture):
    # Empty struct_ref category to simulate missing UniProt accessions
    structure_with_unp = read_structure(sample_cif)
    block_without_struct_ref = structure_with_unp.make_mmcif_block(
        gemmi.MmcifOutputGroups(True, chem_comp=False, struct_ref=False)
    )
    structure = gemmi.make_structure_from_block(block_without_struct_ref)

    accessions = structure2uniprot_accessions(structure)

    assert accessions == set()
    assert "No UniProt accessions found in structure 3JRSB2A" in caplog.text


@pytest.mark.parametrize(
    ("struct_ref_seqs_columns", "expected_records"),
    [
        pytest.param(
            {
                "align_id": [1],
                "ref_id": [1],
                "pdbx_PDB_id_code": ["1ABC"],
                "pdbx_strand_id": ["A"],
                "seq_align_beg": [10],
                "pdbx_seq_align_beg_ins_code": ["?"],
                "seq_align_end": [15],
                "pdbx_seq_align_end_ins_code": ["?"],
                "pdbx_db_accession": ["P12345"],
                "db_align_beg": [10],
                "pdbx_db_align_beg_ins_code": ["?"],
                "db_align_end": [15],
                "pdbx_db_align_end_ins_code": ["?"],
                "pdbx_auth_seq_align_beg": [10],
                "pdbx_auth_seq_align_end": [15],
            },
            [StructRefSeq("P12345", 10, 15, "A", 1.0, 6)],
            id="single-span",
        ),
        pytest.param(
            {
                "align_id": [1, 2],
                "ref_id": [1, 2],
                "pdbx_PDB_id_code": ["1ABC", "1ABC"],
                "pdbx_strand_id": ["A", "A"],
                "seq_align_beg": [10, 20],
                "pdbx_seq_align_beg_ins_code": ["?", "?"],
                "seq_align_end": [15, 24],
                "pdbx_seq_align_end_ins_code": ["?", "?"],
                "pdbx_db_accession": ["P12345", "P12345"],
                "db_align_beg": [10, 20],
                "pdbx_db_align_beg_ins_code": ["?", "?"],
                "db_align_end": [15, 24],
                "pdbx_db_align_end_ins_code": ["?", "?"],
                "pdbx_auth_seq_align_beg": [10, 20],
                "pdbx_auth_seq_align_end": [15, 24],
            },
            [StructRefSeq("P12345", 10, 24, "A", 11 / 15, 11)],
            id="multi-span-single-chain",
        ),
        pytest.param(
            {
                "align_id": [1, 2],
                "ref_id": [1, 2],
                "pdbx_PDB_id_code": ["1ABC", "1ABC"],
                "pdbx_strand_id": ["A", "B"],
                "seq_align_beg": [10, 30],
                "pdbx_seq_align_beg_ins_code": ["?", "?"],
                "seq_align_end": [15, 35],
                "pdbx_seq_align_end_ins_code": ["?", "?"],
                "pdbx_db_accession": ["P12345", "P12345"],
                "db_align_beg": [10, 30],
                "pdbx_db_align_beg_ins_code": ["?", "?"],
                "db_align_end": [15, 35],
                "pdbx_db_align_end_ins_code": ["?", "?"],
                "pdbx_auth_seq_align_beg": [10, 30],
                "pdbx_auth_seq_align_end": [15, 35],
            },
            [
                StructRefSeq("P12345", 10, 15, "A", 1.0, 6),
                StructRefSeq("P12345", 30, 35, "B", 1.0, 6),
            ],
            id="multi-chain",
        ),
        pytest.param(
            {
                "align_id": [1, 2, 3],
                "ref_id": [1, 1, 2],
                "pdbx_PDB_id_code": ["1ABC", "1ABC", "1ABC"],
                "pdbx_strand_id": ["A", "A", "A"],
                "seq_align_beg": [10, 20, 100],
                "pdbx_seq_align_beg_ins_code": ["?", "?", "?"],
                "seq_align_end": [12, 22, 105],
                "pdbx_seq_align_end_ins_code": ["?", "?", "?"],
                "pdbx_db_accession": ["PAAAAA", "PAAAAA", "PBBBBB"],
                "db_align_beg": [10, 20, 100],
                "pdbx_db_align_beg_ins_code": ["?", "?", "?"],
                "db_align_end": [12, 22, 105],
                "pdbx_db_align_end_ins_code": ["?", "?", "?"],
                "pdbx_auth_seq_align_beg": [10, 20, 100],
                "pdbx_auth_seq_align_end": [12, 22, 105],
            },
            [
                StructRefSeq("PAAAAA", 10, 22, "A", 6 / 13, 6),
                StructRefSeq("PBBBBB", 100, 105, "A", 1.0, 6),
            ],
            id="equal-aligned-residues-same-chain",
        ),
        pytest.param(
            {
                "align_id": [1],
                "ref_id": [1],
                "pdbx_PDB_id_code": ["1ABC"],
                "pdbx_strand_id": ["A"],
                "seq_align_beg": [42],
                "pdbx_seq_align_beg_ins_code": ["?"],
                "seq_align_end": [42],
                "pdbx_seq_align_end_ins_code": ["?"],
                "pdbx_db_accession": ["P11111"],
                "db_align_beg": [42],
                "pdbx_db_align_beg_ins_code": ["?"],
                "db_align_end": [42],
                "pdbx_db_align_end_ins_code": ["?"],
                "pdbx_auth_seq_align_beg": [42],
                "pdbx_auth_seq_align_end": [42],
            },
            [StructRefSeq("P11111", 42, 42, "A", 1.0, 1)],
            id="single-residue-inclusive-count",
        ),
    ],
)
def test_struct_ref_seqs_columns_to_records(
    struct_ref_seqs_columns: dict[str, list[int | str]], expected_records: list[StructRefSeq]
):
    records = struct_ref_seqs_columns_to_records(struct_ref_seqs_columns)

    assert sorted(records, key=lambda r: (r.uniprot_accession, r.chain_id)) == sorted(
        expected_records,
        key=lambda r: (r.uniprot_accession, r.chain_id),
    )


class TestStructureToUniprot:
    def test_from_sifts(self, nmr_cif: Path):
        structure = read_structure(nmr_cif)

        result = structure_to_uniprot(structure)

        expected: Pdb2UniprotMapping = {
            "1AMB": {
                ("A", "P05067"),
            }
        }
        assert result == expected

    def test_from_struct_ref(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)

        result = structure_to_uniprot(structure)

        expected: Pdb2UniprotMapping = {"2Y29": {("A", "P05067")}}
        assert result == expected

class TestVerifyInjectUniprotRef:
    def test_none(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)

        new_structure = add_uniprot_accessions2structure(structure, None)
        assert structure == new_structure

    def test_missing_id(self, sample2_cif: Path, caplog: pytest.LogCaptureFixture):
        structure = read_structure(sample2_cif)
        pdb2uniprot = {"1AAA": {("A", "P12345")}}  # wrong PDB ID
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        assert structure == new_structure
        assert (
            "PDB ID 2Y29 not found in pdb2uniprot mapping. Leaving structure unverified and unchanged." in caplog.text
        )

    def test_verify_ok(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)
        pdb2uniprot = {"2Y29": {("A", "P05067")}}
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert structure == new_structure

    def test_inject_into_nostructref(self, no_uniprot_cif: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = {"2Y29": {("A", "P12345")}}
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        result2 = structure_to_uniprot(new_structure)

        expected: Pdb2UniprotMapping = {"2Y29": {("A", "P12345")}}
        assert result2 == expected

    def test_inject_into_existing_sifts(self, nmr_cif: Path):
        structure = read_structure(nmr_cif)
        pdb2uniprot = {"1AMB": {("A", "P12345")}}
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        result2 = structure_to_uniprot(new_structure)

        expected: Pdb2UniprotMapping = {"1AMB": {("A", "P05067"), ("A", "P12345")}}
        assert result2 == expected
