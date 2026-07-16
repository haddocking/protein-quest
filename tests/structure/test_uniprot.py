from pathlib import Path

import gemmi
import pytest

from protein_quest.structure.chains import ChainIdSystem, write_single_chain_structure_file
from protein_quest.structure.formats import read_structure
from protein_quest.structure.types import Pdb2UniprotMapping, StructRefSeq
from protein_quest.structure.uniprot import (
    UniprotSource,
    add_uniprot_accessions2structure,
    selected_struct_ref_seqs_by_chain,
    selected_struct_ref_seqs_from_sifts_by_chain,
    struct_ref_seqs_columns_to_records,
    structure2uniprot_accessions,
    structure_to_uniprot,
)
from protein_quest.uniprot_chains import Pdb2UniprotChainsMapping, UniprotChainMapping, parse_uniprot_chains


def _mapping(pdb_id: str, uniprot_accession: str, uniprot_chains: str) -> Pdb2UniprotChainsMapping:
    return {
        pdb_id: {
            UniprotChainMapping(
                uniprot_accession=uniprot_accession,
                chain_ranges=parse_uniprot_chains(uniprot_chains),
            )
        }
    }


@pytest.mark.parametrize(
    ("fixture_name", "chain2keep", "expected"),
    [
        pytest.param(
            "cif_3jrs",
            "B",
            {"3JRS": {("A", "Q8VZS8")}},
            id="cif_3jrs",
        ),
        pytest.param(
            "multi_entity_cif",
            "B",
            {"1F66": {("A", "P62806")}},
            id="multi_entity_cif",
        ),
    ],
)
def test_write_single_chain_structure_file_preserves_uniprot_mapping(
    request: pytest.FixtureRequest, tmp_path: Path, fixture_name: str, chain2keep: str, expected: Pdb2UniprotMapping
):
    input_file = request.getfixturevalue(fixture_name)
    output_file = write_single_chain_structure_file(input_file=input_file, chain2keep=chain2keep, output_dir=tmp_path)
    structure = read_structure(output_file)

    result = structure_to_uniprot(structure)

    assert result == expected


def test_structure2uniprot_accessions_present(sample2_cif: Path):
    structure = read_structure(sample2_cif)
    accessions = structure2uniprot_accessions(structure)

    assert accessions == {"P05067"}


def test_structure2uniprot_accessions_multiple(multi_accession_cif: Path):
    structure = read_structure(multi_accession_cif)
    accessions = structure2uniprot_accessions(structure)

    assert accessions == {"Q13469", "P01100", "P05412"}


def test_structure2uniprot_accessions_missing(sample_cif: Path):
    # Empty struct_ref category to simulate missing UniProt accessions
    structure_with_unp = read_structure(sample_cif)
    block_without_struct_ref = structure_with_unp.make_mmcif_block(
        gemmi.MmcifOutputGroups(True, chem_comp=False, struct_ref=False)
    )
    structure = gemmi.make_structure_from_block(block_without_struct_ref)

    accessions = structure2uniprot_accessions(structure)

    assert accessions == set()


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
    @pytest.mark.parametrize(
        ("fixture_name", "source", "expected"),
        [
            pytest.param(
                "nmr_cif",
                "sifts",
                {
                    "1AMB": {
                        ("A", "P05067"),
                    }
                },
                id="sifts-only",
            ),
            pytest.param(
                "sample2_cif",
                "struct_ref_seq",
                {"2Y29": {("A", "P05067")}},
                id="struct-ref-seq-only",
            ),
            pytest.param(
                "multi_entity_cif",
                "both",
                {
                    "1F66": {
                        (
                            "A",
                            "P84233",
                        ),
                        (
                            "A",
                            "Q7ZT64",
                        ),
                        (
                            "B",
                            "P62806",
                        ),
                        (
                            "C",
                            "P0C0S5",
                        ),
                        (
                            "C",
                            "P17317",
                        ),
                        (
                            "D",
                            "P02281",
                        ),
                        (
                            "E",
                            "P84233",
                        ),
                        (
                            "E",
                            "Q7ZT64",
                        ),
                        (
                            "F",
                            "P62806",
                        ),
                        (
                            "G",
                            "P0C0S5",
                        ),
                        (
                            "G",
                            "P17317",
                        ),
                        (
                            "H",
                            "P02281",
                        ),
                    }
                },
                id="multi-both-explicit",
            ),
            pytest.param(
                "multi_entity_cif",
                "sifts",
                {
                    "1F66": {
                        (
                            "A",
                            "P84233",
                        ),
                        (
                            "B",
                            "P62806",
                        ),
                        (
                            "C",
                            "P0C0S5",
                        ),
                        (
                            "D",
                            "P02281",
                        ),
                        (
                            "E",
                            "P84233",
                        ),
                        (
                            "F",
                            "P62806",
                        ),
                        (
                            "G",
                            "P0C0S5",
                        ),
                        (
                            "H",
                            "P02281",
                        ),
                    }
                },
                id="multi-sifts",
            ),
            pytest.param(
                "multi_entity_cif",
                "struct_ref_seq",
                {
                    "1F66": {
                        (
                            "A",
                            "Q7ZT64",
                        ),
                        (
                            "B",
                            "P62806",
                        ),
                        (
                            "C",
                            "P17317",
                        ),
                        (
                            "D",
                            "P02281",
                        ),
                        (
                            "E",
                            "Q7ZT64",
                        ),
                        (
                            "F",
                            "P62806",
                        ),
                        (
                            "G",
                            "P17317",
                        ),
                        (
                            "H",
                            "P02281",
                        ),
                    }
                },
                id="multi-struct_ref_seq",
            ),
            pytest.param(
                "multi_entity_cif",
                "fallback",
                {
                    "1F66": {
                        (
                            "A",
                            "P84233",
                        ),
                        (
                            "B",
                            "P62806",
                        ),
                        (
                            "C",
                            "P0C0S5",
                        ),
                        (
                            "D",
                            "P02281",
                        ),
                        (
                            "E",
                            "P84233",
                        ),
                        (
                            "F",
                            "P62806",
                        ),
                        (
                            "G",
                            "P0C0S5",
                        ),
                        (
                            "H",
                            "P02281",
                        ),
                    }
                },
                id="multi-fallback-prefers-sifts",
            ),
            pytest.param(
                "sample2_cif",
                "fallback",
                {"2Y29": {("A", "P05067")}},
                id="fallback-uses-struct-ref-seq-when-sifts-empty",
            ),
            pytest.param(
                "multi_entity_cif",
                None,
                {
                    "1F66": {
                        (
                            "A",
                            "P84233",
                        ),
                        (
                            "B",
                            "P62806",
                        ),
                        (
                            "C",
                            "P0C0S5",
                        ),
                        (
                            "D",
                            "P02281",
                        ),
                        (
                            "E",
                            "P84233",
                        ),
                        (
                            "F",
                            "P62806",
                        ),
                        (
                            "G",
                            "P0C0S5",
                        ),
                        (
                            "H",
                            "P02281",
                        ),
                    }
                },
                id="multi-default",
            ),
            pytest.param(
                "cif_2fui",
                "sifts",
                {"2FUI": {("A", "Q12830")}},
                id="just sift",
            ),
        ],
    )
    def test_sources(
        self,
        request: pytest.FixtureRequest,
        fixture_name: str,
        source: UniprotSource | None,
        expected: Pdb2UniprotMapping,
    ):
        structure_path = request.getfixturevalue(fixture_name)
        structure = read_structure(structure_path)

        result = structure_to_uniprot(structure) if source is None else structure_to_uniprot(structure, source=source)

        assert result == expected


class TestAddUniprotAccessions2Structure:
    def test_none(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)

        new_structure = add_uniprot_accessions2structure(structure, None)
        assert structure == new_structure

    def test_missing_id(self, sample2_cif: Path, caplog: pytest.LogCaptureFixture):
        structure = read_structure(sample2_cif)
        pdb2uniprot = _mapping("1AAA", "P12345", "A=1-10")  # wrong PDB ID
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        assert structure == new_structure
        assert (
            "PDB ID 2Y29 not found in pdb2uniprot mapping. Leaving structure unverified and unchanged." in caplog.text
        )

    def test_verify_ok(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)
        pdb2uniprot = _mapping("2Y29", "P05067", "A=1-770")
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert structure == new_structure

    def test_inject_into_nostructref(self, no_uniprot_cif: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = _mapping("2Y29", "P12345", "A=10-20")
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        result2 = structure_to_uniprot(new_structure)

        expected: Pdb2UniprotMapping = {"2Y29": {("A", "P12345")}}
        assert result2 == expected

        block = new_structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
        struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
        injected_records = [
            record
            for record in struct_ref_seqs_columns_to_records(struct_ref_seq)
            if record.uniprot_accession == "P12345"
        ]
        assert injected_records == [StructRefSeq("P12345", 10, 20, "A", 1.0, 11)]

    def test_inject_into_existing_sifts(self, nmr_cif: Path):
        structure = read_structure(nmr_cif)
        pdb2uniprot = _mapping("1AMB", "P12345", "A=1-100")
        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        result2 = structure_to_uniprot(new_structure)

        expected: Pdb2UniprotMapping = {"1AMB": {("A", "P05067"), ("A", "P12345")}}
        assert result2 == expected

    def test_inject_multiple_ranges_into_nostructref(self, no_uniprot_cif: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = _mapping("2Y29", "P12345", "A=10-20,A=30-35")

        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)

        block = new_structure.make_mmcif_block(gemmi.MmcifOutputGroups(False, struct_ref=True))
        struct_ref = block.get_mmcif_category("_struct_ref.")
        # The uniprot accession is injected twice, once for each range,
        # would have _struct_ref.pdbx_align_begin=[10,30]
        # but gemmi does not preserve that column, so we just check the other columns are correct
        expected_struct_ref = {
            "id": ["1", "2"],
            "entity_id": ["1", "1"],
            "db_name": ["UNP", "UNP"],
            "db_code": [False, False],
            "pdbx_db_accession": ["P12345", "P12345"],
            "pdbx_db_isoform": [None, None],
        }
        assert struct_ref == expected_struct_ref
        struct_ref_seq = block.get_mmcif_category("_struct_ref_seq.")
        injected_records = [
            record
            for record in struct_ref_seqs_columns_to_records(struct_ref_seq)
            if record.uniprot_accession == "P12345"
        ]
        assert injected_records == [StructRefSeq("P12345", 10, 35, "A", 17 / 26, 17)]

    def test_on_single_chain_written(self, multi_accession_cif: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        caplog.set_level("INFO")
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_file = write_single_chain_structure_file(
            input_file=multi_accession_cif, output_dir=input_dir, chain2keep="F"
        )

        structure = read_structure(input_file)
        pdb2uniprot = _mapping("1A02", "P01111", "F=1-100")

        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)
        result = structure_to_uniprot(new_structure)

        expected: Pdb2UniprotMapping = {"1A02": {("A", "P01100"), ("A", "P01111")}}
        assert result == expected

        log = caplog.text
        assert "Structure 1A02 has provenance information indicating it was extracted from chain F to A" in log
        assert "Using this information to verify/add UniProt accessions." in log

    @pytest.mark.parametrize(
        ("chain_system", "chain_id"),
        [
            ("auth", "C"),
            # label chain E is same entity as auth chain C
            ("label", "E"),
        ],
    )
    def test_multi_entity_cif(self, multi_entity_cif: Path, chain_system: ChainIdSystem, chain_id: str):
        structure = read_structure(multi_entity_cif)
        pdb2uniprot = _mapping("1F66", "P12345", f"{chain_id}=1-100")

        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot, chain_system=chain_system)

        result = structure_to_uniprot(new_structure, source="both")

        expected: Pdb2UniprotMapping = {
            "1F66": {
                (
                    "A",
                    "Q7ZT64",
                ),
                (
                    "B",
                    "P62806",
                ),
                (
                    "C",
                    # as requested
                    "P12345",
                ),
                (
                    "C",
                    "P17317",
                ),
                (
                    "D",
                    "P02281",
                ),
                (
                    "E",
                    "Q7ZT64",
                ),
                (
                    "F",
                    "P62806",
                ),
                (
                    "G",
                    # Also updates P0C0S5 to requested as G and C auth chains are same entity
                    "P12345",
                ),
                (
                    "G",
                    "P17317",
                ),
                (
                    "H",
                    "P02281",
                ),
            },
        }
        assert result == expected

    def test_multi_entity_cif_multiple_auth_chains_in_single_mapping(self, multi_entity_cif: Path):
        structure = read_structure(multi_entity_cif)
        pdb2uniprot = _mapping("1F66", "P12345", "A/B=1-100")

        new_structure = add_uniprot_accessions2structure(structure, pdb2uniprot)
        result = structure_to_uniprot(new_structure, source="both")

        assert ("A", "P12345") in result["1F66"]
        assert ("B", "P12345") in result["1F66"]


def test_selected_struct_ref_seqs_by_chain_returns_auth_system(cif_8rw8: Path):
    structure = read_structure(cif_8rw8)
    accessions = structure2uniprot_accessions(structure)

    result = selected_struct_ref_seqs_by_chain(structure, accessions)

    # cif_8rw8 fixture has auth chain B and label chain A, we expect auth chain B
    expected: dict[str, StructRefSeq] = {
        "B": StructRefSeq(
            uniprot_accession="O00327",
            uniprot_start=337,
            uniprot_end=449,
            chain_id="B",
            sequence_identity=1.0,
            aligned_residue_count=113,
        )
    }
    assert result == expected


@pytest.mark.parametrize(
    "cif_fixture, expected",
    [
        pytest.param(
            "cif_2fui",
            {
                "A": StructRefSeq(
                    uniprot_accession="Q12830",
                    uniprot_start=2865,
                    uniprot_end=2921,
                    chain_id="A",
                    sequence_identity=0.3903,
                    aligned_residue_count=1140,
                )
            },
            id="single-chain",
        ),
        pytest.param(
            "sample2_cif",
            {},
            id="struct-ref-seq-only",
        ),
        pytest.param(
            "multi_entity_cif",
            {
                "A": StructRefSeq(
                    uniprot_accession="P84233",
                    uniprot_start=37,
                    uniprot_end=136,
                    chain_id="A",
                    sequence_identity=0.7353,
                    aligned_residue_count=100,
                ),
                "B": StructRefSeq(
                    uniprot_accession="P62806",
                    uniprot_start=24,
                    uniprot_end=103,
                    chain_id="B",
                    sequence_identity=0.7767,
                    aligned_residue_count=80,
                ),
                "C": StructRefSeq(
                    uniprot_accession="P0C0S5",
                    uniprot_start=17,
                    uniprot_end=119,
                    chain_id="C",
                    sequence_identity=0.8655,
                    aligned_residue_count=103,
                ),
                "D": StructRefSeq(
                    uniprot_accession="P02281",
                    uniprot_start=32,
                    uniprot_end=126,
                    chain_id="D",
                    sequence_identity=0.754,
                    aligned_residue_count=95,
                ),
                "E": StructRefSeq(
                    uniprot_accession="P84233",
                    uniprot_start=34,
                    uniprot_end=136,
                    chain_id="E",
                    sequence_identity=0.7574,
                    aligned_residue_count=103,
                ),
                "F": StructRefSeq(
                    uniprot_accession="P62806",
                    uniprot_start=18,
                    uniprot_end=103,
                    chain_id="F",
                    sequence_identity=0.835,
                    aligned_residue_count=86,
                ),
                "G": StructRefSeq(
                    uniprot_accession="P0C0S5",
                    uniprot_start=17,
                    uniprot_end=123,
                    chain_id="G",
                    sequence_identity=0.8699,
                    aligned_residue_count=107,
                ),
                "H": StructRefSeq(
                    uniprot_accession="P02281",
                    uniprot_start=32,
                    uniprot_end=126,
                    chain_id="H",
                    sequence_identity=0.754,
                    aligned_residue_count=95,
                ),
            },
            id="multi-entity",
        ),
        pytest.param(
            "multi_accession_chain_cif",
            {},
            id="multi-accession-chain",
        ),
    ],
)
def test_selected_struct_ref_seqs_from_sifts_by_chain(
    cif_fixture: str, expected: dict[str, StructRefSeq], request: pytest.FixtureRequest
):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)
    accessions = structure2uniprot_accessions(structure)

    actual = selected_struct_ref_seqs_from_sifts_by_chain(structure, accessions)
    assert actual == expected
