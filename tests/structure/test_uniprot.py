import gzip
import logging
from pathlib import Path

import pytest

from protein_quest.structure.chains import ChainExtractionProvenance, ChainIdSystem, write_single_chain_structure_file
from protein_quest.structure.formats import read_structure
from protein_quest.structure.uniprot import (
    FlattenedUniprotChainMapping,
    UniprotSource,
    add_uniprot_accessions2structure,
    apply_chain_provenance_to_uniprot_mappings,
    best_uniprot_per_chain,
    flatten_uniprot_chain_mappings,
    structure2uniprot_accessions,
    structure_to_uniprot,
    uniprot_chain_mappings_from_sifts,
    uniprot_chain_mappings_from_struct_ref_seq,
)
from protein_quest.uniprot_chains import (
    Pdb2UniprotChainsMapping,
    UniprotChainMapping,
    UniprotChainMappings,
    parse_uniprot_chains,
)


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
    "cif_fixture, expected",
    [
        pytest.param("atomless_cif", set(), id="atomless"),
        pytest.param("no_uniprot_cif", set(), id="no-uniprot"),
        pytest.param("sample2_cif", {"P05067"}, id="siftless"),
        pytest.param(
            "cif_2fui",
            {
                "Q12830",  # from sift
                "Q7Z7D6",  # from struct_ref_seq
            },
            id="1chain",
        ),
        pytest.param(
            "multi_accession_cif",
            {"Q13469", "P01100", "P05412"},
            id="multi-accession-separate-chains",
        ),
        pytest.param(
            "multi_accession_chain_cif",
            {"P03950", "P00656"},
            id="multi-accession-same-chain",
        ),
        pytest.param(
            "multi_entity_cif",
            {"P02281", "P0C0S5", "P62806", "P17317", "P84233", "Q7ZT64"},
            id="multiaccession-multichain",
        ),
    ],
)
def test_structure2uniprot_accessions(cif_fixture: str, expected: set[str], request: pytest.FixtureRequest):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)
    accessions = structure2uniprot_accessions(structure)

    assert accessions == expected


class TestAddUniprotAccessions2Structure:
    def test_none(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, None)
        assert structure == new_structure
        assert injected is False
        assert uniprot_chain_mappings == set()

    def test_missing_id(self, sample2_cif: Path, caplog: pytest.LogCaptureFixture):
        structure = read_structure(sample2_cif)
        pdb2uniprot = _mapping("1AAA", "P12345", "A=1-10")  # wrong PDB ID
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)

        assert structure == new_structure
        assert injected is False
        assert uniprot_chain_mappings == set()
        assert (
            "PDB ID 2Y29 not found in pdb2uniprot mapping. Leaving structure unverified and unchanged." in caplog.text
        )

    def test_verify_ok(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)
        pdb2uniprot = _mapping("2Y29", "P05067", "A=1-770")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert structure == new_structure
        assert injected is False
        assert uniprot_chain_mappings == set()

    def test_inject_into_nostructref(self, no_uniprot_cif: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = _mapping("2Y29", "P12345", "A=10-20")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A=10-20"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = structure_to_uniprot(new_structure)
        slim_result = {(r.chain_id, r.uniprot_accession) for r in result}

        expected = {("A", "P12345")}
        assert slim_result == expected

    def test_inject_into_existing_sifts(self, nmr_cif: Path):
        structure = read_structure(nmr_cif)
        pdb2uniprot = _mapping("1AMB", "P12345", "A=1-1000")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = structure_to_uniprot(new_structure)
        slim_result = {(r.chain_id, r.uniprot_accession) for r in result}
        expected = {("A", "P12345")}
        assert slim_result == expected

    def test_inject_multiple_ranges_into_nostructref(self, no_uniprot_cif: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = _mapping("2Y29", "P12345", "A=10-20,A=30-35")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A=10-20,A=30-35"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = uniprot_chain_mappings_from_struct_ref_seq(new_structure)
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A=10-20,A=30-35"),
            ),
        }
        assert result == expected

    def test_on_single_chain_written(self, multi_accession_cif: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        caplog.set_level("INFO")
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_file = write_single_chain_structure_file(
            input_file=multi_accession_cif, output_dir=input_dir, chain2keep="F"
        )

        structure = read_structure(input_file)
        pdb2uniprot = _mapping("1A02", "P01111", "F=1-1000")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P01111",
                chain_ranges=parse_uniprot_chains("A=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected
        result = structure_to_uniprot(new_structure)
        slim_result = {(r.chain_id, r.uniprot_accession) for r in result}
        expected = {("A", "P01111")}
        assert slim_result == expected

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
        pdb2uniprot = _mapping("1F66", "P12345", f"{chain_id}=1-1000")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            structure, pdb2uniprot, chain_system=chain_system
        )
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("C=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = structure_to_uniprot(new_structure, source="both")
        slim_result = {(r.chain_id, r.uniprot_accession) for r in result}
        expected = {
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
                # as requested, sifts entry P17317 gets ignored because P12345 is longer
                "P12345",
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
                "H",
                "P02281",
            ),
        }
        assert slim_result == expected

    def test_multi_entity_cif_multiple_auth_chains_in_single_mapping(self, multi_entity_cif: Path):
        structure = read_structure(multi_entity_cif)
        pdb2uniprot = _mapping("1F66", "P12345", "A/B=1-1000")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(structure, pdb2uniprot)
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A/B=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = structure_to_uniprot(new_structure, source="both")
        slim_result = {(r.chain_id, r.uniprot_accession) for r in result}
        expected = {("A", "P12345"), ("B", "P12345")}
        assert expected <= slim_result


@pytest.mark.parametrize(
    "cif_fixture, expected",
    [
        pytest.param("atomless_cif", set(), id="atomless"),
        pytest.param(
            "no_uniprot_cif",
            set(),
            id="no-uniprot",
        ),
        pytest.param(
            "cif_2fui",
            {
                UniprotChainMapping(
                    uniprot_accession="Q7Z7D6",
                    chain_ranges=parse_uniprot_chains("A=2583-2639"),
                )
            },
            id="1chain",
        ),
        pytest.param(
            "sample_multispan_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="O00255",
                    chain_ranges=parse_uniprot_chains("A=1-53,A=74-386,A=399-459,A=537-593"),
                )
            },
            id="multispan",
        ),
        pytest.param(
            "multi_accession_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="Q13469",
                    chain_ranges=parse_uniprot_chains("N=396-678"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P01100",
                    chain_ranges=parse_uniprot_chains("F=138-193"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05412",
                    chain_ranges=parse_uniprot_chains("J=253-308"),
                ),
            },
            id="multi-accession-separate-chains",
        ),
        pytest.param(
            "multi_accession_chain_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="P00656",
                    chain_ranges=parse_uniprot_chains("A=64-68"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P03950",
                    chain_ranges=parse_uniprot_chains("A=25-61,A=66-147"),
                ),
            },
            id="multi-accession-same-chain",
        ),
        pytest.param(
            "cif_8rw8",
            {
                UniprotChainMapping(
                    uniprot_accession="O00327",
                    chain_ranges=parse_uniprot_chains("B=337-449"),  # B is auth chain
                ),
            },
            id="authchain",
        ),
    ],
)
def test_uniprot_chain_mappings_from_struct_ref_seq(
    request: pytest.FixtureRequest, cif_fixture: str, expected: UniprotChainMappings
):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)

    actual = uniprot_chain_mappings_from_struct_ref_seq(structure)

    assert actual == expected


def test_uniprot_chain_mappings_from_struct_ref_seq_none_pos(
    sample2_cif: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    caplog.set_level(logging.INFO)
    new_fn = tmp_path / sample2_cif.stem
    with gzip.open(sample2_cif, "rt") as f_in:
        body = f_in.read()
        body = body.replace(
            "_struct_ref_seq.db_align_beg                  687", "_struct_ref_seq.db_align_beg                  ?"
        ).replace(
            "_struct_ref_seq.db_align_end                  692", "_struct_ref_seq.db_align_end                  ?"
        )
        new_fn.write_text(body)
    new_structure = read_structure(new_fn)

    actual = uniprot_chain_mappings_from_struct_ref_seq(new_structure)
    assert actual == set()
    assert "Skipping struct_ref_seq row with align_id" in caplog.text


@pytest.mark.parametrize(
    "cif_fixture, expected",
    [
        pytest.param("atomless_cif", set(), id="atomless"),
        pytest.param("sample2_cif", set(), id="siftless"),
        pytest.param(
            "cif_2fui",
            {
                UniprotChainMapping(
                    uniprot_accession="Q12830",
                    chain_ranges=parse_uniprot_chains("A=2865-2921"),
                )
            },
            id="1chain",
        ),
        pytest.param(
            "cif_3jrs",
            {
                UniprotChainMapping(
                    uniprot_accession="Q8VZS8",
                    chain_ranges=parse_uniprot_chains(
                        "A=31-158,A=164-209,B=30-157,B=165-209,C=38-51,C=55-138,C=142-157,C=164-209"
                    ),
                ),
            },
            id="multispan-multichain",
        ),
        pytest.param(
            "multi_entity_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="P0C0S5",
                    chain_ranges=parse_uniprot_chains("C=17-119,G=17-123"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P84233",
                    chain_ranges=parse_uniprot_chains("A=37-136,E=34-136"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P62806",
                    chain_ranges=parse_uniprot_chains("B=24-103,F=18-103"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P02281",
                    chain_ranges=parse_uniprot_chains("D=32-126,H=32-126"),
                ),
            },
            id="multiaccession-multichain",
        ),
        pytest.param(
            "cif_6o5i_updated",
            {
                UniprotChainMapping(
                    uniprot_accession="O00255",
                    chain_ranges=parse_uniprot_chains("A=2-53,A=74-386,A=399-459,A=549-588"),
                ),
            },
            id="multispan",
        ),
    ],
)
def test_uniprot_chain_mappings_from_sifts(
    request: pytest.FixtureRequest, cif_fixture: str, expected: UniprotChainMappings
):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)

    actual = uniprot_chain_mappings_from_sifts(structure)

    assert actual == expected


@pytest.mark.parametrize(
    "mappings, expected",
    [
        pytest.param(
            set(),
            set(),
            id="empty",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="Q12830",
                    chain_ranges=parse_uniprot_chains("A=2865-2921"),
                )
            },
            {
                FlattenedUniprotChainMapping("Q12830", 2865, 2921, "A", 1.0, 57),
            },
            id="1chain",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="Q8VZS8",
                    chain_ranges=parse_uniprot_chains(
                        "A=31-158,A=164-209,B=30-157,B=165-209,C=38-51,C=55-138,C=142-157,C=164-209"
                    ),
                ),
            },
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q8VZS8",
                    uniprot_start=31,
                    uniprot_end=209,
                    chain_id="A",
                    sequence_identity=0.9720670391061452,
                    aligned_residue_count=174,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q8VZS8",
                    uniprot_start=30,
                    uniprot_end=209,
                    chain_id="B",
                    sequence_identity=0.9611111111111111,
                    aligned_residue_count=173,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q8VZS8",
                    uniprot_start=38,
                    uniprot_end=209,
                    chain_id="C",
                    sequence_identity=0.9302325581395349,
                    aligned_residue_count=160,
                ),
            },
            id="multispan-multichain",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="P0C0S5",
                    chain_ranges=parse_uniprot_chains("C=17-119,G=17-123"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P84233",
                    chain_ranges=parse_uniprot_chains("A=37-136,E=34-136"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P62806",
                    chain_ranges=parse_uniprot_chains("B=24-103,F=18-103"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P02281",
                    chain_ranges=parse_uniprot_chains("D=32-126,H=32-126"),
                ),
            },
            {
                FlattenedUniprotChainMapping("P0C0S5", 17, 119, "C", 1.0, 103),
                FlattenedUniprotChainMapping("P0C0S5", 17, 123, "G", 1.0, 107),
                FlattenedUniprotChainMapping("P84233", 37, 136, "A", 1.0, 100),
                FlattenedUniprotChainMapping("P84233", 34, 136, "E", 1.0, 103),
                FlattenedUniprotChainMapping("P62806", 24, 103, "B", 1.0, 80),
                FlattenedUniprotChainMapping("P62806", 18, 103, "F", 1.0, 86),
                FlattenedUniprotChainMapping("P02281", 32, 126, "D", 1.0, 95),
                FlattenedUniprotChainMapping("P02281", 32, 126, "H", 1.0, 95),
            },
            id="multiaccession-multichain",
        ),
    ],
)
def test_flatten_uniprot_chain_mappings(mappings: UniprotChainMappings, expected: set[FlattenedUniprotChainMapping]):
    actual = flatten_uniprot_chain_mappings(mappings)

    assert actual == expected


@pytest.mark.parametrize(
    "mappings, expected",
    [
        pytest.param(
            set(),
            set(),
            id="empty",
        ),
        pytest.param(
            {
                FlattenedUniprotChainMapping("P12345", 10, 15, "A", 1.0, 6),
            },
            {
                FlattenedUniprotChainMapping("P12345", 10, 15, "A", 1.0, 6),
            },
            id="solo",
        ),
        pytest.param(
            {
                FlattenedUniprotChainMapping("P12345", 10, 15, "A", 1.0, 6),
                FlattenedUniprotChainMapping("P67879", 20, 24, "A", 1.0, 5),
            },
            {
                FlattenedUniprotChainMapping("P12345", 10, 15, "A", 1.0, 6),
            },
            id="longest",
        ),
        pytest.param(
            {
                FlattenedUniprotChainMapping("P12345", 10, 14, "A", 1.0, 5),
                FlattenedUniprotChainMapping("P67879", 20, 24, "A", 1.0, 5),
            },
            {
                FlattenedUniprotChainMapping("P12345", 10, 14, "A", 1.0, 5),
            },
            id="same-length-first-alpha",
        ),
        pytest.param(
            {
                FlattenedUniprotChainMapping("P12345", 10, 15, "A", 1.0, 6),
                FlattenedUniprotChainMapping("P67879", 20, 24, "B", 1.0, 5),
            },
            {
                FlattenedUniprotChainMapping("P12345", 10, 15, "A", 1.0, 6),
                FlattenedUniprotChainMapping("P67879", 20, 24, "B", 1.0, 5),
            },
            id="2chain-1each",
        ),
    ],
)
def test_best_uniprot_per_chain(
    mappings: set[FlattenedUniprotChainMapping], expected: set[FlattenedUniprotChainMapping]
):
    actual = best_uniprot_per_chain(mappings)

    assert actual == expected


@pytest.mark.parametrize(
    "cif_fixture, source, expected",
    [
        pytest.param(
            "atomless_cif",
            "both",
            set(),
            id="atomless",
        ),
        pytest.param(
            "sample2_cif",
            "both",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="P05067",
                    uniprot_start=687,
                    uniprot_end=692,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=6,
                ),
            },
            id="1chain-1uniprot",
        ),
        pytest.param(
            "cif_2fui",
            "sifts",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q12830",
                    uniprot_start=2865,
                    uniprot_end=2921,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=57,
                ),
            },
            id="1chain-1uniprot-sifts",
        ),
        pytest.param(
            "cif_2fui",
            "struct_ref_seq",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q7Z7D6",
                    uniprot_start=2583,
                    uniprot_end=2639,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=57,
                ),
            },
            id="1chain-1uniprot-struct_ref_seq",
        ),
        pytest.param(
            "cif_2fui",
            "both",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q12830",  # sifts wins
                    uniprot_start=2865,
                    uniprot_end=2921,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=57,
                ),
            },
            id="1chain-1uniprot-both",
        ),
        pytest.param(
            "cif_2fui",
            "fallback",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q12830",  # sifts wins
                    uniprot_start=2865,
                    uniprot_end=2921,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=57,
                ),
            },
            id="1chain-1uniprot-fallback",
        ),
        pytest.param(
            "multi_entity_cif",
            "both",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q7ZT64",
                    uniprot_start=1,
                    uniprot_end=136,
                    chain_id="A",
                    sequence_identity=2.0,
                    aligned_residue_count=272,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P62806",
                    uniprot_start=1,
                    uniprot_end=102,
                    chain_id="F",
                    sequence_identity=1.0,
                    aligned_residue_count=102,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P62806",
                    uniprot_start=1,
                    uniprot_end=102,
                    chain_id="B",
                    sequence_identity=1.0,
                    aligned_residue_count=102,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P17317",
                    uniprot_start=1,
                    uniprot_end=127,
                    chain_id="C",
                    sequence_identity=1.0,
                    aligned_residue_count=127,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P02281",
                    uniprot_start=1,
                    uniprot_end=125,
                    chain_id="D",
                    sequence_identity=1.0,
                    aligned_residue_count=125,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="Q7ZT64",
                    uniprot_start=1,
                    uniprot_end=136,
                    chain_id="E",
                    sequence_identity=2.0,
                    aligned_residue_count=272,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P02281",
                    uniprot_start=1,
                    uniprot_end=125,
                    chain_id="H",
                    sequence_identity=1.0,
                    aligned_residue_count=125,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P17317",
                    uniprot_start=1,
                    uniprot_end=127,
                    chain_id="G",
                    sequence_identity=1.0,
                    aligned_residue_count=127,
                ),
            },
            id="multiaccession-multichain",
        ),
        pytest.param(
            "multi_accession_chain_cif",
            "both",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="P03950",
                    uniprot_start=25,
                    uniprot_end=147,
                    chain_id="A",
                    sequence_identity=0.967479674796748,
                    aligned_residue_count=119,
                ),
            },
            id="multi-accession-same-chain",
        ),
    ],
)
def test_structure_to_uniprot(
    request: pytest.FixtureRequest, cif_fixture: str, source: UniprotSource, expected: set[FlattenedUniprotChainMapping]
):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)

    actual = structure_to_uniprot(structure, source=source)

    assert actual == expected


def test_structure_to_uniprot_allow_multiple_accessions_per_chain(multi_accession_chain_cif: Path):
    structure = read_structure(multi_accession_chain_cif)

    actual = structure_to_uniprot(structure, source="both", one_uniprot_per_chain=False)

    expected = {
        FlattenedUniprotChainMapping(
            uniprot_accession="P00656",
            uniprot_start=64,
            uniprot_end=68,
            chain_id="A",
            sequence_identity=1.0,
            aligned_residue_count=5,
        ),
        FlattenedUniprotChainMapping(
            uniprot_accession="P03950",
            uniprot_start=25,
            uniprot_end=147,
            chain_id="A",
            sequence_identity=0.967479674796748,
            aligned_residue_count=119,
        ),
    }

    assert actual == expected


def test_structure_to_uniprot_bad_source(sample2_cif: Path):
    structure = read_structure(sample2_cif)

    with pytest.raises(ValueError, match="Invalid source 'badsource'"):
        # pyrefly: ignore [bad-argument-type]
        structure_to_uniprot(structure, source="badsource")


@pytest.mark.parametrize(
    "mappings, chain_provenance, expected",
    [
        pytest.param(
            set(),
            ChainExtractionProvenance("F", "A"),
            set(),
            id="empty",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("F=1-10"),
                )
            },
            ChainExtractionProvenance("F", "A"),
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("A=1-10"),
                )
            },
            id="rename",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("F=1-10"),
                )
            },
            ChainExtractionProvenance("G", "A"),
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("F=1-10"),
                )
            },
            id="no-rename",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("B/F=1-10"),
                )
            },
            ChainExtractionProvenance("F", "A"),
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("A/B=1-10"),
                )
            },
            id="rename-half",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("F=1-10"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P6789",
                    chain_ranges=parse_uniprot_chains("B=1-10"),
                ),
            },
            ChainExtractionProvenance("F", "A"),
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("A=1-10"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P6789",
                    chain_ranges=parse_uniprot_chains("B=1-10"),
                ),
            },
            id="ignores-other",
        ),
        pytest.param(
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("F=1-10,F=20-30"),
                )
            },
            ChainExtractionProvenance("F", "A"),
            {
                UniprotChainMapping(
                    uniprot_accession="P12345",
                    chain_ranges=parse_uniprot_chains("A=1-10,A=20-30"),
                )
            },
            id="multi-range",
        ),
    ],
)
def test_apply_chain_provenance_to_uniprot_mappings(
    mappings: UniprotChainMappings, chain_provenance: ChainExtractionProvenance, expected: UniprotChainMappings
):
    actual = apply_chain_provenance_to_uniprot_mappings(mappings, chain_provenance)
    assert actual == expected
