import gzip
import logging
from pathlib import Path

import pytest

from protein_quest.structure.formats import read_structure
from protein_quest.structure.uniprot_extraction import (
    FlattenedUniprotChainMapping,
    UniprotSource,
    best_uniprot_per_chain,
    flatten_uniprot_chain_mappings,
    structure2uniprot_accessions,
    structure_to_uniprot,
    uniprot_chain_mappings_from_struct_ref_seq,
)
from protein_quest.uniprot_chains import UniprotChainMapping, UniprotChainMappings, parse_uniprot_chains


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
    accessions = structure2uniprot_accessions(path)

    assert accessions == expected


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
                    chain_ranges=parse_uniprot_chains("B=337-449"),
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
                    uniprot_accession="Q12830",
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
                    uniprot_accession="Q12830",
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
                    uniprot_accession="P02281",
                    uniprot_start=1,
                    uniprot_end=126,
                    chain_id="D",
                    sequence_identity=1.0,
                    aligned_residue_count=126,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P62806",
                    uniprot_start=1,
                    uniprot_end=103,
                    chain_id="B",
                    sequence_identity=1.0,
                    aligned_residue_count=103,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P02281",
                    uniprot_start=1,
                    uniprot_end=126,
                    chain_id="H",
                    sequence_identity=1.0,
                    aligned_residue_count=126,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P0C0S5",
                    uniprot_start=1,
                    uniprot_end=128,
                    chain_id="C",
                    sequence_identity=1.0,
                    aligned_residue_count=128,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P0C0S5",
                    uniprot_start=1,
                    uniprot_end=128,
                    chain_id="G",
                    sequence_identity=1.0,
                    aligned_residue_count=128,
                ),
                FlattenedUniprotChainMapping(
                    uniprot_accession="P62806",
                    uniprot_start=1,
                    uniprot_end=103,
                    chain_id="F",
                    sequence_identity=1.0,
                    aligned_residue_count=103,
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
                    uniprot_accession="Q7ZT64",
                    uniprot_start=1,
                    uniprot_end=136,
                    chain_id="A",
                    sequence_identity=2.0,
                    aligned_residue_count=272,
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
        pytest.param(
            "cif_4gpq",
            "sifts",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="O00255",
                    uniprot_start=1,
                    uniprot_end=598,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=598,
                ),
            },
            id="4gpq-both",
        ),
        pytest.param(
            "em_cif",
            "both",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="P0ABE7",
                    uniprot_start=4,
                    uniprot_end=127,
                    chain_id="A",
                    sequence_identity=0.9274193548387096,
                    aligned_residue_count=115,
                )
            },
            id="em-cif-both",
        ),
        pytest.param(
            "em_cif",
            "fallback",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="P0ABE7",
                    uniprot_start=4,
                    uniprot_end=127,
                    chain_id="A",
                    sequence_identity=0.9274193548387096,
                    aligned_residue_count=115,
                )
            },
            id="em-cif-fallback",
        ),
        pytest.param(
            "em_cif",
            "sifts",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="P0ABE7",
                    uniprot_start=4,
                    uniprot_end=127,
                    chain_id="A",
                    sequence_identity=0.9274193548387096,
                    aligned_residue_count=115,
                )
            },
            id="em-cif-sifts",
        ),
        pytest.param(
            "em_cif",
            "struct_ref_seq",
            {
                FlattenedUniprotChainMapping(
                    uniprot_accession="P0ABE7",
                    uniprot_start=23,
                    uniprot_end=127,
                    chain_id="A",
                    sequence_identity=1.0,
                    aligned_residue_count=105,
                )
            },
            id="em-cif-struct_ref_seq",
        ),
    ],
)
def test_structure_to_uniprot(
    request: pytest.FixtureRequest, cif_fixture: str, source: UniprotSource, expected: set[FlattenedUniprotChainMapping]
):
    path = request.getfixturevalue(cif_fixture)
    structure = read_structure(path)

    actual = structure_to_uniprot(path, structure=structure, source=source)

    assert actual == expected


def test_structure_to_uniprot_allow_multiple_accessions_per_chain(multi_accession_chain_cif: Path):
    structure = read_structure(multi_accession_chain_cif)

    actual = structure_to_uniprot(
        multi_accession_chain_cif, structure=structure, source="both", one_uniprot_per_chain=False
    )

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


def test_structure_to_uniprot_uses_raw_cif_sifts_segments_when_path_provided(em_cif: Path):
    structure = read_structure(em_cif)

    actual = structure_to_uniprot(em_cif, structure=structure, source="sifts")

    assert actual == {
        FlattenedUniprotChainMapping(
            uniprot_accession="P0ABE7",
            uniprot_start=4,
            uniprot_end=127,
            chain_id="A",
            sequence_identity=0.9274193548387096,
            aligned_residue_count=115,
        )
    }


def test_structure_to_uniprot_logs_isoform_stripping(cif_4gpq: Path, caplog: pytest.LogCaptureFixture):
    caplog.set_level(logging.INFO)
    structure = read_structure(cif_4gpq)

    _ = structure_to_uniprot(cif_4gpq, structure=structure, source="sifts")

    assert "_pdbx_sifts_unp_segments.unp_acc is isoform, stripping to base accession: O00255-1 -> O00255" in caplog.text


def test_structure_to_uniprot_issue_155_keeps_all_accessions_for_single_chain(cif_3plz: Path):
    """Regression test for issue #155: mixed SIFTS + struct_ref_seq data should keep both accessions."""
    structure = read_structure(cif_3plz)

    actual = structure_to_uniprot(cif_3plz, structure=structure, source="both", one_uniprot_per_chain=False)

    assert {(mapping.chain_id, mapping.uniprot_accession) for mapping in actual if mapping.chain_id == "A"} == {
        ("A", "Q9UEC0"),
        ("A", "O00482"),
    }


def test_structure_to_uniprot_bad_source(sample2_cif: Path):
    structure = read_structure(sample2_cif)

    with pytest.raises(ValueError, match="Invalid source 'badsource'"):
        # pyrefly: ignore [bad-argument-type]
        structure_to_uniprot(sample2_cif, structure=structure, source="badsource")
