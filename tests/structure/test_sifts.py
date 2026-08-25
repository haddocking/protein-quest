import logging
from pathlib import Path

import gemmi
import pytest

from protein_quest.structure.formats import read_structure_as_cif_block
from protein_quest.structure.sifts import uniprot_chain_mappings_from_cif, uniprot_chain_mappings_to_cif
from protein_quest.uniprot_chains import UniprotChainMapping, UniprotChainMappings, parse_uniprot_chains


@pytest.mark.parametrize(
    "cif_fixture, expected",
    [
        pytest.param("atomless_cif", set(), id="atomless"),
        pytest.param("af_pdb", set(), id="not-a-cif"),
        pytest.param("comments_only_cif", set(), id="no-blocks"),
        pytest.param("sample2_cif", set(), id="no-sifts-segments"),
        pytest.param(
            "em_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="P0ABE7",
                    chain_ranges=parse_uniprot_chains("A=4-13,A=23-127"),
                )
            },
            id="multispan-em",
        ),
        pytest.param(
            "nmr_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="P05067",
                    chain_ranges=parse_uniprot_chains("A=672-699"),
                )
            },
            id="best-mapping-only-nmr-isoforms",
        ),
        pytest.param(
            "cif_2fui",
            {
                UniprotChainMapping(
                    uniprot_accession="Q12830",
                    chain_ranges=parse_uniprot_chains("A=2865-2921"),
                ),
            },
            id="best-mapping-only",
        ),
        pytest.param(
            "cif_3jrs",
            {
                UniprotChainMapping(
                    uniprot_accession="Q8VZS8",
                    chain_ranges=parse_uniprot_chains("A=8-211,B=8-211,C=8-211"),
                ),
            },
            id="single-span-multichain",
        ),
        pytest.param(
            "cif_6o5i_updated",
            {
                UniprotChainMapping(
                    uniprot_accession="O00255",
                    chain_ranges=parse_uniprot_chains("A=1-593"),
                ),
            },
            id="best-mapping-only-multisegment-entry",
        ),
        pytest.param(
            "cif_8rw8",
            {
                UniprotChainMapping(
                    uniprot_accession="O00327",
                    chain_ranges=parse_uniprot_chains("A=337-449"),
                ),
            },
            id="best-mapping-only-label-auth-mismatch-raw-asym-id",
        ),
        pytest.param(
            "cif_3plz",
            {
                UniprotChainMapping(
                    uniprot_accession="Q15596",
                    chain_ranges=parse_uniprot_chains("B=740-753,D=740-753"),
                ),
                UniprotChainMapping(
                    uniprot_accession="O00482",
                    chain_ranges=parse_uniprot_chains("A=300-541,C=300-541"),
                ),
            },
            id="best-mapping-only-multiple-accessions-multichain",
        ),
    ],
)
def test_uniprot_chain_mappings_from_cif(
    request: pytest.FixtureRequest, cif_fixture: str, expected: UniprotChainMappings
):
    path = request.getfixturevalue(cif_fixture)
    block = read_structure_as_cif_block(path)

    actual = uniprot_chain_mappings_from_cif(block) if block is not None else set()

    assert actual == expected


@pytest.mark.parametrize(
    "cif_fixture, expected",
    [
        pytest.param(
            "cif_3plz",
            {
                UniprotChainMapping(
                    uniprot_accession="Q15596",
                    chain_ranges=parse_uniprot_chains("B=740-753,D=740-753"),
                ),
                UniprotChainMapping(
                    uniprot_accession="O00482-2",
                    chain_ranges=parse_uniprot_chains("A=254-495,C=254-495"),
                ),
                UniprotChainMapping(
                    uniprot_accession="O00482",
                    chain_ranges=parse_uniprot_chains("A=300-541,C=300-541"),
                ),
                UniprotChainMapping(
                    uniprot_accession="O00482-3",
                    chain_ranges=parse_uniprot_chains("A=198-369,C=198-369"),
                ),
                UniprotChainMapping(
                    uniprot_accession="O00482-4",
                    chain_ranges=parse_uniprot_chains("A=228-469,C=228-469"),
                ),
            },
            id="all-segments-multiple-accessions-multichain",
        ),
        pytest.param(
            "nmr_cif",
            {
                UniprotChainMapping(
                    uniprot_accession="P05067-3",
                    chain_ranges=parse_uniprot_chains("A=579-606"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-4",
                    chain_ranges=parse_uniprot_chains("A=597-624"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-5",
                    chain_ranges=parse_uniprot_chains("A=598-625"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-8",
                    chain_ranges=parse_uniprot_chains("A=653-680"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067",
                    chain_ranges=parse_uniprot_chains("A=672-699"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-10",
                    chain_ranges=parse_uniprot_chains("A=541-568"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-6",
                    chain_ranges=parse_uniprot_chains("A=616-643"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-7",
                    chain_ranges=parse_uniprot_chains("A=635-662"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-11",
                    chain_ranges=parse_uniprot_chains("A=648-675"),
                ),
                UniprotChainMapping(
                    uniprot_accession="P05067-9",
                    chain_ranges=parse_uniprot_chains("A=654-681"),
                ),
            },
            id="all-segments-nmr-isoforms",
        ),
    ],
)
def test_uniprot_chain_mappings_from_cif_best_only_false(
    request: pytest.FixtureRequest, cif_fixture: str, expected: UniprotChainMappings
):
    path = request.getfixturevalue(cif_fixture)
    block = read_structure_as_cif_block(path)

    actual = (
        uniprot_chain_mappings_from_cif(block, best_only=False, strip_isoform=False) if block is not None else set()
    )

    assert actual == expected


def test_uniprot_chain_mappings_from_cif_skips_non_numeric_segment_bounds(
    malformed_sifts_segments_cif: Path, caplog: pytest.LogCaptureFixture
):
    caplog.set_level(logging.INFO)
    block = read_structure_as_cif_block(malformed_sifts_segments_cif)

    assert block is not None

    actual = uniprot_chain_mappings_from_cif(block)

    assert actual == {
        UniprotChainMapping(
            uniprot_accession="P0ABE7",
            chain_ranges=parse_uniprot_chains("A=23-127"),
        )
    }
    assert "Skipping pdbx_sifts_unp_segments row with accession P0ABE7" in caplog.text
    assert block.name in caplog.text


def test_uniprot_chain_mappings_to_cif_writes_sifts_segments():
    block = gemmi.cif.Document().add_new_block("test")
    mappings = {
        UniprotChainMapping(
            uniprot_accession="Q12830",
            chain_ranges=parse_uniprot_chains("B/A=2865-2921"),
        ),
        UniprotChainMapping(
            uniprot_accession="P05067",
            chain_ranges=parse_uniprot_chains("C=672-699"),
        ),
    }

    uniprot_chain_mappings_to_cif(mappings, block)

    assert block.get_mmcif_category("_pdbx_sifts_unp_segments.") == {
        "asym_id": ["C", "B", "A"],
        "unp_acc": ["P05067", "Q12830", "Q12830"],
        "unp_start": ["672", "2865", "2865"],
        "unp_end": ["699", "2921", "2921"],
        "best_mapping": ["y", "y", "y"],
    }
    assert uniprot_chain_mappings_from_cif(block) == {
        UniprotChainMapping(
            uniprot_accession="Q12830",
            chain_ranges=parse_uniprot_chains("B=2865-2921,A=2865-2921"),
        ),
        UniprotChainMapping(
            uniprot_accession="P05067",
            chain_ranges=parse_uniprot_chains("C=672-699"),
        ),
    }
