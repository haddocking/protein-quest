from pathlib import Path

import gemmi
import pytest

from protein_quest.structure.chains import ChainExtractionProvenance, ChainIdSystem, write_single_chain_structure_file
from protein_quest.structure.formats import read_structure, write_structure
from protein_quest.structure.uniprot_extraction import (
    FlattenedUniprotChainMapping,
    UniprotSource,
    structure_to_uniprot,
    uniprot_chain_mappings_from_struct_ref_seq,
)
from protein_quest.structure.uniprot_injection import (
    add_uniprot_accessions2structure,
    apply_chain_provenance_to_uniprot_mappings,
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


def _write_and_extract_uniprots(
    tmp_path: Path,
    output_name: str,
    structure: gemmi.Structure,
    *,
    source: UniprotSource = "both",
) -> set[FlattenedUniprotChainMapping]:
    output_file = tmp_path / output_name
    write_structure(structure, output_file)
    return structure_to_uniprot(output_file, source=source)


class TestAddUniprotAccessions2Structure:
    def test_none(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            sample2_cif, structure=structure
        )
        assert structure == new_structure
        assert injected is False
        assert uniprot_chain_mappings == set()

    def test_missing_id(self, sample2_cif: Path, caplog: pytest.LogCaptureFixture):
        structure = read_structure(sample2_cif)
        pdb2uniprot = _mapping("1AAA", "P12345", "A=1-10")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            sample2_cif, pdb2uniprot, structure=structure
        )

        assert structure == new_structure
        assert injected is False
        assert uniprot_chain_mappings == set()
        assert (
            "PDB ID 2Y29 not found in pdb2uniprot mapping. Leaving structure unverified and unchanged." in caplog.text
        )

    def test_verify_ok(self, sample2_cif: Path):
        structure = read_structure(sample2_cif)
        pdb2uniprot = _mapping("2Y29", "P05067", "A=1-770")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            sample2_cif, pdb2uniprot, structure=structure
        )
        assert structure == new_structure
        assert injected is False
        assert uniprot_chain_mappings == set()

    def test_inject_into_nostructref(self, no_uniprot_cif: Path, tmp_path: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = _mapping("2Y29", "P12345", "A=10-20")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            no_uniprot_cif, pdb2uniprot, structure=structure
        )
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A=10-20"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = _write_and_extract_uniprots(tmp_path, "no_uniprot_injected.cif.gz", new_structure)
        slim_result = {(mapping.chain_id, mapping.uniprot_accession) for mapping in result}

        expected = {("A", "P12345")}
        assert slim_result == expected

    def test_inject_into_existing_sifts(self, nmr_cif: Path, tmp_path: Path):
        structure = read_structure(nmr_cif)
        pdb2uniprot = _mapping("1AMB", "P12345", "A=1-1000")
        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            nmr_cif, pdb2uniprot, structure=structure
        )
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = _write_and_extract_uniprots(tmp_path, "nmr_injected.cif.gz", new_structure)
        slim_result = {(mapping.chain_id, mapping.uniprot_accession) for mapping in result}
        expected = {("A", "P12345")}
        assert slim_result == expected

    def test_inject_multiple_ranges_into_nostructref(self, no_uniprot_cif: Path):
        structure = read_structure(no_uniprot_cif)
        pdb2uniprot = _mapping("2Y29", "P12345", "A=10-20,A=30-35")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            no_uniprot_cif, pdb2uniprot, structure=structure
        )
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

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            input_file, pdb2uniprot, structure=structure
        )
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P01111",
                chain_ranges=parse_uniprot_chains("A=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected
        result = _write_and_extract_uniprots(tmp_path, "single_chain_injected.cif.gz", new_structure)
        slim_result = {(mapping.chain_id, mapping.uniprot_accession) for mapping in result}
        expected = {("A", "P01111")}
        assert slim_result == expected

        log = caplog.text
        assert "Structure 1A02 has provenance information indicating it was extracted from chain F to A" in log
        assert "Using this information to verify/add UniProt accessions." in log

    @pytest.mark.parametrize(
        ("chain_system", "chain_id"),
        [
            ("auth", "C"),
            ("label", "E"),
        ],
    )
    def test_multi_entity_cif(self, multi_entity_cif: Path, chain_system: ChainIdSystem, chain_id: str, tmp_path: Path):
        structure = read_structure(multi_entity_cif)
        pdb2uniprot = _mapping("1F66", "P12345", f"{chain_id}=1-1000")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            multi_entity_cif, pdb2uniprot, structure=structure, chain_system=chain_system
        )
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("C=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = _write_and_extract_uniprots(tmp_path, f"multi_entity_{chain_system}.cif.gz", new_structure)
        slim_result = {(mapping.chain_id, mapping.uniprot_accession) for mapping in result}
        expected = {
            ("A", "Q7ZT64"),
            ("B", "P62806"),
            ("C", "P12345"),
            ("D", "P02281"),
            ("E", "Q7ZT64"),
            ("F", "P62806"),
            ("G", "P12345"),
            ("H", "P02281"),
        }
        assert slim_result == expected

    def test_multi_entity_cif_multiple_auth_chains_in_single_mapping(self, multi_entity_cif: Path, tmp_path: Path):
        structure = read_structure(multi_entity_cif)
        pdb2uniprot = _mapping("1F66", "P12345", "A/B=1-1000")

        new_structure, injected, uniprot_chain_mappings = add_uniprot_accessions2structure(
            multi_entity_cif, pdb2uniprot, structure=structure
        )
        assert injected is True
        expected = {
            UniprotChainMapping(
                uniprot_accession="P12345",
                chain_ranges=parse_uniprot_chains("A/B=1-1000"),
            ),
        }
        assert uniprot_chain_mappings == expected

        result = _write_and_extract_uniprots(tmp_path, "multi_entity_multi_auth.cif.gz", new_structure)
        slim_result = {(mapping.chain_id, mapping.uniprot_accession) for mapping in result}
        expected = {("A", "P12345"), ("B", "P12345")}
        assert expected <= slim_result


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
