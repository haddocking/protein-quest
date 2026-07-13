import csv
from dataclasses import replace
from pathlib import Path

import gemmi
import pytest
from cyclopts.types import StdioPath

from protein_quest.alphafold.confidence import ConfidenceFilterQuery
from protein_quest.filters.combined import (
    CombinedFilterQuery,
    CombinedFilterResult,
    CombinedPartitions,
    combined_filter,
    combined_filter_stats,
    combined_filter_summary,
)
from protein_quest.filters.quality import Scores
from protein_quest.structure.formats import read_structure
from protein_quest.structure.metadata import StructureMetadata


class TestCombinedFilterQuery:
    def test_confidence_filter_query_with_defaults(self):
        actual = CombinedFilterQuery().confidence_filter_query()

        expected = ConfidenceFilterQuery()
        assert actual == expected


class TestCombinedPartitions:
    def test_extend_empty(self):
        root = CombinedPartitions()
        others = [CombinedPartitions(), CombinedPartitions()]

        actual = root.extend(others)

        expected = CombinedPartitions()
        assert actual == expected


def strip_resolution(input_file: Path, output_file: Path):
    s = read_structure(input_file)
    doc = s.make_mmcif_document(gemmi.MmcifOutputGroups(True, chem_comp=False))
    block = doc.sole_block()
    block.set_pair("_refine.ls_d_res_high", "?")
    body = doc.as_string()
    output_file.write_text(body)


@pytest.mark.skip("Remnant for developing strip_resolution")
def test_strip_resolution(tmp_path: Path, sample2_cif: Path):
    output_file = tmp_path / "output.cif"

    strip_resolution(sample2_cif, output_file)

    s = read_structure(output_file)
    assert s.resolution == 0.0


class TestCombinedFilter:
    def test_with_all_good_cifs_without_scores(self, all_cifs: list[Path], tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_files = []
        for cif in all_cifs:
            input_file = input_dir / cif.name
            input_file.hardlink_to(cif)
            input_files.append(input_file)
        scores: dict[str, Scores] = {}
        query = CombinedFilterQuery(min_sequence_identity=0.8)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        results = combined_filter(input_files, scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz",
                pdb_id="AF-A0A0C5B5G6-F1",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=10,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "AF-A0A0C5B5G6-F1-model_v6.cif.gz",
            ),
            CombinedFilterResult(
                input_file=input_dir / "1amb_updated.cif.gz",
                pdb_id="1AMB",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason="No UniProt accession, no resolution, no geometry quality",
                output_file=None,
            ),
            CombinedFilterResult(
                input_file=input_dir / "3jrs_updated_B2A.cif.gz",
                pdb_id="3JRS",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "3jrs_updated_B2A.cif.gz",
            ),
            CombinedFilterResult(
                input_file=input_dir / "2Y29.cif.gz",
                pdb_id="2Y29",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "2Y29.cif.gz",
            ),
            CombinedFilterResult(
                input_file=input_dir / "6O5I.cif.gz",
                pdb_id="6O5I",
                metadata=StructureMetadata(
                    id="6O5I",
                    uniprot_accession="O00255",
                    resolution=1.24,
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "6O5I.cif.gz",
            ),
            CombinedFilterResult(
                input_file=input_dir / "1un5.cif.gz",
                pdb_id="1UN5",
                metadata=StructureMetadata(
                    id="1UN5",
                    uniprot_accession="P03950",
                    resolution=2.6,
                    total_residue_count=131,
                    is_alphafold=False,
                    uniprot_start=25,
                    uniprot_end=147,
                    sequence_identity=0.967,
                    chain_length=131,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "1un5.cif.gz",
            ),
            CombinedFilterResult(
                input_file=input_dir / "8w77_updated.cif.gz",
                pdb_id="8W77",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "8w77_updated.cif.gz",
            ),
        ]
        assert results == expected

    def test_with_bad_cifs_and_defaults(self, bad_cifs: list[Path], tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_files = []
        for cif in bad_cifs:
            input_file = input_dir / cif.name
            input_file.hardlink_to(cif)
            input_files.append(input_file)
        scores: dict[str, Scores] = {}
        query = CombinedFilterQuery()
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        results = combined_filter(input_files, scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_dir / "unreadable.cif",
                pdb_id=None,
                metadata=None,
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason=f"Failed to read structure: {input_dir / 'unreadable.cif'}:1:0(0): expected block header (data_)",
                output_file=None,
            ),
            CombinedFilterResult(
                input_file=input_dir / "atomless.cif",
                pdb_id=" ",
                metadata=None,
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason="Failed to extract metadata: No chains found in structure  ",
                output_file=None,
            ),
            CombinedFilterResult(
                input_file=input_dir / "no_uniprot.cif",
                pdb_id="2Y29",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason="Sorted index 1 > 0",
                output_file=None,
            ),
        ]
        assert results == expected

    def test_with_geometry_quality(self, nmr_cif: Path, tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_file = input_dir / nmr_cif.name
        input_file.hardlink_to(nmr_cif)
        scores: dict[str, Scores] = {
            "1amb": Scores(
                geometry_quality=70.3,
                data_quality=71,
                experiment_data_available=False,
                overall_quality=70.3,
            )
        }
        query = CombinedFilterQuery()
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        results = combined_filter([input_file], scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_file,
                pdb_id="1AMB",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=70.3,
                passed=True,
                reason=None,
                output_file=output_dir / nmr_cif.name,
            )
        ]
        assert results == expected

    def test_af_too_short(self, af_cif: Path, tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_file = input_dir / af_cif.name
        input_file.hardlink_to(af_cif)
        scores: dict[str, Scores] = {}
        query = CombinedFilterQuery(
            min_residues=20,
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        results = combined_filter([input_file], scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_file,
                pdb_id="AF-A0A0C5B5G6-F1",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=10,
                geometry_quality=None,
                passed=False,
                reason="Low confidence or too few or too many residues (10)",
                output_file=None,
            )
        ]
        assert results == expected

    def test_xray_too_short(self, sample2_cif: Path, tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_file = input_dir / sample2_cif.name
        input_file.hardlink_to(sample2_cif)
        scores: dict[str, Scores] = {}
        query = CombinedFilterQuery(
            min_residues=20,
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        results = combined_filter([input_file], scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_file,
                pdb_id="2Y29",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason="Chain length 8 not in range [20, 10000000]",
                output_file=None,
            )
        ]
        assert results == expected

    def test_top_uniprot_cluster_with_top1_and_xrays(self, xray_p05067_cifs: list[Path], tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_files = []
        for cif in xray_p05067_cifs:
            input_file = input_dir / cif.name
            input_file.hardlink_to(cif)
            input_files.append(input_file)
        scores: dict[str, Scores] = {}
        query = CombinedFilterQuery(
            top_uniprot_cluster=1,
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        actual = combined_filter(input_files, scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_dir / "2y2a_updated.cif.gz",
                pdb_id="2Y2A",
                metadata=StructureMetadata(
                    id="2Y2A",
                    uniprot_accession="P05067",
                    resolution=1.91,
                    total_residue_count=10,
                    is_alphafold=False,
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=10,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "2y2a_updated.cif.gz",
            ),
            CombinedFilterResult(
                input_file=input_dir / "2Y29.cif.gz",
                pdb_id="2Y29",
                metadata=StructureMetadata(
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
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason="Sorted index 2 > 1",
                output_file=None,
            ),
        ]
        assert actual == expected

    def setup_resolutionless_xrays(self, xray_p05067_cifs: list[Path], tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_files = []
        for cif in xray_p05067_cifs:
            # Clear resolution so geometry_quality branch is taken
            input_file = input_dir / cif.stem
            strip_resolution(cif, input_file)
            input_files.append(input_file)
        scores: dict[str, Scores] = {
            "8t89": Scores(
                geometry_quality=70.3,
                data_quality=71,
                experiment_data_available=False,
                overall_quality=70.3,
            ),
            "2y2a": Scores(
                geometry_quality=68.5,
                data_quality=70,
                experiment_data_available=False,
                overall_quality=68.5,
            ),
            "2y29": Scores(
                geometry_quality=85.0,
                data_quality=89,
                experiment_data_available=False,
                overall_quality=85.0,
            ),
        }
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        return input_dir, input_files, scores, output_dir

    def test_top_uniprot_cluster_with_top1_and_resolutionless_xrays(self, xray_p05067_cifs: list[Path], tmp_path: Path):
        input_dir, input_files, scores, output_dir = self.setup_resolutionless_xrays(xray_p05067_cifs, tmp_path)
        query = CombinedFilterQuery(
            top_uniprot_cluster=1,
        )

        actual = combined_filter(input_files, scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_dir / "2Y29.cif",
                pdb_id="2Y29",
                metadata=StructureMetadata(
                    id="2Y29",
                    uniprot_accession="P05067",
                    resolution=0.0,
                    total_residue_count=8,
                    is_alphafold=False,
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=8,
                    auth_chain="A",
                    label_chain="A",
                    method="Other",
                ),
                high_confidence_residues_count=None,
                geometry_quality=85.0,
                passed=True,
                reason=None,
                output_file=output_dir / "2Y29.cif",
            ),
            CombinedFilterResult(
                input_file=input_dir / "2y2a_updated.cif",
                pdb_id="2Y2A",
                metadata=StructureMetadata(
                    id="2Y2A",
                    uniprot_accession="P05067",
                    resolution=0.0,
                    total_residue_count=10,
                    is_alphafold=False,
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=10,
                    auth_chain="A",
                    label_chain="A",
                    method="Other",
                ),
                high_confidence_residues_count=None,
                geometry_quality=68.5,
                passed=False,
                reason="Sorted index 2 > 1",
                output_file=None,
            ),
        ]
        assert actual == expected

    def test_top_uniprot_cluster_with_high_qual_and_resolutionless_xrays(
        self, xray_p05067_cifs: list[Path], tmp_path: Path
    ):
        input_dir, input_files, scores, output_dir = self.setup_resolutionless_xrays(xray_p05067_cifs, tmp_path)
        query = CombinedFilterQuery(
            min_geometry_quality=80.0,
        )

        actual = combined_filter(input_files, scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_dir / "2Y29.cif",
                pdb_id="2Y29",
                metadata=StructureMetadata(
                    id="2Y29",
                    uniprot_accession="P05067",
                    resolution=0.0,
                    total_residue_count=8,
                    is_alphafold=False,
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=8,
                    auth_chain="A",
                    label_chain="A",
                    method="Other",
                ),
                high_confidence_residues_count=None,
                geometry_quality=85.0,
                passed=True,
                reason=None,
                output_file=output_dir / "2Y29.cif",
            ),
            CombinedFilterResult(
                input_file=input_dir / "2y2a_updated.cif",
                pdb_id="2Y2A",
                metadata=StructureMetadata(
                    id="2Y2A",
                    uniprot_accession="P05067",
                    resolution=0.0,
                    total_residue_count=10,
                    is_alphafold=False,
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=10,
                    auth_chain="A",
                    label_chain="A",
                    method="Other",
                ),
                high_confidence_residues_count=None,
                geometry_quality=68.5,
                passed=False,
                reason="Geometry quality 68.5 below minimum 80.0",
                output_file=None,
            ),
        ]
        assert actual == expected

    def setup_uniprotless_resolution_less(self, no_uniprot_cif: Path, tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_file = input_dir / no_uniprot_cif.name
        strip_resolution(no_uniprot_cif, input_file)
        scores: dict[str, Scores] = {
            "2y29": Scores(
                geometry_quality=85.0,
                data_quality=89,
                experiment_data_available=False,
                overall_quality=85.0,
            ),
        }
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        expected = [
            CombinedFilterResult(
                input_file=input_file,
                pdb_id="2Y29",
                metadata=StructureMetadata(
                    id="2Y29",
                    uniprot_accession=None,
                    resolution=0.0,
                    total_residue_count=8,
                    is_alphafold=False,
                    uniprot_start=0,
                    uniprot_end=0,
                    sequence_identity=0.0,
                    chain_length=8,
                    auth_chain="A",
                    label_chain="A",
                    method="Other",
                ),
                high_confidence_residues_count=None,
                geometry_quality=85.0,
                passed=False,
                reason="Sorted index 1 > 0",
                output_file=None,
            )
        ]
        return input_file, scores, output_dir, expected

    def test_uniprotless_resolutionless_with_defaults(self, no_uniprot_cif: Path, tmp_path: Path):
        input_file, scores, output_dir, expected = self.setup_uniprotless_resolution_less(no_uniprot_cif, tmp_path)
        query = CombinedFilterQuery()

        actual = combined_filter([input_file], scores, query, output_dir, scheduler_address="sequential")

        assert actual == expected

    def test_uniprotless_resolutionless_with_scores_in_top(self, no_uniprot_cif: Path, tmp_path: Path):
        input_file, scores, output_dir, expected_base = self.setup_uniprotless_resolution_less(no_uniprot_cif, tmp_path)
        query = CombinedFilterQuery(
            top_non_uniprot=1,
        )

        actual = combined_filter([input_file], scores, query, output_dir, scheduler_address="sequential")

        expected = [
            replace(
                expected_base[0],
                passed=True,
                reason=None,
                output_file=output_dir / "no_uniprot.cif",
            )
        ]
        assert actual == expected

    def test_uniprotless_resolutionless_with_geometry_quality_too_low(self, no_uniprot_cif: Path, tmp_path: Path):
        input_file, scores, output_dir, expected_base = self.setup_uniprotless_resolution_less(no_uniprot_cif, tmp_path)
        query = CombinedFilterQuery(
            top_non_uniprot=1,
            min_geometry_quality=90.0,
        )

        actual = combined_filter([input_file], scores, query, output_dir, scheduler_address="sequential")

        expected = [
            replace(
                expected_base[0],
                passed=False,
                reason="Geometry quality 85.0 < 90.0",
                output_file=None,
            )
        ]
        assert actual == expected

    def test_sequence_identity(self, sample_multispan_cif: Path, multi_accession_chain_cif: Path, tmp_path: Path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        input_files = []
        for cif in [sample_multispan_cif, multi_accession_chain_cif]:
            input_file = input_dir / cif.name
            input_file.hardlink_to(cif)
            input_files.append(input_file)
        scores: dict[str, Scores] = {}
        query = CombinedFilterQuery(min_sequence_identity=0.9)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        results = combined_filter(input_files, scores, query, output_dir, scheduler_address="sequential")

        expected = [
            CombinedFilterResult(
                input_file=input_dir / "6O5I.cif.gz",
                pdb_id="6O5I",
                metadata=StructureMetadata(
                    id="6O5I",
                    uniprot_accession="O00255",
                    resolution=1.24,
                    total_residue_count=1346,
                    is_alphafold=False,
                    uniprot_start=1,
                    uniprot_end=593,
                    sequence_identity=0.816,  # is below 0.9
                    chain_length=1346,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=False,
                reason="Sequence identity 0.816 < 0.9",
                output_file=None,
            ),
            CombinedFilterResult(
                input_file=input_dir / "1un5.cif.gz",
                pdb_id="1UN5",
                metadata=StructureMetadata(
                    id="1UN5",
                    uniprot_accession="P03950",
                    resolution=2.6,
                    total_residue_count=131,
                    is_alphafold=False,
                    uniprot_start=25,
                    uniprot_end=147,
                    sequence_identity=0.967,  # is above 0.9
                    chain_length=131,
                    auth_chain="A",
                    label_chain="A",
                    method="X-ray",
                ),
                high_confidence_residues_count=None,
                geometry_quality=None,
                passed=True,
                reason=None,
                output_file=output_dir / "1un5.cif.gz",
            ),
        ]
        assert results == expected


@pytest.fixture
def fake_results() -> list[CombinedFilterResult]:
    return [
        CombinedFilterResult(
            input_file=Path("file1.cif"),
            pdb_id="1ABC",
            metadata=StructureMetadata(
                id="1ABC",
                uniprot_accession="P12345",
                resolution=2.0,
                total_residue_count=100,
                is_alphafold=False,
                uniprot_start=1,
                uniprot_end=100,
                sequence_identity=1.0,
                chain_length=100,
                auth_chain="A",
                label_chain="A",
                method="X-ray",
            ),
            high_confidence_residues_count=None,
            geometry_quality=0.5,
            passed=True,
            reason=None,
            output_file=Path("output/file1.cif"),
        ),
        CombinedFilterResult(
            input_file=Path("file2.cif"),
            pdb_id="2DEF",
            metadata=StructureMetadata(
                id="2DEF",
                uniprot_accession=None,
                resolution=0.0,
                total_residue_count=50,
                is_alphafold=False,
                uniprot_start=0,
                uniprot_end=0,
                sequence_identity=0.0,
                chain_length=50,
                auth_chain="B",
                label_chain="B",
                method="NMR",
            ),
            high_confidence_residues_count=None,
            geometry_quality=None,
            passed=False,
            reason="No UniProt accession, no resolution, no geometry quality",
            output_file=None,
        ),
        CombinedFilterResult(
            input_file=Path("file3.cif"),
            pdb_id="AF-3GHI",
            metadata=StructureMetadata(
                id="AF-3GHI",
                uniprot_accession="A0ABCD",
                resolution=0.0,
                total_residue_count=200,
                is_alphafold=True,
                uniprot_start=1,
                uniprot_end=200,
                sequence_identity=1.0,
                chain_length=200,
                auth_chain="A",
                label_chain="A",
                method="Predicted",
            ),
            high_confidence_residues_count=150,
            geometry_quality=None,
            passed=True,
            reason=None,
            output_file=Path("output/file3.cif"),
        ),
        CombinedFilterResult(
            input_file=Path("file4.cif"),
            pdb_id="4JKL",
            metadata=StructureMetadata(
                id="4JKL",
                uniprot_accession="P99999",
                resolution=4.5,
                total_residue_count=80,
                is_alphafold=False,
                uniprot_start=10,
                uniprot_end=90,
                sequence_identity=0.7,
                chain_length=80,
                auth_chain="C",
                label_chain="C",
                method="X-ray",
            ),
            high_confidence_residues_count=None,
            geometry_quality=25.0,
            passed=False,
            reason="Resolution 4.5 > 4.0",
            output_file=None,
        ),
        CombinedFilterResult(
            input_file=Path("file5.cif"),
            pdb_id="AF-5MNO",
            metadata=StructureMetadata(
                id="AF-5MNO",
                uniprot_accession="A1BCDE",
                resolution=0.0,
                total_residue_count=180,
                is_alphafold=True,
                uniprot_start=1,
                uniprot_end=180,
                sequence_identity=1.0,
                chain_length=180,
                auth_chain="A",
                label_chain="A",
                method="Predicted",
            ),
            high_confidence_residues_count=10,
            geometry_quality=15.0,
            passed=False,
            reason="High confidence residue ratio too low",
            output_file=None,
        ),
    ]


class TestCombinedFilterSummary:
    def test_zero(self):
        results = []

        actual = combined_filter_summary(results)

        expected = [
            "Total structures: 0, passed: 0, discarded: 0",
            "",
            "AlphaFold structures: 0, passed: 0, discarded: 0",
            "",
            "Structures with UniProt accession: 0, passed: 0, discarded: 0",
            "",
            "Structures with resolution: 0, passed: 0, discarded: 0",
            "",
            "Structures with geometry quality: 0, passed: 0, discarded: 0",
        ]
        assert actual == expected

    def test_nonzero(self, fake_results: list[CombinedFilterResult]):
        actual = combined_filter_summary(fake_results)

        expected = [
            "Total structures: 5, passed: 2, discarded: 3",
            "",
            "AlphaFold structures: 2, passed: 1, discarded: 1",
            "",
            "Structures with UniProt accession: 4, passed: 2, discarded: 2",
            "",
            "Structures with resolution: 2, passed: 1, discarded: 1",
            "",
            "Structures with geometry quality: 3, passed: 1, discarded: 2",
        ]
        assert actual == expected


def test_combined_filter_stats(fake_results: list[CombinedFilterResult], tmp_path: Path):
    path = StdioPath(tmp_path / "stats.csv")

    combined_filter_stats(fake_results, path)

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        actual = list(reader)
        expected = [
            {
                "chain_length": "100",
                "geometry_quality": "0.5",
                "high_confidence_residues_count": "",
                "input_file": "file1.cif",
                "is_alphafold": "False",
                "method": "X-ray",
                "output_file": "output/file1.cif",
                "passed": "True",
                "pdb_id": "1ABC",
                "reason": "",
                "resolution": "2.0",
                "sequence_identity": "1.0",
                "total_residue_count": "100",
                "uniprot_accession": "P12345",
                "uniprot_end": "100",
                "uniprot_start": "1",
            },
            {
                "chain_length": "50",
                "geometry_quality": "",
                "high_confidence_residues_count": "",
                "input_file": "file2.cif",
                "is_alphafold": "False",
                "method": "NMR",
                "output_file": "",
                "passed": "False",
                "pdb_id": "2DEF",
                "reason": "No UniProt accession, no resolution, no geometry quality",
                "resolution": "0.0",
                "sequence_identity": "0.0",
                "total_residue_count": "50",
                "uniprot_accession": "",
                "uniprot_end": "0",
                "uniprot_start": "0",
            },
            {
                "chain_length": "200",
                "geometry_quality": "",
                "high_confidence_residues_count": "150",
                "input_file": "file3.cif",
                "is_alphafold": "True",
                "method": "Predicted",
                "output_file": "output/file3.cif",
                "passed": "True",
                "pdb_id": "AF-3GHI",
                "reason": "",
                "resolution": "0.0",
                "sequence_identity": "1.0",
                "total_residue_count": "200",
                "uniprot_accession": "A0ABCD",
                "uniprot_end": "200",
                "uniprot_start": "1",
            },
            {
                "chain_length": "80",
                "geometry_quality": "25.0",
                "high_confidence_residues_count": "",
                "input_file": "file4.cif",
                "is_alphafold": "False",
                "method": "X-ray",
                "output_file": "",
                "passed": "False",
                "pdb_id": "4JKL",
                "reason": "Resolution 4.5 > 4.0",
                "resolution": "4.5",
                "sequence_identity": "0.7",
                "total_residue_count": "80",
                "uniprot_accession": "P99999",
                "uniprot_end": "90",
                "uniprot_start": "10",
            },
            {
                "chain_length": "180",
                "geometry_quality": "15.0",
                "high_confidence_residues_count": "10",
                "input_file": "file5.cif",
                "is_alphafold": "True",
                "method": "Predicted",
                "output_file": "",
                "passed": "False",
                "pdb_id": "AF-5MNO",
                "reason": "High confidence residue ratio too low",
                "resolution": "0.0",
                "sequence_identity": "1.0",
                "total_residue_count": "180",
                "uniprot_accession": "A1BCDE",
                "uniprot_end": "180",
                "uniprot_start": "1",
            },
        ]
        assert actual == expected

        expected_fieldnames = [
            "input_file",
            "pdb_id",
            "uniprot_accession",
            "resolution",
            "high_confidence_residues_count",
            "total_residue_count",
            "method",
            "is_alphafold",
            "uniprot_start",
            "uniprot_end",
            "sequence_identity",
            "chain_length",
            "geometry_quality",
            "passed",
            "output_file",
            "reason",
        ]
        assert reader.fieldnames == expected_fieldnames
