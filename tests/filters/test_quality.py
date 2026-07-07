import csv
from pathlib import Path

import pytest
from cyclopts.types import StdioPath

from protein_quest.filters.quality import (
    FilterQualityResult,
    QualityClusteringPartitions,
    QualityStructure,
    UnclusteredStructure,
    filter_unclustered_structures,
    partition_structures_for_quality_clustering,
    write_quality_stats_csv,
)
from protein_quest.pdbe.ws import Scores


@pytest.mark.parametrize(
    "results, expected",
    [
        pytest.param(
            [],
            [],
            id="empty",
        ),
        pytest.param(
            [
                FilterQualityResult(
                    pdb_id="1AAA",
                    input_file=Path("/fake/input/1AAA.pdb"),
                    geometry_quality=0.8,
                    passed=True,
                )
            ],
            [
                {
                    "pdb_id": "1AAA",
                    "input_file": "/fake/input/1AAA.pdb",
                    "geometry_quality": "0.8",
                    "passed": "True",
                    "output_file": "/fake/output/1AAA.pdb",
                    "reason": "",
                }
            ],
            id="passed",
        ),
        pytest.param(
            [
                FilterQualityResult(
                    pdb_id="1AAA",
                    input_file=Path("/fake/input/1AAA.pdb"),
                    geometry_quality=0.3,
                    passed=False,
                    reason="Geometry quality score 0.3 < 0.5",
                )
            ],
            [
                {
                    "pdb_id": "1AAA",
                    "input_file": "/fake/input/1AAA.pdb",
                    "geometry_quality": "0.3",
                    "passed": "False",
                    "output_file": "",
                    "reason": "Geometry quality score 0.3 < 0.5",
                }
            ],
            id="gq_to_low",
        ),
    ],
)
def test_write_quality_stats_csv(
    tmp_path: Path,
    results: list[FilterQualityResult],
    expected: list[dict[str, str]],
):
    output_file = StdioPath(tmp_path / "quality_stats.csv")

    write_quality_stats_csv(results, output_file, Path("/fake/output"))

    assert output_file.exists()

    with output_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert rows == expected


class TestPartitionStructuresForQualityClustering:
    @staticmethod
    def _scores() -> dict[str, Scores]:
        return {
            "1un5": Scores(
                geometry_quality=45.23,
                data_quality=None,
                overall_quality=45.23,
                experiment_data_available=False,
            ),
            "2y29": Scores(
                geometry_quality=55.9,
                data_quality=81.09,
                overall_quality=63.83,
                experiment_data_available=True,
            ),
            "3jrs": Scores(
                geometry_quality=31.58,
                data_quality=22.4,
                overall_quality=27.13,
                experiment_data_available=True,
            ),
            "1amb": Scores(
                geometry_quality=None,
                data_quality=None,
                overall_quality=None,
                experiment_data_available=False,
            ),
            "6o5i": Scores(
                geometry_quality=64.38,
                data_quality=61.45,
                overall_quality=63.17,
                experiment_data_available=True,
            ),
            "8w77": Scores(
                geometry_quality=75.0,
                data_quality=75.0,
                overall_quality=75.0,
                experiment_data_available=True,
            ),
        }

    def test_pass_given_resolution(self, sample_cif: Path):
        """Structures with valid resolution bypass quality checks when pass_given_resolution=True."""
        partitions = partition_structures_for_quality_clustering(
            [sample_cif],
            {
                "3jrs": Scores(
                    geometry_quality=1.0,
                    data_quality=None,
                    overall_quality=None,
                    experiment_data_available=False,
                )
            },
            pass_given_resolution=True,
            scheduler_address="sequential",
        )

        expected = QualityClusteringPartitions(
            clusterable_structures=[],
            unclustered_structures=[],
            no_quality_results=[],
            resolution_passed_results=[
                FilterQualityResult(
                    pdb_id="3JRS",
                    input_file=sample_cif,
                    geometry_quality=1.0,
                    passed=True,
                    reason="Passed due to valid resolution 2.05",
                )
            ],
        )
        assert partitions == expected

    @pytest.mark.parametrize(
        "scheduler_address",
        [
            pytest.param("sequential", id="sequential"),
            pytest.param(None, id="default-parallel"),
        ],
    )
    def test_scheduler_address(
        self,
        sample_cif: Path,
        sample2_cif: Path,
        sample_multispan_cif: Path,
        multi_accession_chain_cif: Path,
        af_cif: Path,
        nmr_cif: Path,
        em_cif: Path,
        no_uniprot_cif: Path,
        unreadable_cif: Path,
        atomless_cif: Path,
        scheduler_address: str | None,
    ):
        scores = self._scores()

        partitions = partition_structures_for_quality_clustering(
            [
                sample_cif,
                sample2_cif,
                sample_multispan_cif,
                multi_accession_chain_cif,
                af_cif,
                nmr_cif,
                em_cif,
                no_uniprot_cif,
                unreadable_cif,
                atomless_cif,
            ],
            scores,
            scheduler_address=scheduler_address,
        )

        expected = QualityClusteringPartitions(
            clusterable_structures=[
                QualityStructure(
                    id="3JRS",
                    uniprot_accession="Q8VZS8",
                    uniprot_start=8,
                    uniprot_end=211,
                    sequence_identity=1.0,
                    chain_length=173,
                    geometry_quality=31.58,
                    input_file=sample_cif,
                ),
                QualityStructure(
                    id="2Y29",
                    uniprot_accession="P05067",
                    uniprot_start=687,
                    uniprot_end=692,
                    sequence_identity=1.0,
                    chain_length=8,
                    geometry_quality=55.9,
                    input_file=sample2_cif,
                ),
                QualityStructure(
                    id="6O5I",
                    uniprot_accession="O00255",
                    uniprot_start=1,
                    uniprot_end=593,
                    sequence_identity=0.816,
                    chain_length=1346,
                    geometry_quality=64.38,
                    input_file=sample_multispan_cif,
                ),
                QualityStructure(
                    id="1UN5",
                    uniprot_accession="P03950",
                    uniprot_start=25,
                    uniprot_end=147,
                    sequence_identity=0.967,
                    chain_length=131,
                    geometry_quality=45.23,
                    input_file=multi_accession_chain_cif,
                ),
                QualityStructure(
                    id="8W77",
                    uniprot_accession="P0ABE7",
                    uniprot_start=23,
                    uniprot_end=127,
                    sequence_identity=1.0,
                    chain_length=260,
                    geometry_quality=75.0,
                    input_file=em_cif,
                ),
            ],
            unclustered_structures=[
                UnclusteredStructure(
                    input_file=no_uniprot_cif,
                    pdb_id="2Y29",
                    geometry_quality=55.9,
                    chain_length=8,
                    sequence_identity=0.0,
                )
            ],
            no_quality_results=[
                FilterQualityResult(
                    pdb_id="1AMB",
                    input_file=nmr_cif,
                    geometry_quality=None,
                    passed=False,
                    reason="Missing geometry quality",
                ),
                FilterQualityResult(
                    pdb_id=None,
                    input_file=unreadable_cif,
                    geometry_quality=None,
                    passed=False,
                    reason=f"Failed to read structure: {unreadable_cif}:1:0(0): expected block header (data_)",
                ),
                FilterQualityResult(
                    pdb_id=" ",
                    input_file=atomless_cif,
                    geometry_quality=None,
                    passed=False,
                    reason="Failed to extract metadata: No chains found in structure  ",
                ),
            ],
            resolution_passed_results=[
                FilterQualityResult(
                    pdb_id="AF-A0A0C5B5G6-F1",
                    input_file=af_cif,
                    geometry_quality=None,
                    passed=True,
                    reason="AlphaFold structure passes quality filter",
                )
            ],
        )

        assert partitions == expected


def test_filter_unclustered_structures_passed_but_outside_top():
    """Structures that pass quality but exceed the top limit are rejected."""
    structures = [
        UnclusteredStructure(
            input_file=Path("/a.cif"),
            pdb_id="1AAA",
            geometry_quality=80.0,
            chain_length=100,
            sequence_identity=1.0,
        ),
        UnclusteredStructure(
            input_file=Path("/b.cif"),
            pdb_id="2BBB",
            geometry_quality=70.0,
            chain_length=100,
            sequence_identity=1.0,
        ),
    ]

    results = filter_unclustered_structures(
        structures,
        minimal_geometry_quality=50.0,
        top=1,
    )

    expected = [
        FilterQualityResult(
            pdb_id="1AAA",
            input_file=Path("/a.cif"),
            geometry_quality=80.0,
            passed=True,
            reason=None,
        ),
        FilterQualityResult(
            pdb_id="2BBB",
            input_file=Path("/b.cif"),
            geometry_quality=70.0,
            passed=False,
            reason="Excluded by top 1 limit for unclustered structures",
        ),
    ]
    assert results == expected
