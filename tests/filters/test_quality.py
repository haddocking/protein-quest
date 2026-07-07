import csv
from pathlib import Path

import pytest
from cyclopts.types import StdioPath

from protein_quest.filters.quality import (
    FilterQualityResult,
    QualityClusteringPartitions,
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


def test_partition_structures_pass_given_resolution(sample_cif: Path):
    """Structures with valid resolution bypass quality checks when pass_given_resolution=True."""
    # sample_cif is 3JRS with resolution=2.05; low geometry_quality should not matter
    scores = {
        "3jrs": Scores(
            geometry_quality=1.0,
            data_quality=None,
            overall_quality=None,
            experiment_data_available=False,
        )
    }

    partitions = partition_structures_for_quality_clustering(
        [sample_cif],
        scores,
        pass_given_resolution=True,
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
            reason="",
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
