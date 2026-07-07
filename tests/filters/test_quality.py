import csv
from pathlib import Path

import pytest
from cyclopts.types import StdioPath

from protein_quest.filters.quality import FilterQualityResult, write_quality_stats_csv


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
