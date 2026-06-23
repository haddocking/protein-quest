import csv
from pathlib import Path

import pytest
from cyclopts.types import StdioPath

from protein_quest.filters.quality import FilterQualityResult, filter_by_pdbe_quality, write_quality_stats_csv
from protein_quest.pdbe.ws import Scores
from protein_quest.structure.files import LocateStructureFilesByIdResult


class TestFilterByPdbeQuality:
    @pytest.mark.parametrize(
        "scores, located_ids, top, expected",
        [
            pytest.param({}, LocateStructureFilesByIdResult(), None, [], id="empty"),
            pytest.param(
                {
                    "1AAA": Scores(
                        geometry_quality=80.0, data_quality=90.0, overall_quality=85.0, experiment_data_available=True
                    )
                },
                LocateStructureFilesByIdResult(found={("1AAA", Path("/fake/input/1AAA.pdb"))}),
                None,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=Path("/fake/input/1AAA.pdb"),
                        geometry_quality=80.0,
                        passed=True,
                    )
                ],
                id="passed",
            ),
            pytest.param(
                {
                    "1AAA": Scores(
                        geometry_quality=80.0, data_quality=90.0, overall_quality=85.0, experiment_data_available=True
                    )
                },
                LocateStructureFilesByIdResult(not_found={"1AAA"}),
                None,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=None,
                        geometry_quality=80.0,
                        passed=False,
                        reason="File not found",
                    )
                ],
                id="not_found_on_disk",
            ),
            pytest.param(
                {},
                LocateStructureFilesByIdResult(
                    extras={Path("/fake/input/1AAA.pdb")},
                ),
                None,
                [
                    FilterQualityResult(
                        input_file=Path("/fake/input/1AAA.pdb"),
                        passed=False,
                        reason="File not found in quality scores",
                    )
                ],
                id="not_found_in_scores",
            ),
            pytest.param(
                {
                    "1AAA": Scores(
                        geometry_quality=None, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    )
                },
                LocateStructureFilesByIdResult(
                    found={("1AAA", Path("/fake/input/1AAA.pdb"))},
                ),
                None,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=Path("/fake/input/1AAA.pdb"),
                        geometry_quality=None,
                        passed=False,
                        reason="No geometry quality score",
                    )
                ],
                id="no_geometry_quality",
            ),
            pytest.param(
                {
                    "1AAA": Scores(
                        geometry_quality=10.0, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    )
                },
                LocateStructureFilesByIdResult(
                    found={("1AAA", Path("/fake/input/1AAA.pdb"))},
                ),
                None,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=Path("/fake/input/1AAA.pdb"),
                        geometry_quality=10.0,
                        passed=False,
                        reason="Geometry quality score 10.0 < 50.0",
                    )
                ],
                id="gq_to_low",
            ),
            pytest.param(
                {
                    "3CCC": Scores(
                        geometry_quality=90.0, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    ),
                    "2BBB": Scores(
                        geometry_quality=80.0, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    ),
                    "1AAA": Scores(
                        geometry_quality=100.0, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    ),
                },
                LocateStructureFilesByIdResult(
                    found={
                        ("1AAA", Path("/fake/input/1AAA.pdb")),
                        ("2BBB", Path("/fake/input/2BBB.pdb")),
                        ("3CCC", Path("/fake/input/3CCC.pdb")),
                    },
                ),
                None,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=Path("/fake/input/1AAA.pdb"),
                        geometry_quality=100.0,
                        passed=True,
                    ),
                    FilterQualityResult(
                        pdb_id="3CCC",
                        input_file=Path("/fake/input/3CCC.pdb"),
                        geometry_quality=90.0,
                        passed=True,
                    ),
                    FilterQualityResult(
                        pdb_id="2BBB",
                        input_file=Path("/fake/input/2BBB.pdb"),
                        geometry_quality=80.0,
                        passed=True,
                    ),
                ],
                id="sorted_by_quality",
            ),
            pytest.param(
                {
                    "3CCC": Scores(
                        geometry_quality=90, data_quality=0.6, overall_quality=0.55, experiment_data_available=True
                    ),
                    "2BBB": Scores(
                        geometry_quality=80, data_quality=0.6, overall_quality=0.55, experiment_data_available=True
                    ),
                    "1AAA": Scores(
                        geometry_quality=100, data_quality=0.6, overall_quality=0.55, experiment_data_available=True
                    ),
                },
                LocateStructureFilesByIdResult(
                    found={
                        ("1AAA", Path("/fake/input/1AAA.pdb")),
                        ("2BBB", Path("/fake/input/2BBB.pdb")),
                        ("3CCC", Path("/fake/input/3CCC.pdb")),
                    },
                ),
                2,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=Path("/fake/input/1AAA.pdb"),
                        geometry_quality=100,
                        passed=True,
                    ),
                    FilterQualityResult(
                        pdb_id="3CCC",
                        input_file=Path("/fake/input/3CCC.pdb"),
                        geometry_quality=90,
                        passed=True,
                    ),
                ],
                id="top2",
            ),
            pytest.param(
                {
                    "3CCC": Scores(
                        geometry_quality=90.0, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    ),
                    "2BBB": Scores(
                        geometry_quality=None, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    ),
                    "1AAA": Scores(
                        geometry_quality=100.0, data_quality=60.0, overall_quality=55.0, experiment_data_available=True
                    ),
                },
                LocateStructureFilesByIdResult(
                    found={
                        ("1AAA", Path("/fake/input/1AAA.pdb")),
                        ("2BBB", Path("/fake/input/2BBB.pdb")),
                        ("3CCC", Path("/fake/input/3CCC.pdb")),
                    },
                ),
                None,
                [
                    FilterQualityResult(
                        pdb_id="1AAA",
                        input_file=Path("/fake/input/1AAA.pdb"),
                        geometry_quality=100.0,
                        passed=True,
                    ),
                    FilterQualityResult(
                        pdb_id="3CCC",
                        input_file=Path("/fake/input/3CCC.pdb"),
                        geometry_quality=90.0,
                        passed=True,
                    ),
                    FilterQualityResult(
                        pdb_id="2BBB",
                        input_file=Path("/fake/input/2BBB.pdb"),
                        geometry_quality=None,
                        passed=False,
                        reason="No geometry quality score",
                    ),
                ],
                id="no_geometry_quality_worst",
            ),
        ],
    )
    def test(
        self,
        scores: dict[str, Scores],
        located_ids: LocateStructureFilesByIdResult,
        top: int | None,
        expected: list[FilterQualityResult],
    ):
        result = filter_by_pdbe_quality(
            scores,
            located_ids,
            minimal_geometry_quality=50.0,
            top=top,
        )
        assert result == expected

    def test_bypass_with_set_resolution(self, tmp_path: Path, sample_cif: Path):
        scores = {
            "1AAA": Scores(geometry_quality=0.1, data_quality=0.6, overall_quality=0.55, experiment_data_available=True)
        }
        located_ids = LocateStructureFilesByIdResult(
            found={("1AAA", sample_cif)},
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        result = filter_by_pdbe_quality(
            scores,
            located_ids,
            minimal_geometry_quality=0.5,
            pass_given_resolution=True,
        )

        expected = [
            FilterQualityResult(
                pdb_id="1AAA",
                input_file=sample_cif,
                geometry_quality=0.1,
                passed=True,
                reason="Passed due to valid resolution 2.05",
            )
        ]
        assert result == expected

    def test_bypass_with_unset_resolution(self, tmp_path: Path, nmr_cif: Path):
        scores = {
            "1AAA": Scores(geometry_quality=0.1, data_quality=0.6, overall_quality=0.55, experiment_data_available=True)
        }
        located_ids = LocateStructureFilesByIdResult(
            found={("1AAA", nmr_cif)},
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        result = filter_by_pdbe_quality(
            scores,
            located_ids,
            minimal_geometry_quality=0.5,
            pass_given_resolution=True,
        )

        expected = [
            FilterQualityResult(
                pdb_id="1AAA",
                input_file=nmr_cif,
                geometry_quality=0.1,
                passed=False,
                reason="Geometry quality score 0.1 < 0.5",
            )
        ]
        assert result == expected


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
