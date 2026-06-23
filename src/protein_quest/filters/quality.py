import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

from cyclopts.types import StdioPath
from tqdm.rich import tqdm

from protein_quest.pdbe.ws import Scores
from protein_quest.structure.files import LocateStructureFilesByIdResult
from protein_quest.structure.formats import read_structure

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class FilterQualityResult:
    """Result of filtering a structure file based on PDBe quality scores.

    Attributes:
        pdb_id: The PDB ID of the structure if available.
        input_file: The path to the input structure file if available.
        geometry_quality: The geometry quality score if available.
        passed: A boolean indicating whether the structure passed the filter.
        reason: The reason for discarding or passing the structure.
            None if passed without any special reason.
    """

    pdb_id: str | None = None
    input_file: Path | None = None
    geometry_quality: float | None = None
    passed: bool = False
    reason: str | None = None


class PdbIdGeometryQualityPair(TypedDict):
    pdb_id: str
    geometry_quality: float | None


def _associate_files_with_sorted_scores(
    scores: dict[str, Scores], located_ids: LocateStructureFilesByIdResult, top: int | None = None
) -> tuple[list[PdbIdGeometryQualityPair], dict[str, list[Path]]]:
    flattened_scores: list[PdbIdGeometryQualityPair] = []
    found_ids = {pdb_id for pdb_id, _ in located_ids.found}
    for pdb_id, score in scores.items():
        if pdb_id not in found_ids:
            continue
        flattened_scores.append(
            PdbIdGeometryQualityPair(
                pdb_id=pdb_id,
                geometry_quality=score.geometry_quality,
            )
        )

    def sorter(x: PdbIdGeometryQualityPair) -> tuple[bool, float | None]:
        # None as worst and highest as best with reverse sort
        return (x["geometry_quality"] is not None, x["geometry_quality"])

    sorted_scores = sorted(flattened_scores, key=sorter, reverse=True)
    top_sorted_scores = sorted_scores[:top] if top is not None else sorted_scores

    # Convert (('a', file_a), ('b', file_b), ('a', file_c), ('b', file_c)) to {'a': [file_a, file_c], 'b': [file_b]}
    found_id2files: dict[str, list[Path]] = {}
    found_files2ids: dict[Path, str] = {}
    for pdb_id, file in located_ids.found:
        if file in found_files2ids:
            logger.info(
                f"File {file} already associated with PDB ID {found_files2ids[file]}, "
                f"skipping association with PDB ID {pdb_id}"
            )
            continue
        found_id2files.setdefault(pdb_id, []).append(file)
        found_files2ids[file] = pdb_id
    return top_sorted_scores, found_id2files


def filter_by_pdbe_quality(
    scores: dict[str, Scores],
    located_ids: LocateStructureFilesByIdResult,
    /,
    *,
    minimal_geometry_quality: float = 0.0,
    top: int | None = None,
    pass_given_resolution: bool = False,
) -> list[FilterQualityResult]:
    """Filter structure files based on PDBe quality scores.

    Args:
        scores: A dictionary mapping PDB IDs to their corresponding Scores objects.
        located_ids: A dictionary containing located structure files by PDB ID.
        minimal_geometry_quality: Minimum geometry quality score to pass the filter.
        top: If set, only consider the top N structures based on geometry quality.
        pass_given_resolution: If set, structures with a valid resolution will pass regardless of other criteria.

    Returns:
        A list of QualityResult objects.
        No files are written to disk; the caller is responsible for copying or moving the files as needed.
    """
    top_sorted_scores, found_id2files = _associate_files_with_sorted_scores(scores, located_ids, top)

    results: list[FilterQualityResult] = []
    for score_dict in tqdm(top_sorted_scores, desc="Filtering on PDBe quality", unit="file"):
        pdb_id = score_dict["pdb_id"]
        score = scores[pdb_id]
        input_files = found_id2files[pdb_id]
        for input_file in input_files:
            if pass_given_resolution and (
                # Only read structure when needed as checking score is cheaper than reading structure
                score.geometry_quality is None or score.geometry_quality < minimal_geometry_quality
            ):
                structure = read_structure(input_file)
                if structure.resolution != 0.0:
                    results.append(
                        FilterQualityResult(
                            pdb_id=pdb_id,
                            input_file=input_file,
                            geometry_quality=score.geometry_quality,
                            passed=True,
                            reason=f"Passed due to valid resolution {structure.resolution}",
                        )
                    )
                    logger.debug(f"Passing {pdb_id} due to valid resolution {structure.resolution}")
                    continue
            if score.geometry_quality is None:
                results.append(
                    FilterQualityResult(
                        pdb_id=pdb_id,
                        input_file=input_file,
                        reason="No geometry quality score",
                    )
                )
                continue
            if score.geometry_quality >= minimal_geometry_quality:
                results.append(
                    FilterQualityResult(
                        pdb_id=pdb_id,
                        input_file=input_file,
                        geometry_quality=score.geometry_quality,
                        passed=True,
                    )
                )
            else:
                results.append(
                    FilterQualityResult(
                        pdb_id=pdb_id,
                        input_file=input_file,
                        geometry_quality=score.geometry_quality,
                        reason=f"Geometry quality score {score.geometry_quality} < {minimal_geometry_quality}",
                    )
                )

    results.extend(
        FilterQualityResult(
            pdb_id=pdb_id,
            geometry_quality=scores[pdb_id].geometry_quality,
            reason="File not found",
        )
        for pdb_id in located_ids.not_found
    )
    results.extend(
        FilterQualityResult(
            input_file=extra,
            reason="File not found in quality scores",
        )
        for extra in located_ids.extras
    )

    return results


def write_quality_stats_csv(
    results: list[FilterQualityResult],
    write_stats: StdioPath,
    output_dir: Path,
):
    """Writes a CSV file containing quality statistics for PDB IDs.

    Args:
        results: A list of QualityResult objects containing the quality results for each PDB ID.
        write_stats: Path to the output CSV file to write statistics.
        output_dir: The directory where the output files are located.

    The CSV file will contain the following columns:
        - pdb_id: The PDB ID of the structure.
        - input_file: The path to the input structure file.
        - geometry_quality: The geometry quality score from the Scores object.
        - passed: A boolean indicating whether the structure passed the filter.
        - output_file: The path to the output structure file if it passed the filter.
        - reason: The reason for discarding or passing the structure.
    """
    with write_stats.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["pdb_id", "input_file", "geometry_quality", "passed", "output_file", "reason"]
        )
        writer.writeheader()
        for result in results:
            output_file = ""
            if result.passed and result.input_file is not None:
                output_file = output_dir / result.input_file.name

            writer.writerow(
                {
                    "pdb_id": result.pdb_id,
                    "input_file": result.input_file,
                    "geometry_quality": result.geometry_quality,
                    "passed": result.passed,
                    "output_file": output_file,
                    "reason": result.reason,
                }
            )
