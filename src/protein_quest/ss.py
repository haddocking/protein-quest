"""Module for dealing with secondary structure."""

import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path

from gemmi import Structure, read_structure, set_leak_warnings

logger = logging.getLogger(__name__)

# TODO remove once v0.7.4 of gemmi is released,
# as uv pip install git+https://github.com/project-gemmi/gemmi.git installs 0.7.4.dev0 which does not print leaks
# Swallow gemmi leaked function warnings
set_leak_warnings(False)

# TODO if a structure has no secondary structure information, calculate it with `gemmi ss`.
# https://github.com/MonomerLibrary/monomers/wiki/Installation as --monomers dir
# gemmi executable is in https://pypi.org/project/gemmi-program/
# `gemmi ss` only prints secondary structure to stdout with `-v` flag.


def nr_of_residues_in_total(structure: Structure) -> int:
    """Count the total number of residues in the structure.

    Args:
        structure: The gemmi Structure object to analyze.

    Returns:
        The total number of residues in the structure.
    """
    count = 0
    for model in structure:
        for chain in model:
            count += len(chain)
    return count


def nr_of_residues_in_helix(structure: Structure) -> int:
    """Count the number of residues in alpha helices.

    Requires structure to have secondary structure information.

    Args:
        structure: The gemmi Structure object to analyze.

    Returns:
        The number of residues in alpha helices.
    """
    # For cif files from AlphaFold the helix.length is set to -1
    # so use resid instead
    count = 0
    for helix in structure.helices:
        end = helix.end.res_id.seqid.num
        start = helix.start.res_id.seqid.num
        if end is None or start is None:
            logger.warning(f"Invalid helix coordinates: {helix.end} or {helix.start}")
            continue
        length = end - start + 1
        count += length
    return count


def nr_of_residues_in_sheet(structure: Structure) -> int:
    """Count the number of residues in beta sheets.

    Requires structure to have secondary structure information.

    Args:
        structure: The gemmi Structure object to analyze.

    Returns:
        The number of residues in beta sheets.
    """
    count = 0
    for sheet in structure.sheets:
        for strand in sheet.strands:
            end = strand.end.res_id.seqid.num
            start = strand.start.res_id.seqid.num
            if end is None or start is None:
                logger.warning(f"Invalid strand coordinates: {strand.end} or {strand.start}")
                continue
            length = end - start + 1
            count += length
    return count


@dataclass
class SecondaryStructureFilterQuery:
    """Query object to filter on secondary structure.

    Parameters:
        abs_min_helix_residues: Minimum number of residues in helices (absolute).
        abs_max_helix_residues: Maximum number of residues in helices (absolute).
        abs_min_sheet_residues: Minimum number of residues in sheets (absolute).
        abs_max_sheet_residues: Maximum number of residues in sheets (absolute).
        ratio_min_helix_residues: Minimum number of residues in helices (relative).
        ratio_max_helix_residues: Maximum number of residues in helices (relative).
        ratio_min_sheet_residues: Minimum number of residues in sheets (relative).
        ratio_max_sheet_residues: Maximum number of residues in sheets (relative).
    """

    abs_min_helix_residues: int | None = None
    abs_max_helix_residues: int | None = None
    abs_min_sheet_residues: int | None = None
    abs_max_sheet_residues: int | None = None
    ratio_min_helix_residues: float | None = None
    ratio_max_helix_residues: float | None = None
    ratio_min_sheet_residues: float | None = None
    ratio_max_sheet_residues: float | None = None


@dataclass
class SecondaryStructureStats:
    nr_residues: int
    nr_helix_residues: int
    nr_sheet_residues: int
    helix_ratio: float
    sheet_ratio: float


@dataclass
class SecondaryStructureFilterResult:
    stats: SecondaryStructureStats
    passed: bool = False


def _gather_stats(structure: Structure) -> SecondaryStructureStats:
    nr_total_residues = nr_of_residues_in_total(structure)
    nr_helix_residues = nr_of_residues_in_helix(structure)
    nr_sheet_residues = nr_of_residues_in_sheet(structure)
    if nr_total_residues == 0:
        msg = "Structure has zero residues; cannot compute secondary structure ratios."
        raise ValueError(msg)
    helix_ratio = nr_helix_residues / nr_total_residues
    sheet_ratio = nr_sheet_residues / nr_total_residues
    return SecondaryStructureStats(
        nr_residues=nr_total_residues,
        nr_helix_residues=nr_helix_residues,
        nr_sheet_residues=nr_sheet_residues,
        helix_ratio=helix_ratio,
        sheet_ratio=sheet_ratio,
    )


def filter_on_secondary_structure(
    structure: Structure,
    query: SecondaryStructureFilterQuery,
) -> SecondaryStructureFilterResult:
    """Filter a structure based on secondary structure criteria.

    Args:
        structure: The gemmi Structure object to analyze.
        query: The filtering criteria to apply.

    Returns:
        Filtering statistics and whether structure passed.
    """
    stats = _gather_stats(structure)
    conditions: list[bool] = []

    # Helix absolute thresholds
    if query.abs_min_helix_residues is not None:
        conditions.append(stats.nr_helix_residues >= query.abs_min_helix_residues)
    if query.abs_max_helix_residues is not None:
        conditions.append(stats.nr_helix_residues <= query.abs_max_helix_residues)

    # Helix ratio thresholds
    if query.ratio_min_helix_residues is not None:
        conditions.append(stats.helix_ratio >= query.ratio_min_helix_residues)
    if query.ratio_max_helix_residues is not None:
        conditions.append(stats.helix_ratio <= query.ratio_max_helix_residues)

    # Sheet absolute thresholds
    if query.abs_min_sheet_residues is not None:
        conditions.append(stats.nr_sheet_residues >= query.abs_min_sheet_residues)
    if query.abs_max_sheet_residues is not None:
        conditions.append(stats.nr_sheet_residues <= query.abs_max_sheet_residues)

    # Sheet ratio thresholds
    if query.ratio_min_sheet_residues is not None:
        conditions.append(stats.sheet_ratio >= query.ratio_min_sheet_residues)
    if query.ratio_max_sheet_residues is not None:
        conditions.append(stats.sheet_ratio <= query.ratio_max_sheet_residues)

    if not conditions:
        msg = "No filtering conditions provided. Please specify at least one condition."
        raise ValueError(msg)
    passed = all(conditions)
    return SecondaryStructureFilterResult(stats=stats, passed=passed)


def filter_file_on_secondary_structure(
    file_path: Path,
    query: SecondaryStructureFilterQuery,
) -> SecondaryStructureFilterResult:
    """Filter a structure file based on secondary structure criteria.

    Args:
        file_path: The path to the structure file to analyze.
        query: The filtering criteria to apply.

    Returns:
        Filtering statistics and whether file passed.
    """
    structure = read_structure(str(file_path))
    return filter_on_secondary_structure(structure, query)


def filter_files_on_secondary_structure(
    file_paths: Iterable[Path],
    query: SecondaryStructureFilterQuery,
) -> Generator[tuple[Path, SecondaryStructureFilterResult]]:
    """Filter multiple structure files based on secondary structure criteria.

    Args:
        file_paths: A list of paths to the structure files to analyze.
        query: The filtering criteria to apply.

    Yields:
        For each file returns the filtering statistics and whether structure passed.
    """
    # TODO check if quick enough in serial mode, if not switch to dask map
    for file_path in file_paths:
        result = filter_file_on_secondary_structure(file_path, query)
        yield file_path, result
