"""Module for dealing with secondary structure."""

import logging
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from pathlib import Path

from gemmi import Structure

from protein_quest.converter import PositiveInt, Ratio, converter
from protein_quest.io import read_structure

logger = logging.getLogger(__name__)

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

    abs_min_helix_residues: PositiveInt | None = None
    abs_max_helix_residues: PositiveInt | None = None
    abs_min_sheet_residues: PositiveInt | None = None
    abs_max_sheet_residues: PositiveInt | None = None
    ratio_min_helix_residues: Ratio | None = None
    ratio_max_helix_residues: Ratio | None = None
    ratio_min_sheet_residues: Ratio | None = None
    ratio_max_sheet_residues: Ratio | None = None

    def is_actionable(self) -> bool:
        """Check if the secondary structure query has any actionable filters.

        Returns:
            True if any of the filters are set, False otherwise.
        """
        return any(
            field is not None
            for field in [
                self.abs_min_helix_residues,
                self.abs_max_helix_residues,
                self.abs_min_sheet_residues,
                self.abs_max_sheet_residues,
                self.ratio_min_helix_residues,
                self.ratio_max_helix_residues,
                self.ratio_min_sheet_residues,
                self.ratio_max_sheet_residues,
            ]
        )


def _check_range(min_val, max_val, label):
    if min_val is not None and max_val is not None and min_val >= max_val:
        msg = f"Invalid {label} range: min {min_val} must be smaller than max {max_val}"
        raise ValueError(msg)


base_query_hook = converter.get_structure_hook(SecondaryStructureFilterQuery)


@converter.register_structure_hook
def secondary_structure_filter_query_hook(value, _type) -> SecondaryStructureFilterQuery:
    result: SecondaryStructureFilterQuery = base_query_hook(value, _type)
    _check_range(result.abs_min_helix_residues, result.abs_max_helix_residues, "absolute helix residue")
    _check_range(result.abs_min_sheet_residues, result.abs_max_sheet_residues, "absolute sheet residue")
    _check_range(result.ratio_min_helix_residues, result.ratio_max_helix_residues, "ratio helix residue")
    _check_range(result.ratio_min_sheet_residues, result.ratio_max_sheet_residues, "ratio sheet residue")
    return result


@dataclass
class SecondaryStructureStats:
    """Statistics about the secondary structure of a protein.

    Parameters:
        nr_residues: Total number of residues in the structure.
        nr_helix_residues: Number of residues in helices.
        nr_sheet_residues: Number of residues in sheets.
        helix_ratio: Ratio of residues in helices.
        sheet_ratio: Ratio of residues in sheets.
    """

    nr_residues: PositiveInt
    nr_helix_residues: PositiveInt
    nr_sheet_residues: PositiveInt
    helix_ratio: Ratio
    sheet_ratio: Ratio


@dataclass
class SecondaryStructureFilterResult:
    """Result of filtering on secondary structure.

    Parameters:
        stats: The secondary structure statistics.
        passed: Whether the structure passed the filtering criteria.
    """

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
    structure = read_structure(file_path)
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
