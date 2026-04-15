"""Filter structure files by number of residues."""

from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path

from tqdm.auto import tqdm

from protein_quest.structure import nr_residues_in_chain
from protein_quest.utils import CopyMethod, copyfile


@dataclass
class ResidueFilterStatistics:
    """Statistics for filtering files based on residue count in a specific chain.

    Parameters:
        input_file: The path to the input file.
        residue_count: The number of residues.
        passed: Whether the file passed the filtering criteria.
        output_file: The path to the output file, if passed.
    """

    input_file: Path
    residue_count: int
    passed: bool
    output_file: Path | None


def filter_files_on_residues(
    input_files: list[Path],
    output_dir: Path,
    min_residues: int,
    max_residues: int,
    chain: str = "A",
    copy_method: CopyMethod = "copy",
) -> Generator[ResidueFilterStatistics]:
    """Filter PDB/mmCIF files by number of residues in given chain.

    Args:
        input_files: The list of input PDB/mmCIF files.
        output_dir: The directory where the filtered files will be written.
        min_residues: The minimum number of residues in chain.
        max_residues: The maximum number of residues in chain.
        chain: The chain to count residues of.
        copy_method: The method used to copy passed files to the output directory.

    Yields:
        Objects containing information about the filtering process for each input file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for input_file in tqdm(input_files, unit="file"):
        residue_count = nr_residues_in_chain(input_file, chain=chain)
        passed = min_residues <= residue_count <= max_residues
        if passed:
            output_file = output_dir / input_file.name
            copyfile(input_file, output_file, copy_method)
            yield ResidueFilterStatistics(input_file, residue_count, True, output_file)
        else:
            yield ResidueFilterStatistics(input_file, residue_count, False, None)
