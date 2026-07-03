"""Filter subcommands for protein-quest."""

import csv
import logging
import os
from typing import TYPE_CHECKING, Annotated

from cyclopts import App, Parameter
from cyclopts.types import NormFloat, PositiveInt
from rich.panel import Panel

from protein_quest.alphafold.confidence import ConfidenceFilterQuery, filter_files_on_confidence
from protein_quest.cli.common import (
    CacheParameter,
    Common,
    InputDir,
    InputFile,
    MaxResidues,
    MinResidues,
    OutputDir,
    OutputFile,
    console,
    write_lines,
)
from protein_quest.converter import converter
from protein_quest.filters.chain import filter_files_on_chain
from protein_quest.filters.quality import (
    filter_by_pdbe_quality,
    filter_by_pdbe_quality_clustered,
    write_quality_stats_csv,
)
from protein_quest.filters.residues import filter_files_on_residues
from protein_quest.filters.resolution import (
    filter_files_on_resolution,
    write_resolution_stats,
)
from protein_quest.filters.ss import SecondaryStructureFilterQuery, filter_files_on_secondary_structure
from protein_quest.pdbe.ws import Scores
from protein_quest.structure.chains import ChainIdSystem
from protein_quest.structure.files import glob_structure_files, locate_structure_file, locate_structure_files_by_id
from protein_quest.utils import copyfile

if TYPE_CHECKING:
    from pathlib import Path

rprint = console.print
logger = logging.getLogger(__name__)

filter_app = App(name="filter", help="Filter files")


@filter_app.command
def confidence(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    filters: ConfidenceFilterQuery | None = None,
    write_stats: OutputFile | None = None,
    scheduler_address: str | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter AlphaFold mmcif/PDB files by confidence (plDDT).

    Filter AlphaFold mmcif/PDB files by confidence (plDDT).
    Passed files are written with residues below threshold removed.

    Args:
        input_dir: Directory with AlphaFold mmcif/PDB files.
        output_dir: Directory to write filtered mmcif/PDB files.
        filters: Confidence filtering criteria.
        write_stats: Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    if filters is None:
        filters = ConfidenceFilterQuery()
    output_dir.mkdir(parents=True, exist_ok=True)
    input_files = sorted(glob_structure_files(input_dir))
    nr_input_files = len(input_files)
    rprint(f"Starting confidence filtering of {nr_input_files} mmcif/PDB files in {input_dir} directory.")

    if write_stats and str(write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    cache = cache or CacheParameter()

    passed_count = 0
    results = filter_files_on_confidence(
        input_files, filters, output_dir, copy_method=cache.copy_method, scheduler_address=scheduler_address
    )

    stats_lines = ["input_file,residue_count,passed,output_file"]
    for r in results:
        if r.filtered_file:
            passed_count += 1
        stats_lines.append(f"{r.input_file},{r.count},{r.filtered_file is not None},{r.filtered_file}")

    if write_stats:
        write_lines(write_stats, stats_lines)
    rprint(f"Filtered {passed_count} mmcif/PDB files by confidence, written to {output_dir} directory")
    if str(write_stats) != "-":
        rprint(f"Statistics written to {write_stats}")


@filter_app.command
def chain(
    chains: InputFile,
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    chain_system: ChainIdSystem = "auth",
    scheduler_address: str | None = None,
    cache: CacheParameter | None = None,
    force: Annotated[bool, Parameter(negative="")] = False,
    _: Common | None = None,
) -> None:
    """Filter on chain.

    For each input PDB/mmCIF and chain combination write a PDB/mmCIF file with just the given chain
    and rename it to chain `A`. Filtering is done in parallel using a Dask cluster.

    Args:
        chains: CSV file with `pdb_id` and `chain` columns. Other columns are ignored.
        input_dir: Directory with PDB/mmCIF files.
            Expected filenames are `{pdb_id}.cif.gz`, `{pdb_id}.cif`, `{pdb_id}.pdb.gz` or `{pdb_id}.pdb`.
        output_dir: Directory to write the single-chain PDB/mmCIF files.
            Output files are in same format as input files.
        chain_system: System of chain ids in the input CSV.
            Set to 'label' to use chain ids assigned by PDB.
            See [docs](https://www.bonvinlab.org/protein_quest/autoapi/protein_quest/structure/chains.html#protein_quest.structure.chains.ChainIdSystem)
            for more information on chain id system.
        force: If not set and given chain exists and is only one,
            then file is copied unchanged. If set always rewrites and overwrites output files.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache or CacheParameter()

    with chains.open("r", encoding="utf-8") as chains_file:
        rows = list(csv.DictReader(chains_file))
    file2chain: set[tuple[Path, str]] = set()
    errors: list[FileNotFoundError] = []

    for row in rows:
        pdb_id = row["pdb_id"]
        chain = row["chain"]
        try:
            f = locate_structure_file(input_dir, pdb_id)
            file2chain.add((f, chain))
        except FileNotFoundError as e:
            errors.append(e)

    if errors:
        msg = f"Some structure files could not be found ({len(errors)} missing), skipping them"
        rprint(Panel(os.linesep.join(map(str, errors)), title=msg, style="red"))

    if not file2chain:
        msg = "No valid structure files found."
        raise ValueError(msg)

    results = filter_files_on_chain(
        file2chain,
        output_dir,
        chain_system=chain_system,
        force=force,
        scheduler_address=scheduler_address,
        copy_method=cache.copy_method,
    )

    nr_written = len([r for r in results if r.passed])

    rprint(f"Wrote {nr_written} single-chain PDB/mmCIF files to {output_dir}.")

    for result in results:
        if result.discard_reason:
            rprint(f"[red]Discarding {result.input_file} ({result.discard_reason})[/red]")


@filter_app.command
def residue(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    min_residues: MinResidues = 0,
    max_residues: MaxResidues = 10_000_000,
    write_stats: OutputFile | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter PDB/mmCIF files by number of residues in chain A.

    Filter PDB/mmCIF files by number of residues in chain A.

    Args:
        input_dir: Directory with PDB/mmCIF files (for example from 'filter chain').
        output_dir: Directory to write filtered PDB/mmCIF files. Files are copied without modification.
        min_residues: Min residues in chain A.
        max_residues: Max residues in chain A.
        write_stats: Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache or CacheParameter()

    if write_stats and str(write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    stats_lines = ["input_file,residue_count,passed,output_file"]

    nr_passed = 0
    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)
    rprint(f"Filtering {nr_total} files in {input_dir} directory by number of residues in chain A.")

    for r in filter_files_on_residues(
        input_files,
        output_dir,
        min_residues=min_residues,
        max_residues=max_residues,
        copy_method=cache.copy_method,
    ):
        stats_lines.append(f"{r.input_file},{r.residue_count},{r.passed},{r.output_file or ''}")
        if r.passed:
            nr_passed += 1

    if write_stats:
        write_lines(write_stats, stats_lines)
    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if write_stats and str(write_stats) != "-":
        rprint(f"Statistics written to {write_stats}")


@filter_app.command
def resolution(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    no_group_by_uniprot_accession: Annotated[bool, Parameter(negative="")] = False,
    top: PositiveInt = 1_000,
    no_coverage: Annotated[bool, Parameter(negative="")] = False,
    min_sequence_identity: NormFloat = 1.0,
    lax: Annotated[bool, Parameter(negative="")] = False,
    scheduler_address: str | None = None,
    write_stats: OutputFile | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter structure files by best resolution.

    * Structures with no resolution and no Uniprot accession are discarded.
    * Structures with lower resolution are preferred.
    * If resolution is the same, structures with more residues are preferred.
    * If resolution is missing, those structures are undesirable.
    * Structures with low sequence identity (smaller than 1) are undesirable.
    * Structures are grouped by Uniprot accession and then clustered by their coverage of the Uniprot sequence
      (can be disabled with `--no-group-by-uniprot-accession` and `--no-coverage`).

    To see how clustering was done use `protein-quest convert clusters` command.

    Args:
        input_dir: Directory structure files.
        output_dir: Directory to write the selected structure files.
        no_group_by_uniprot_accession: Disable grouping by Uniprot accession
            and use global top-N ranking across all files.
        top: Maximum number of files to keep.
        no_coverage: If not set, will take top by first grouping by uniprot accession
            and then clustering files by their coverage and then take the top.
            See
            [clustering documentation](https://www.bonvinlab.org/protein-quest/autoapi/protein_quest/clustering.html#protein_quest.pdbe.clustering.filter_pdbs_on_clustered_resolution)
            for details on the clustering and ordering criteria.
            If set will sort files per Uniprot accession by just their resolution.
        min_sequence_identity: Minimum sequence identity ratio to the Uniprot sequence for a structure to be passed.
            If not set then discards structures that are not fully identical to the Uniprot sequence.
            For example if set to 0.8 then structures that have sequence identity below 0.8 are discarded.
        lax: If set will passthrough files that do not have valid resolution, regardless of other flags like --top.
            By default filter is applied strictly and those files are discarded.
            Useful if you have a directory with files from different sources/methods and
            want to call multiple filters sequentially to discard files gradually.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        write_stats: Write filter statistics to file.
            In CSV format with columns:
            `<input_file>,<id>,<uniprot_accession>,<resolution>,<total_residue_count>,<is_alphafold>,<uniprot_start>,<uniprot_end>,<sequence_identity>,<chain_length>,<passed>,<output_file>,<discard_reason>,<discard_reason_type>`.
            Use `-` for stdout.
        cache: Cache options
        _: Common CLI options.

    """
    coverage = not no_coverage
    group_by = not no_group_by_uniprot_accession
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = cache or CacheParameter()

    if write_stats and str(write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)
    if not group_by:
        rprint(f"Filtering {nr_total} files in {input_dir} directory by global resolution ranking (no grouping).")
    else:
        rprint(f"Filtering {nr_total} files in {input_dir} directory by resolution grouped by uniprot accession.")

    results = list(
        filter_files_on_resolution(
            input_files,
            output_dir,
            top=top,
            coverage=coverage,
            group_by=group_by,
            min_sequence_identity=min_sequence_identity,
            lax=lax,
            copy_method=cache.copy_method,
            scheduler_address=scheduler_address,
        )
    )

    nr_passed = sum(1 for r in results if r.passed)
    if lax:
        nr_passed_due_to_lax = sum(1 for r in results if r.passed and r.discard_reason is not None)
        rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
        rprint(f"Additionally wrote {nr_passed_due_to_lax} files to {output_dir} directory due to lax mode.")
        if not write_stats:
            rprint(
                "[yellow]Note: You can use --write-stats to see which files were passed due to lax mode "
                "and their discard reasons.[/yellow]"
            )
    else:
        rprint(f"Wrote {nr_passed} files to {output_dir} directory.")

    if write_stats:
        write_resolution_stats(results, write_stats)
        if str(write_stats) != "-":
            rprint(f"Statistics written to {write_stats}")


@filter_app.command
def pdbe_quality(
    input_dir: InputDir,
    quality_json: InputFile,
    output_dir: OutputDir,
    /,
    *,
    minimal_geometry_quality: float = 50.0,
    top: PositiveInt | None = None,
    pass_given_resolution: Annotated[bool, Parameter(negative="")] = False,
    cluster_by_uniprot_accession_and_coverage: PositiveInt | None = None,
    write_stats: OutputFile | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
):
    """Filter PDB/mmCIF files by PDBe quality scores.

    Args:
        input_dir: Directory with PDB/mmCIF files.
        quality_json: JSON file with PDBe quality scores.
            Can be made with `protein-quest search pdbe-quality` command.
        output_dir: Directory to write filtered PDB/mmCIF files. Files are copied without modification.
        minimal_geometry_quality: Minimal geometry quality score to pass the filter.
        top: Maximum number of files to keep. If not given, top is same as number of files that pass the filter.
        pass_given_resolution: If set will passthrough files that have a valid resolution.
        cluster_by_uniprot_accession_and_coverage: Number of top structures to keep per UniProt cluster.
            Structures are grouped by UniProt accession, then clustered by residue-range overlap.
        write_stats: Write filter statistics to file.
            In CSV format with columns:
            `<pdb_id>,<input_file>,<geometry_quality>,<passed>,<output_file>,<reason>`.
            Use `-` for stdout.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.

    """
    cache = cache or CacheParameter()
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Reading PDBe quality scores from "%s"', quality_json)
    with quality_json.open("r", encoding="utf-8") as f:
        scores = converter.loads(f.read(), dict[str, Scores])

    logger.info('Finding structure files in "%s" directory that match ids from PDBe quality scores', input_dir)
    located_ids = locate_structure_files_by_id(set(scores.keys()), input_dir)

    logger.info("Filtering structure files by PDBe quality scores")
    if cluster_by_uniprot_accession_and_coverage is not None:
        results = filter_by_pdbe_quality_clustered(
            scores,
            located_ids,
            minimal_geometry_quality=minimal_geometry_quality,
            top=top,
            pass_given_resolution=pass_given_resolution,
            cluster_by_uniprot_accession_and_coverage=cluster_by_uniprot_accession_and_coverage,
        )
    else:
        results = filter_by_pdbe_quality(
            scores,
            located_ids,
            minimal_geometry_quality=minimal_geometry_quality,
            top=top,
            pass_given_resolution=pass_given_resolution,
        )

    for result in results:
        if result.passed and result.input_file is not None:
            copyfile(result.input_file, output_dir / result.input_file.name, cache.copy_method)

    rprint(f"Wrote {len([r for r in results if r.passed])} files to {output_dir} directory.")
    rprint(f"Discarded {len([r for r in results if not r.passed])} files due to any reason.")
    rprint(f"Discarded {len(located_ids.extras)} files in {input_dir} directory that were not in {quality_json}.")

    if write_stats:
        if str(write_stats) != "-":
            write_stats.parent.mkdir(parents=True, exist_ok=True)
        write_quality_stats_csv(results, write_stats, output_dir)
        rprint(f"Statistics written to {write_stats}")


@filter_app.command
def secondary_structure(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    filters: SecondaryStructureFilterQuery | None = None,
    write_stats: OutputFile | None = None,
    cache: CacheParameter | None = None,
    _: Common | None = None,
) -> None:
    """Filter PDB/mmCIF files by secondary structure.

    Filter PDB/mmCIF files by secondary structure.

    Args:
        input_dir: Directory with PDB/mmCIF files.
        output_dir: Directory to write filtered PDB/mmCIF files. Files are copied without modification.
        filters: Secondary structure filtering criteria.
        write_stats: Write filter statistics to file.
            In CSV format with columns:
            `<input_file>,<nr_residues>,<nr_helix_residues>,<nr_sheet_residues>,
            <helix_ratio>,<sheet_ratio>,<passed>,<output_file>`.
            Use `-` for stdout.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _: Common CLI options.
    """
    if filters is None:
        filters = SecondaryStructureFilterQuery()
    cache = cache or CacheParameter()
    output_dir.mkdir(parents=True, exist_ok=True)

    if write_stats and str(write_stats) != "-":
        write_stats.parent.mkdir(parents=True, exist_ok=True)

    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)

    stats_lines = [
        "input_file,nr_residues,nr_helix_residues,nr_sheet_residues,helix_ratio,sheet_ratio,passed,output_file"
    ]

    rprint(f"Filtering {nr_total} files in {input_dir} directory by secondary structure.")
    nr_passed = 0
    for input_file, result in filter_files_on_secondary_structure(input_files, query=filters):
        output_file: Path | None = None
        if result.passed:
            output_file = output_dir / input_file.name
            copyfile(input_file, output_file, cache.copy_method)
            nr_passed += 1
        stats_lines.append(
            f"{input_file},{result.stats.nr_residues},{result.stats.nr_helix_residues},"
            f"{result.stats.nr_sheet_residues},{round(result.stats.helix_ratio, 3)},"
            f"{round(result.stats.sheet_ratio, 3)},{result.passed},{output_file or ''}"
        )
        # TODO when some discard reason are resolvable by the user then
        # make cli have exit code of non-zero and raise a ExceptionGroup with all those errors
        # a user resolvable discard reason is for example multi-chain accession ambiguity
    if write_stats:
        write_lines(write_stats, stats_lines)
    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if str(write_stats) != "-":
        rprint(f"Statistics written to {write_stats}")
