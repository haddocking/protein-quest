"""Convert subcommands for protein-quest."""

from typing import Annotated

from cyclopts import App, Parameter
from tqdm import tqdm

from protein_quest.cli.common import CacheParameter, Common, InputDir, OutputDir, OutputFile, console, write_lines
from protein_quest.clustering_io import (
    cluster_results_by_accession,
    write_clusters_csv,
    write_condensed_distances_csv,
    write_dendrogram_nwk,
    write_linkage_matrix_csv,
    write_stats_csv,
)
from protein_quest.filters.resolution import load_resolution_statistics
from protein_quest.io import CifOutputFormat, convert_to_cif_files, glob_structure_files, read_structure
from protein_quest.structure import structure2uniprot_accessions

rprint = console.print


convert_app = App(name="convert", help="Convert files between formats")


@convert_app.command
def uniprot(
    input_dir: InputDir,
    output: OutputFile,
    /,
    *,
    grouped: Annotated[
        bool,
        Parameter(negative=""),
    ] = False,
    _common: Common | None = None,
) -> None:
    """Convert structure files to list of UniProt accessions.

    UniProt accessions are read from database reference of each structure.

    Args:
        input_dir: Directory with structure files. Supported extensions are .cif, .cif.gz, .pdb, .pdb.gz.
        output: Output text file with UniProt accessions (one per line). Use '-' for stdout.
        grouped: Whether to group accessions by structure file.
            If set output changes to `<structure_file1>,<acc1>\\n<structure_file1>,<acc2>` format.
        _common: Common CLI options.
    """
    input_files = sorted(glob_structure_files(input_dir))

    if grouped:
        lines = []
        for input_file in tqdm(input_files, unit="file"):
            s = read_structure(input_file)
            uniprot_accessions = structure2uniprot_accessions(s)
            lines.extend([f"{input_file},{uniprot_accession}" for uniprot_accession in sorted(uniprot_accessions)])
        write_lines(output, lines)
    else:
        uniprot_accessions: set[str] = set()
        for input_file in tqdm(input_files, unit="file"):
            s = read_structure(input_file)
            uniprot_accessions.update(structure2uniprot_accessions(s))
        write_lines(output, sorted(uniprot_accessions))


@convert_app.command
def structures(
    input_dir: InputDir,
    /,
    *,
    output_dir: OutputDir | None = None,
    output_format: CifOutputFormat = ".cif",
    cache: CacheParameter | None = None,
    _common: Common | None = None,
) -> None:
    """Convert structure files between formats.

    Convert structure files between formats.

    Args:
        input_dir: Directory with structure files.
            Supported extensions are .pdb, .pdb.gz, .ent, .ent.gz, .cif,
            .cif.gz, .bcif, .bcif.gz.
        output_dir: Directory to write converted structure files.
            If not given, files are written to input_dir.
        output_format: Output format for converted files. Supported values are .cif and .cif.gz.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _common: Common CLI options.
    """
    cache = cache or CacheParameter()
    output_dir = output_dir if output_dir is not None else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    input_files = sorted(glob_structure_files(input_dir))
    rprint(f"Converting {len(input_files)} files in {input_dir} directory to {output_format} format.")

    for _ in tqdm(
        convert_to_cif_files(
            input_files,
            output_dir,
            copy_method=cache.copy_method,
            output_format=output_format,
        ),
        total=len(input_files),
        unit="file",
    ):
        pass

    rprint(f"Converted {len(input_files)} files into {output_dir}.")


@convert_app.command
def clusters(
    input_dir: InputDir,
    output_dir: OutputDir,
    /,
    *,
    condensed_distances: Annotated[bool, Parameter(negative="")] = False,
    linkage_matrix: Annotated[bool, Parameter(negative="")] = False,
    dendrogram: Annotated[bool, Parameter(negative="")] = False,
    scheduler_address: str | None = None,
    _common: Common | None = None,
) -> None:
    """Cluster structures per UniProt accession and write clustering outputs.

    Always writes two CSV files into ``output_dir``:

    - ``clusters.csv``: one row per structure with cluster assignment.
    - ``stats.csv``: one row per UniProt accession with cluster summary.

    Can be used to investigate why `protein-quest filter resolution ...` or
    `protein-quest search pdbe --top_clustered_resolution_per_uniprot_accession ...` or
    `protein-quest search structure --top_clustered_resolution_per_uniprot_accession ...`
    keeps or discards certain structures
    by checking their intermediate cluster assignments and statistics.

    Args:
        input_dir: Directory with structure files.
        output_dir: Directory where clustering outputs are written.
        condensed_distances: Write per-accession condensed distance CSV files.
            Writes file per uniprot_accession with `<accession>_distances.csv` name pattern.
        linkage_matrix: Write per-accession linkage matrix CSV files.
            Writes file per uniprot_accession with `<accession>_linkage.csv` name pattern.
        dendrogram: Write per-accession dendrogram files in Newick format.
            Writes file per uniprot_accession with `<accession>_dendrogram.nwk` name pattern.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        _common: Common CLI options.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    input_files = sorted(glob_structure_files(input_dir))

    rprint(f"Clustering {len(input_files)} files in {input_dir} directory by UniProt accession.")
    stats = load_resolution_statistics(input_files, scheduler_address)
    results = cluster_results_by_accession(stats)

    write_clusters_csv(results, output_dir / "clusters.csv")
    write_stats_csv(results, output_dir / "stats.csv")

    if condensed_distances or linkage_matrix or dendrogram:
        for result in results:
            uniprot_accession = result.uniprot_accession

            structures = result.structures
            if len(structures) < 2:
                continue

            if condensed_distances:
                write_condensed_distances_csv(
                    result.condensed_distances,
                    structures,
                    output_dir / f"{uniprot_accession}_distances.csv",
                )

            if linkage_matrix or dendrogram:
                local_linkage_matrix = result.linkage_matrix
                if local_linkage_matrix is None:
                    continue

                if linkage_matrix:
                    write_linkage_matrix_csv(
                        local_linkage_matrix, structures, output_dir / f"{uniprot_accession}_linkage.csv"
                    )
                if dendrogram:
                    write_dendrogram_nwk(
                        local_linkage_matrix, structures, output_dir / f"{uniprot_accession}_dendrogram.nwk"
                    )

    rprint(f"Wrote clusters for {len(results)} uniprot accessions to {output_dir}")
