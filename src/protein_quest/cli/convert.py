"""Convert subcommands for protein-quest."""

import csv
from pathlib import Path
from typing import Annotated

from cyclopts import App, Parameter
from tqdm import tqdm

from protein_quest.cli.common import (
    CacheParameter,
    Common,
    InputDir,
    InputFile,
    OutputDir,
    OutputFile,
    console,
    write_lines,
)
from protein_quest.clustering_io import (
    cluster_results_by_accession,
    write_clusters_csv,
    write_condensed_distances_by_accession_csv,
    write_dendrogram_nwk,
    write_linkage_matrix_by_accession_csv,
    write_stats_csv,
)
from protein_quest.filters.resolution import load_resolution_statistics
from protein_quest.structure.chains import ChainIdSystem
from protein_quest.structure.convert import convert_to_cif_files
from protein_quest.structure.files import glob_structure_files
from protein_quest.structure.formats import read_structure
from protein_quest.structure.types import CifOutputFormat
from protein_quest.structure.uniprot import structure2uniprot_accessions
from protein_quest.uniprot_chains import Pdb2UniprotChainsMapping, UniprotChainMapping, parse_uniprot_chains

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


def _read_pdb2uniprot_csv(uniprots: Path | None) -> Pdb2UniprotChainsMapping:
    """Read CSV file with PDB id to UniProt chain/range mappings.

    Expects 3 columns: `pdb_id,uniprot_accession,uniprot_chains`.

    Args:
        uniprots: CSV file with PDB to UniProt mappings. If None, returns empty dictionary.

    Returns:
        Dictionary mapping PDB ID to UniProt chain mappings.
    """
    uniprot_ref_dict: Pdb2UniprotChainsMapping = {}
    if uniprots is None:
        return uniprot_ref_dict

    with uniprots.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = row["pdb_id"]
            uniprot_accession = row["uniprot_accession"]
            uniprot_chains = row["uniprot_chains"]
            if pdb_id and uniprot_accession and uniprot_chains:
                if pdb_id not in uniprot_ref_dict:
                    uniprot_ref_dict[pdb_id] = set()
                uniprot_ref_dict[pdb_id].add(
                    UniprotChainMapping(
                        uniprot_accession=uniprot_accession,
                        chain_ranges=parse_uniprot_chains(uniprot_chains),
                    )
                )
    return uniprot_ref_dict


@convert_app.command
def structures(
    input_dir: InputDir,
    /,
    *,
    output_dir: OutputDir | None = None,
    output_format: CifOutputFormat = ".cif.gz",
    uniprots: InputFile | None = None,
    chain_system: ChainIdSystem = "auth",
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
        uniprots: Supply Uniprot to PDB id and chain mappings.
            Adds UniProt accessions to structures that are missing them based on the provided mapping.
            The supplied file must be in CSV format with 3 columns: `pdb_id,uniprot_accession,uniprot_chains`.
            (column order does not matter)
            This CSV file can be generated with `protein-quest search pdbe ...`.
        chain_system: System of chain ids in the input CSV.
            Set to 'label' to use chain ids assigned by PDB.
            See [docs](https://www.bonvinlab.org/protein_quest/autoapi/protein_quest/structure/chains.html#protein_quest.structure.chains.ChainIdSystem)
            for more information on chain id system.
        output_format: Output format for converted files. Supported values are .cif and .cif.gz.
        cache: Cache options including no_cache, cache_dir, and copy_method.
        _common: Common CLI options.
    """
    cache = cache or CacheParameter()
    output_dir = output_dir if output_dir is not None else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    input_files = sorted(glob_structure_files(input_dir))
    rprint(f"Converting {len(input_files)} files in {input_dir} directory to {output_format} format.")

    pdb2uniprot = _read_pdb2uniprot_csv(uniprots)

    for _ in tqdm(
        convert_to_cif_files(
            input_files,
            output_dir,
            copy_method=cache.copy_method,
            output_format=output_format,
            pdb2uniprot=pdb2uniprot,
            chain_system=chain_system,
        ),
        total=len(input_files),
        unit="file",
    ):
        pass

    rprint(f"Converted {len(input_files)} files into {output_dir}.")


@convert_app.command
def clusters(
    input_dir: InputDir,
    output_file: OutputFile,
    /,
    *,
    stats: OutputFile | None = None,
    condensed_distances: OutputFile | None = None,
    linkage_matrix: OutputFile | None = None,
    dendrogram: OutputDir | None = None,
    scheduler_address: str | None = None,
    _common: Common | None = None,
) -> None:
    """Group structures per UniProt accession and cluster based on coverage.

    Always writes one CSV file:

    - ``output_file``: one row per structure with cluster assignment.


    Can be used to investigate why `protein-quest filter resolution ...` or
    `protein-quest search pdbe --top_clustered_resolution_per_uniprot_accession ...` or
    `protein-quest search structure --top_clustered_resolution_per_uniprot_accession ...`
    keeps or discards certain structures
    by checking their intermediate cluster assignments and statistics.

    See
    [clustering documentation](https://www.bonvinlab.org/protein-quest/autoapi/protein_quest/clustering.html#protein_quest.pdbe.clustering.filter_pdbs_on_clustered_resolution)
    for details on the clustering and ordering criteria.

    Args:
        input_dir: Directory with structure files.
        output_file: Output CSV file with cluster assignments.
        stats: Optional output CSV file with per-accession cluster summary.
            Only written when provided.
        condensed_distances: Optional output CSV file with condensed distances from all accessions.
            Writes a single file with ``uniprot_accession`` as first column.
        linkage_matrix: Optional output CSV file with linkage rows from all accessions.
            Writes a single file with ``uniprot_accession`` as first column.
        dendrogram: Optional output directory for per-accession Newick files.
            Writes ``<accession>_dendrogram.nwk`` files into this directory.
        scheduler_address: Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially.
        _common: Common CLI options.
    """
    input_files = sorted(glob_structure_files(input_dir))

    rprint(f"Clustering {len(input_files)} files in {input_dir} directory by UniProt accession.")
    resolution_stats = load_resolution_statistics(input_files, scheduler_address)
    results = cluster_results_by_accession(resolution_stats)

    write_clusters_csv(results, output_file)
    if stats is not None:
        write_stats_csv(results, stats)
    if condensed_distances is not None:
        write_condensed_distances_by_accession_csv(results, condensed_distances)
    if linkage_matrix is not None:
        write_linkage_matrix_by_accession_csv(results, linkage_matrix)

    if dendrogram is not None:
        dendrogram.mkdir(parents=True, exist_ok=True)
        for result in results:
            uniprot_accession = result.uniprot_accession

            structures = result.structures
            if len(structures) < 2:
                continue

            local_linkage_matrix = result.linkage_matrix
            if local_linkage_matrix is None:
                continue

            write_dendrogram_nwk(local_linkage_matrix, structures, dendrogram / f"{uniprot_accession}_dendrogram.nwk")

    rprint(f"Wrote clusters for {len(results)} uniprot accessions to {output_file}")
