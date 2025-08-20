import argparse
import asyncio
import csv
import logging
import os
from collections.abc import Callable, Iterable
from importlib.util import find_spec
from io import TextIOWrapper
from pathlib import Path
from textwrap import dedent

from cattrs import structure
from rich import print as rprint
from rich.logging import RichHandler
from rich_argparse import ArgumentDefaultsRichHelpFormatter
from tqdm.rich import tqdm

from protein_quest.__version__ import __version__
from protein_quest.alphafold.confidence import ConfidenceFilterQuery, filter_files_on_confidence
from protein_quest.alphafold.fetch import DownloadableFormat, downloadable_formats
from protein_quest.alphafold.fetch import fetch_many as af_fetch
from protein_quest.filters import filter_files_on_chain, filter_files_on_residues
from protein_quest.go import Aspect, allowed_aspects, search_gene_ontology_term, write_go_terms_to_csv
from protein_quest.pdbe import fetch as pdbe_fetch
from protein_quest.pdbe.io import glob_structure_files
from protein_quest.uniprot import PdbResult, Query, search4af, search4pdb, search4uniprot

logger = logging.getLogger(__name__)


def _add_search_uniprot_parser(subparsers: argparse._SubParsersAction):
    """Add search uniprot subcommand parser."""
    parser = subparsers.add_parser(
        "uniprot",
        help="Search UniProt accessions",
        description="Search for UniProt accessions based on various criteria in the Uniprot SPARQL endpoint.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "output",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output text file for UniProt accessions (one per line). Use `-` for stdout.",
    )
    parser.add_argument("--taxon-id", type=str, help="NCBI Taxon ID, e.g. 9606 for Homo Sapiens")
    parser.add_argument(
        "--reviewed",
        action=argparse.BooleanOptionalAction,
        help="Reviewed=swissprot, no-reviewed=trembl. Default is uniprot=swissprot+trembl.",
        default=None,
    )
    parser.add_argument(
        "--subcellular-location-uniprot",
        type=str,
        help="Subcellular location label as used by UniProt (e.g. nucleus)",
    )
    parser.add_argument(
        "--subcellular-location-go",
        dest="subcellular_location_go",
        action="append",
        help="GO term(s) for subcellular location (e.g. GO:0005634). Can be given multiple times.",
    )
    parser.add_argument(
        "--molecular-function-go",
        dest="molecular_function_go",
        action="append",
        help="GO term(s) for molecular function (e.g. GO:0003677). Can be given multiple times.",
    )
    parser.add_argument("--limit", type=int, default=10_000, help="Maximum number of uniprot accessions to return")
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")


def _add_search_pdbe_parser(subparsers: argparse._SubParsersAction):
    """Add search pdbe subcommand parser."""
    parser = subparsers.add_parser(
        "pdbe",
        help="Search PDBe structures of given UniProt accessions",
        description="Search for PDB structures of given UniProt accessions in the Uniprot SPARQL endpoint.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "uniprot_accs",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line). Use `-` for stdin.",
    )
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help=dedent("""\
            Output CSV with `uniprot_acc`, `pdb_id`, `method`, `resolution`, `uniprot_chains`, `chain` columns.
            Where `uniprot_chains` is the raw UniProt chain string, for example `A=1-100`.
            and where `chain` is the first chain from `uniprot_chains`, for example `A`.
            Use `-` for stdout.
        """),
    )
    parser.add_argument(
        "--limit", type=int, default=10_000, help="Maximum number of PDB uniprot accessions combinations to return"
    )
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")


def _add_search_alphafold_parser(subparsers: argparse._SubParsersAction):
    """Add search alphafold subcommand parser."""
    parser = subparsers.add_parser(
        "alphafold",
        help="Search AlphaFold structures of given UniProt accessions",
        description="Search for AlphaFold structures of given UniProt accessions in the Uniprot SPARQL endpoint.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "uniprot_accs",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line). Use `-` for stdin.",
    )
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with AlphaFold IDs per UniProt accession. Use `-` for stdout.",
    )
    parser.add_argument(
        "--limit", type=int, default=10_000, help="Maximum number of Alphafold entry identifiers to return"
    )
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")


def _add_search_go_parser(subparsers: argparse._SubParsersAction):
    """Add search go subcommand parser"""
    parser = subparsers.add_parser(
        "go",
        help="Search for Gene Ontology (GO) terms",
        description="Search for Gene Ontology (GO) terms in the EBI QuickGO API.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "term",
        type=str,
        help="GO term to search for. For example `apoptosome`.",
    )
    parser.add_argument("--aspect", type=str, choices=allowed_aspects, help="Filter on aspect.")
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with GO term results. Use `-` for stdout.",
    )
    parser.add_argument("--limit", type=int, default=100, help="Maximum number of GO term results to return")


def _add_retrieve_pdbe_parser(subparsers: argparse._SubParsersAction):
    """Add retrieve pdbe subcommand parser."""
    parser = subparsers.add_parser(
        "pdbe",
        help="Retrieve PDBe gzipped mmCIF files for PDB IDs in CSV.",
        description=dedent("""\
            Retrieve mmCIF files from Protein Data Bank in Europe Knowledge Base (PDBe) website
            for unique PDB IDs listed in a CSV file.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "pdbe_csv",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="CSV file with `pdb_id` column. Other columns are ignored. Use `-` for stdin.",
    )
    parser.add_argument("output_dir", type=Path, help="Directory to store downloaded PDBe mmCIF files")
    parser.add_argument(
        "--max-parallel-downloads",
        type=int,
        default=5,
        help="Maximum number of parallel downloads",
    )


def _add_retrieve_alphafold_parser(subparsers: argparse._SubParsersAction):
    """Add retrieve alphafold subcommand parser."""
    parser = subparsers.add_parser(
        "alphafold",
        help="Retrieve AlphaFold files for IDs in CSV",
        description="Retrieve AlphaFold files from the AlphaFold Protein Structure Database.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "alphafold_csv",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="CSV file with `af_id` column. Other columns are ignored. Use `-` for stdin.",
    )
    parser.add_argument("output_dir", type=Path, help="Directory to store downloaded AlphaFold files")
    parser.add_argument(
        "--what-af-formats",
        type=str,
        action="append",
        choices=sorted(downloadable_formats),
        help=dedent("""AlphaFold formats to retrieve. Can be specified multiple times.
            Default is 'pdb'. Summary is always downloaded as `<entryId>.json`."""),
    )
    parser.add_argument(
        "--max-parallel-downloads",
        type=int,
        default=5,
        help="Maximum number of parallel downloads",
    )


def _add_filter_confidence_parser(subparsers: argparse._SubParsersAction):
    """Add filter confidence subcommand parser."""
    parser = subparsers.add_parser(
        "confidence",
        help="Filter AlphaFold mmcif/PDB files by confidence",
        description=dedent("""\
            Filter AlphaFold mmcif/PDB files by confidence (plDDT).
            Passed files are written with residues below threshold removed."""),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument("input_dir", type=Path, help="Directory with AlphaFold mmcif/PDB files")
    parser.add_argument("output_dir", type=Path, help="Directory to write filtered mmcif/PDB files")
    parser.add_argument("--confidence-threshold", type=float, default=70, help="pLDDT confidence threshold (0-100)")
    parser.add_argument(
        "--min-residues", type=int, default=0, help="Minimum number of high-confidence residues a structure should have"
    )
    parser.add_argument(
        "--max-residues",
        type=int,
        default=10_000_000,
        help="Maximum number of high-confidence residues a structure should have",
    )
    parser.add_argument(
        "--write-stats",
        type=argparse.FileType("w", encoding="UTF-8"),
        help=dedent("""\
            Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout."""),
    )


def _add_filter_chain_parser(subparsers: argparse._SubParsersAction):
    """Add filter chain subcommand parser."""
    parser = subparsers.add_parser(
        "chain",
        help="Filter on chain.",
        description=dedent("""\
            For each input PDB/mmCIF and chain combination
            write a PDB/mmCIF file with just the given chain
            and rename it to chain `A`.
            Filtering is done in parallel using a Dask cluster."""),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "chains",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="CSV file with `pdb_id` and `chain` columns. Other columns are ignored.",
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help=dedent("""\
        Directory with PDB/mmCIF files.
        Expected filenames are `{pdb_id}.cif.gz`, `{pdb_id}.cif`, `{pdb_id}.pdb.gz` or `{pdb_id}.pdb`.
    """),
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help=dedent("""\
        Directory to write the single-chain PDB/mmCIF files. Output files are in same format as input files."""),
    )
    parser.add_argument(
        "--scheduler-address",
        help="Address of the Dask scheduler to connect to. If not provided, will create a local cluster.",
    )


def _add_filter_residue_parser(subparsers: argparse._SubParsersAction):
    """Add filter residue subcommand parser."""
    parser = subparsers.add_parser(
        "residue",
        help="Filter PDB/mmCIF files by number of residues in chain A",
        description=dedent("""\
            Filter PDB/mmCIF files by number of residues in chain A.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument("input_dir", type=Path, help="Directory with PDB/mmCIF files (e.g., from 'filter chain')")
    parser.add_argument(
        "output_dir",
        type=Path,
        help=dedent("""\
        Directory to write filtered PDB/mmCIF files. Files are copied without modification.
    """),
    )
    parser.add_argument("--min-residues", type=int, default=0, help="Min residues in chain A")
    parser.add_argument("--max-residues", type=int, default=10_000_000, help="Max residues in chain A")
    parser.add_argument(
        "--write-stats",
        type=argparse.FileType("w", encoding="UTF-8"),
        help=dedent("""\
            Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout."""),
    )


def _add_search_subcommands(subparsers: argparse._SubParsersAction):
    """Add search command and its subcommands."""
    parser = subparsers.add_parser(
        "search",
        help="Search data sources",
        description="Search various things online.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    subsubparsers = parser.add_subparsers(dest="search_cmd", required=True)

    _add_search_uniprot_parser(subsubparsers)
    _add_search_pdbe_parser(subsubparsers)
    _add_search_alphafold_parser(subsubparsers)
    _add_search_go_parser(subsubparsers)


def _add_retrieve_subcommands(subparsers: argparse._SubParsersAction):
    """Add retrieve command and its subcommands."""
    parser = subparsers.add_parser(
        "retrieve",
        help="Retrieve structure files",
        description="Retrieve structure files from online resources.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    subsubparsers = parser.add_subparsers(dest="retrieve_cmd", required=True)

    _add_retrieve_pdbe_parser(subsubparsers)
    _add_retrieve_alphafold_parser(subsubparsers)


def _add_filter_subcommands(subparsers: argparse._SubParsersAction):
    """Add filter command and its subcommands."""
    parser = subparsers.add_parser("filter", help="Filter files", formatter_class=ArgumentDefaultsRichHelpFormatter)
    subsubparsers = parser.add_subparsers(dest="filter_cmd", required=True)

    _add_filter_confidence_parser(subsubparsers)
    _add_filter_chain_parser(subsubparsers)
    _add_filter_residue_parser(subsubparsers)


def _add_mcp_command(subparsers: argparse._SubParsersAction):
    """Add MCP command."""

    parser = subparsers.add_parser(
        "mcp",
        help="Run Model Context Protocol (MCP) server",
        description=(
            "Run Model Context Protocol (MCP) server. "
            "Can be used by agentic LLMs like Claude Sonnet 4 as a set of tools."
        ),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "--transport", default="stdio", choices=["stdio", "http", "streamable-http"], help="Transport protocol to use"
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind the server to")
    parser.add_argument("--port", default=8000, type=int, help="Port to bind the server to")


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Protein Quest CLI", prog="protein-quest", formatter_class=ArgumentDefaultsRichHelpFormatter
    )
    parser.add_argument("--log-level", default="WARNING", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_search_subcommands(subparsers)
    _add_retrieve_subcommands(subparsers)
    _add_filter_subcommands(subparsers)
    _add_mcp_command(subparsers)

    return parser


def main():
    """Main entry point for the CLI."""
    parser = make_parser()
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level, handlers=[RichHandler(show_level=False)])

    # Dispatch table to reduce complexity
    cmd = args.command
    sub = getattr(args, f"{cmd}_cmd", None)
    handler = HANDLERS.get((cmd, sub))
    if handler is None:
        msg = f"Unknown command: {cmd} {sub}"
        raise SystemExit(msg)
    handler(args)


def _handle_search_uniprot(args):
    taxon_id = args.taxon_id
    reviewed = args.reviewed
    subcellular_location_uniprot = args.subcellular_location_uniprot
    subcellular_location_go = args.subcellular_location_go
    molecular_function_go = args.molecular_function_go
    limit = args.limit
    timeout = args.timeout
    output_file = args.output

    query = structure(
        {
            "taxon_id": taxon_id,
            "reviewed": reviewed,
            "subcellular_location_uniprot": subcellular_location_uniprot,
            "subcellular_location_go": _as_scalar_or_list(subcellular_location_go),
            "molecular_function_go": _as_scalar_or_list(molecular_function_go),
        },
        Query,
    )
    rprint("Searching for UniProt accessions")
    accs = search4uniprot(query=query, limit=limit, timeout=timeout)
    rprint(f"Found {len(accs)} UniProt accessions, written to {output_file.name}")
    _write_lines(output_file, sorted(accs))


def _handle_search_pdbe(args):
    uniprot_accs = args.uniprot_accs
    limit = args.limit
    timeout = args.timeout
    output_csv = args.output_csv

    accs = set(_read_lines(uniprot_accs))
    rprint(f"Finding PDB entries for {len(accs)} uniprot accessions")
    results = search4pdb(accs, limit=limit, timeout=timeout)
    total_pdbs = sum([len(v) for v in results.values()])
    rprint(f"Found {total_pdbs} PDB entries for {len(results)} uniprot accessions")
    rprint(f"Written to {output_csv.name}")
    _write_pdbe_csv(output_csv, results)


def _handle_search_alphafold(args):
    uniprot_accs = args.uniprot_accs
    limit = args.limit
    timeout = args.timeout
    output_csv = args.output_csv

    accs = _read_lines(uniprot_accs)
    rprint(f"Finding AlphaFold entries for {len(accs)} uniprot accessions")
    results = search4af(accs, limit=limit, timeout=timeout)
    rprint(f"Found {len(results)} AlphaFold entries, written to {output_csv.name}")
    _write_alphafold_csv(output_csv, results)


def _handle_search_go(args):
    term = structure(args.term, str)
    aspect: Aspect | None = args.aspect
    limit = structure(args.limit, int)
    output_csv: TextIOWrapper = args.output_csv

    if aspect:
        rprint(f"Searching for GO terms matching '{term}' with aspect '{aspect}'")
    else:
        rprint(f"Searching for GO terms matching '{term}'")
    results = asyncio.run(search_gene_ontology_term(term, aspect=aspect, limit=limit))
    rprint(f"Found {len(results)} GO terms, written to {output_csv.name}")
    write_go_terms_to_csv(results, output_csv)


def _handle_retrieve_pdbe(args):
    pdbe_csv = args.pdbe_csv
    output_dir = args.output_dir
    max_parallel_downloads = args.max_parallel_downloads

    pdb_ids = _read_pdbe_ids_from_csv(pdbe_csv)
    rprint(f"Retrieving {len(pdb_ids)} PDBe entries")
    result = pdbe_fetch.fetch(pdb_ids, output_dir, max_parallel_downloads=max_parallel_downloads)
    rprint(f"Retrieved {len(result)} PDBe entries")


def _handle_retrieve_alphafold(args):
    download_dir = args.output_dir
    what_af_formats = args.what_af_formats
    alphafold_csv = args.alphafold_csv
    max_parallel_downloads = args.max_parallel_downloads

    if what_af_formats is None:
        what_af_formats = {"pdb"}

    # TODO besides `uniprot_acc,af_id\n` csv also allow headless single column format
    #
    af_ids = [r["af_id"] for r in _read_alphafold_csv(alphafold_csv)]
    validated_what: set[DownloadableFormat] = structure(what_af_formats, set[DownloadableFormat])
    rprint(f"Retrieving {len(af_ids)} AlphaFold entries with formats {validated_what}")
    afs = af_fetch(af_ids, download_dir, what=validated_what, max_parallel_downloads=max_parallel_downloads)
    total_nr_files = sum(af.nr_of_files() for af in afs)
    rprint(f"Retrieved {total_nr_files} AlphaFold files and {len(afs)} summaries, written to {download_dir}")


def _handle_filter_confidence(args: argparse.Namespace):
    # we are repeating types here and in add_argument call
    # TODO replace argparse with modern alternative like cyclopts
    # to get rid of duplication
    input_dir = structure(args.input_dir, Path)
    output_dir = structure(args.output_dir, Path)
    confidence_threshold = structure(args.confidence_threshold, float)
    # TODO add min/max
    min_residues = structure(args.min_residues, int)
    max_residues = structure(args.max_residues, int)
    stats_file: TextIOWrapper | None = args.write_stats

    output_dir.mkdir(parents=True, exist_ok=True)
    input_files = sorted(glob_structure_files(input_dir))
    nr_input_files = len(input_files)
    rprint(f"Starting confidence filtering of {nr_input_files} mmcif/PDB files in {input_dir} directory.")
    query = structure(
        {
            "confidence": confidence_threshold,
            "min_threshold": min_residues,
            "max_threshold": max_residues,
        },
        ConfidenceFilterQuery,
    )
    if stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["input_file", "residue_count", "passed", "output_file"])

    passed_count = 0
    for r in tqdm(filter_files_on_confidence(input_files, query, output_dir), total=len(input_files), unit="file"):
        if r.filtered_file:
            passed_count += 1
        if stats_file:
            writer.writerow([r.input_file, r.count, r.filtered_file is not None, r.filtered_file])

    rprint(f"Filtered {passed_count} mmcif/PDB files by confidence, written to {output_dir} directory")
    if stats_file:
        rprint(f"Statistics written to {stats_file.name}")


def _handle_filter_chain(args):
    input_dir = args.input_dir
    output_dir = args.output_dir
    pdb_id2chain_mapping_file = args.chains
    scheduler_address = args.scheduler_address

    rows = list(_iter_pdbe_csv(pdb_id2chain_mapping_file))
    id2chains: dict[str, str] = {row["pdb_id"]: row["chain"] for row in rows}

    new_files = filter_files_on_chain(input_dir, id2chains, output_dir, scheduler_address)

    nr_written = len([r for r in new_files if r[2] is not None])

    rprint(f"Wrote {nr_written} single-chain PDB/mmCIF files to {output_dir}.")


def _handle_filter_residue(args):
    input_dir = structure(args.input_dir, Path)
    output_dir = structure(args.output_dir, Path)
    min_residues = structure(args.min_residues, int)
    max_residues = structure(args.max_residues, int)
    stats_file: TextIOWrapper | None = args.write_stats

    if stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["input_file", "residue_count", "passed", "output_file"])

    nr_passed = 0
    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)
    rprint(f"Filtering {nr_total} files in {input_dir} directory by number of residues in chain A.")
    for r in filter_files_on_residues(input_files, output_dir, min_residues=min_residues, max_residues=max_residues):
        if stats_file:
            writer.writerow([r.input_file, r.residue_count, r.passed, r.output_file])
        if r.passed:
            nr_passed += 1

    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if stats_file:
        rprint(f"Statistics written to {stats_file.name}")


def _handle_mcp(args):
    if find_spec("fastmcp") is None:
        msg = "Unable to start MCP server, please install `protein-quest[mcp]`."
        raise ImportError(msg)

    from protein_quest.mcp_server import mcp  # noqa: PLC0415

    if args.transport == "stdio":
        mcp.run(transport=args.transport)
    else:
        mcp.run(transport=args.transport, host=args.host, port=args.port)


HANDLERS: dict[tuple[str, str | None], Callable] = {
    ("search", "uniprot"): _handle_search_uniprot,
    ("search", "pdbe"): _handle_search_pdbe,
    ("search", "alphafold"): _handle_search_alphafold,
    ("search", "go"): _handle_search_go,
    ("retrieve", "pdbe"): _handle_retrieve_pdbe,
    ("retrieve", "alphafold"): _handle_retrieve_alphafold,
    ("filter", "confidence"): _handle_filter_confidence,
    ("filter", "chain"): _handle_filter_chain,
    ("filter", "residue"): _handle_filter_residue,
    ("mcp", None): _handle_mcp,
}


def _as_scalar_or_list(values: list[str] | None) -> str | list[str] | None:
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    return values


def _read_lines(file: TextIOWrapper) -> list[str]:
    return [line.strip() for line in file]


def _make_sure_parent_exists(file: TextIOWrapper):
    if file.name != "<stdout>":
        Path(file.name).parent.mkdir(parents=True, exist_ok=True)


def _write_lines(file: TextIOWrapper, lines: Iterable[str]):
    _make_sure_parent_exists(file)
    file.writelines(line + os.linesep for line in lines)


def _write_pdbe_csv(path: TextIOWrapper, data: dict[str, set[PdbResult]]):
    _make_sure_parent_exists(path)
    fieldnames = ["uniprot_acc", "pdb_id", "method", "resolution", "uniprot_chains", "chain"]
    writer = csv.DictWriter(path, fieldnames=fieldnames)
    writer.writeheader()
    for uniprot_acc, entries in sorted(data.items()):
        for e in sorted(entries, key=lambda x: (x.id, x.method)):
            writer.writerow(
                {
                    "uniprot_acc": uniprot_acc,
                    "pdb_id": e.id,
                    "method": e.method,
                    "resolution": e.resolution or "",
                    "uniprot_chains": e.uniprot_chains,
                    "chain": e.chain,
                }
            )


def _write_alphafold_csv(file: TextIOWrapper, data: dict[str, set[str]]):
    _make_sure_parent_exists(file)
    fieldnames = ["uniprot_acc", "af_id"]

    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    for uniprot_acc, af_ids in sorted(data.items()):
        for af_id in sorted(af_ids):
            writer.writerow({"uniprot_acc": uniprot_acc, "af_id": af_id})


def _read_alphafold_csv(file: TextIOWrapper):
    reader = csv.DictReader(file)
    yield from reader


def _iter_pdbe_csv(file: TextIOWrapper):
    reader = csv.DictReader(file)
    # Expect columns: uniprot_acc, pdb_id, method, resolution, uniprot_chains
    yield from reader


def _read_pdbe_ids_from_csv(file: TextIOWrapper) -> set[str]:
    return {row["pdb_id"] for row in _iter_pdbe_csv(file)}
