import argparse
import csv
import json
import logging
import os
from collections.abc import Callable, Iterable
from dataclasses import asdict
from io import TextIOWrapper
from pathlib import Path
from shutil import copyfile
from textwrap import dedent

from cattrs import structure
from rich import print as rprint
from rich.logging import RichHandler
from rich_argparse import ArgumentDefaultsRichHelpFormatter
from tqdm.rich import tqdm

from protein_quest.__version__ import __version__
from protein_quest.alphafold import DownloadableFormat, downloadable_formats
from protein_quest.alphafold import fetch_many as af_fetch
from protein_quest.alphafold.confidence import DensityFilterQuery, filter_on_density
from protein_quest.alphafold.entry_summary import EntrySummary
from protein_quest.pdbe import fetch as pdbe_fetch
from protein_quest.pdbe.io import (
    is_chain_in_residues_range,
    write_single_chain_pdb_file,
)
from protein_quest.uniprot import Query, search4af, search4pdb, search4uniprot

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
        help="Output CSV with PDB info per UniProt accession. Use `-` for stdout.",
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


def _add_retrieve_pdbe_parser(subparsers: argparse._SubParsersAction):
    """Add retrieve pdbe subcommand parser."""
    parser = subparsers.add_parser(
        "pdbe",
        help="Retrieve PDBe gzipped mmCIF files for PDB IDs in CSV.",
        description=dedent("""\
            Retrieve mmCIF files from Protein Data Bank in Europe Knowledge Base (PDBe) website
            for PDB IDs listed in a CSV file.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "pdbe_csv",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="CSV file with `pdb_id` column. Use `-` for stdin.",
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
        help="CSV file with `af_id` column. Use `-` for stdin.",
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
        help="Filter AlphaFold PDBs by confidence",
        description=dedent("""\
            Filter AlphaFold PDB files by confidence.
            Passed files are written with residues below threshold removed."""),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument("input_dir", type=Path, help="Directory with AlphaFold PDB files")
    parser.add_argument("output_dir", type=Path, help="Directory to write filtered PDB files")
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


def _add_filter_chain_parser(subparsers: argparse._SubParsersAction):
    """Add filter chain subcommand parser."""
    parser = subparsers.add_parser(
        "chain",
        help="Keep first UniProt chain from PDBe mmCIF",
        description=dedent("""\
            For each input PDB and uniprot combination
            write a PDB file with just the first chain belonging to that Uniprot accession
            and rename it to chain `A`."""),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "pdbe_csv",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="CSV file with `pdb_id`, `uniprot_id` and `uniprot_chain` columns.",
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help=dedent("""\
        Directory with PDB/mmCIF files.
        Allowed formats are `.cif.gz`, `.cif`, `.pdb.gz` or `.pdb`.
    """),
    )
    parser.add_argument("output_dir", type=Path, help="Directory to write single-chain PDB files")


def _add_filter_residue_parser(subparsers: argparse._SubParsersAction):
    """Add filter residue subcommand parser."""
    parser = subparsers.add_parser(
        "residue",
        help="Filter PDBs by number of residues in chain A",
        description=dedent("""\
            Filter PDBs by number of residues in chain A.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    # TODO besides *.pdb also allow *.cif, *.pdb.gz and *.cif.gz
    parser.add_argument("input_dir", type=Path, help="Directory with PDB files (e.g., from 'filter chain')")
    parser.add_argument("output_dir", type=Path, help="Directory to write filtered PDB files")
    parser.add_argument("--min-residues", type=int, default=0, help="Min residues in chain A")
    parser.add_argument("--max-residues", type=int, default=10_000_000, help="Max residues in chain A")


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
    query = structure(
        {
            "taxon_id": args.taxon_id,
            "reviewed": args.reviewed,
            "subcellular_location_uniprot": args.subcellular_location_uniprot,
            "subcellular_location_go": _as_scalar_or_list(args.subcellular_location_go),
            "molecular_function_go": _as_scalar_or_list(args.molecular_function_go),
        },
        Query,
    )
    rprint("Searching for UniProt accessions")
    accs = search4uniprot(query=query, limit=args.limit, timeout=args.timeout)
    rprint(f"Found {len(accs)} UniProt accessions, written to {args.output.name}")
    _write_lines(args.output, sorted(accs))


def _handle_search_pdbe(args):
    accs = set(_read_lines(args.uniprot_accs))
    rprint(f"Finding PDB entries for {len(accs)} uniprot accessions")
    results = search4pdb(accs, limit=args.limit, timeout=args.timeout)
    total_pdbs = sum([len(v) for v in results.values()])
    rprint(f"Found {total_pdbs} PDB entries for {len(results)} uniprot accessions")
    rprint(f"Written to {args.output_csv.name}")
    _write_pdbe_csv(args.output_csv, results)


def _handle_search_alphafold(args):
    accs = _read_lines(args.uniprot_accs)
    rprint(f"Finding AlphaFold entries for {len(accs)} uniprot accessions")
    results = search4af(accs, limit=args.limit, timeout=args.timeout)
    rprint(f"Found {len(results)} AlphaFold entries, written to {args.output_csv.name}")
    _write_alphafold_csv(args.output_csv, results)


def _write_alphafold_summary(summary: EntrySummary, download_dir: Path):
    data = asdict(summary)
    fn = download_dir / f"{summary.entryId}.json"
    fn.write_text(json.dumps(data, indent=2))


def _handle_retrieve_pdbe(args):
    pdb_ids = _read_pdbe_ids_from_csv(args.pdbe_csv)
    rprint(f"Retrieving {len(pdb_ids)} PDBe entries")
    result = pdbe_fetch.fetch(pdb_ids, args.output_dir, max_parallel_downloads=args.max_parallel_downloads)
    rprint(f"Retrieved {len(result)} PDBe entries")


def _handle_retrieve_alphafold(args):
    download_dir = args.output_dir
    what_af_formats = args.what_af_formats
    if what_af_formats is None:
        what_af_formats = {"pdb"}

    # TODO besides `uniprot_acc,af_id\n` csv also allow headless single column format
    #
    af_ids = [r["af_id"] for r in _read_alphafold_csv(args.alphafold_csv)]
    validated_what = structure(what_af_formats, set[DownloadableFormat])
    rprint(f"Retrieving {len(af_ids)} AlphaFold entries with formats {validated_what}")
    afs = af_fetch(af_ids, download_dir, what=validated_what, max_parallel_downloads=args.max_parallel_downloads)

    for af in tqdm(afs, unit="entry", desc="Writing summaries to disk"):
        # TODO move writing of summary.json to af_fetch function and do concurrently
        if af.summary is None:
            continue
        _write_alphafold_summary(af.summary, download_dir)

    total_nr_files = sum(af.nr_of_files() for af in afs)
    rprint(f"Retrieved {total_nr_files} AlphaFold files and {len(afs)} summaries, written to {download_dir}")


def _handle_filter_confidence(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdb_files = sorted(args.input_dir.glob("*.pdb"))
    rprint(f"Filtering {len(pdb_files)} AlphaFold PDB files from {args.input_dir} directory by confidence")
    query = structure(
        {
            "confidence": args.confidence_threshold,
            "min_threshold": args.min_residues,
            "max_threshold": args.max_residues,
        },
        DensityFilterQuery,
    )
    passed_count = 0
    for r in tqdm(list(filter_on_density(pdb_files, query, args.output_dir)), unit="pdb"):
        # TODO log the nr of residues in a csv file if --store-count is given
        if r.density_filtered_file:
            passed_count += 1

    rprint(f"Filtered {passed_count} PDB files by confidence, written to {args.output_dir} directory")


def _handle_filter_chain(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # TODO allow simpler csv, uniprot_chain `A=1-100` can simplified to chain 'A'
    # TODO use range from uniprot_chain to only write residues in that range
    # TODO make handler smaller, by moving code to functions that do not use args.*
    rows = list(_iter_pdbe_csv(args.pdbe_csv))
    rprint(f"Filtering {len(rows)} PDBe entries for uniprot chain")
    nr_written = 0
    for row in tqdm(rows, unit="file"):
        pdb_id = row["pdb_id"]
        input_file = _locate_structure_file(args.input_dir, pdb_id)
        # TODO allow to specify output format, similar to input
        output_file = write_single_chain_pdb_file(
            input_file, row["uniprot_chains"], row["uniprot_acc"], args.output_dir
        )
        if output_file:
            nr_written += 1

    rprint(f"Wrote {nr_written} single-chain PDB files to {args.output_dir}.")


def _locate_structure_file(root: Path, pdb_id: str) -> Path:
    exts = [".cif.gz", ".cif", ".pdb.gz", ".pdb"]
    for ext in exts:
        candidate = root / f"{pdb_id.lower()}{ext}"
        if candidate.exists():
            return candidate
    msg = f"No structure file found for {pdb_id} in {root}"
    raise FileNotFoundError(msg)


def _handle_filter_residue(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    passed_count = 0
    # TODO run in parallel using dask.distributed
    # TODO make handler smaller, by moving code to functions that do not use args.*
    input_files = sorted(args.input_dir.glob("*.pdb"))
    for pdb in tqdm(input_files, unit="file"):
        # TODO log the nr of residues in a csv file if --store-count is given
        if not is_chain_in_residues_range(pdb, args.min_residues, args.max_residues, chain="A"):
            continue
        copyfile(pdb, args.output_dir / pdb.name)
        passed_count += 1

    discarded_count = len(input_files) - passed_count
    rprint(f"Completed filtering on residues. Filtered {passed_count} and discarded {discarded_count} pdb files.")


HANDLERS: dict[tuple[str, str | None], Callable] = {
    ("search", "uniprot"): _handle_search_uniprot,
    ("search", "pdbe"): _handle_search_pdbe,
    ("search", "alphafold"): _handle_search_alphafold,
    ("retrieve", "pdbe"): _handle_retrieve_pdbe,
    ("retrieve", "alphafold"): _handle_retrieve_alphafold,
    ("filter", "confidence"): _handle_filter_confidence,
    ("filter", "chain"): _handle_filter_chain,
    ("filter", "residue"): _handle_filter_residue,
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


def _write_pdbe_csv(path: TextIOWrapper, data: dict[str, set]):
    _make_sure_parent_exists(path)
    fieldnames = ["uniprot_acc", "pdb_id", "method", "resolution", "uniprot_chains"]
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
