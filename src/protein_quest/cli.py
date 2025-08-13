import argparse
import csv
import json
import logging
from collections.abc import Callable, Iterable
from dataclasses import asdict
from pathlib import Path
from textwrap import dedent

from rich import print as rprint
from rich.logging import RichHandler
from tqdm.rich import tqdm

from protein_quest.__version__ import __version__
from protein_quest.alphafold import downloadable_formats
from protein_quest.alphafold import fetch_many as af_fetch
from protein_quest.alphafold.confidence import DensityFilterQuery, filter_on_density
from protein_quest.alphafold.entry_summary import EntrySummary
from protein_quest.pdbe import fetch as pdbe_fetch
from protein_quest.pdbe.io import (
    ProteinPdbRow,
    SingleChainQuery,
    nr_residues_in_chain,
    write_single_chain_pdb_file,
)
from protein_quest.uniprot import Query, search4af, search4pdb, search4uniprot

logger = logging.getLogger(__name__)


def _add_search_uniprot_parser(search_subparsers: argparse._SubParsersAction):
    """Add search uniprot subcommand parser."""
    sp_uniprot = search_subparsers.add_parser("uniprot", help="Search UniProt accessions")
    sp_uniprot.add_argument("output", type=Path, help="Output text file for UniProt accessions (one per line)")
    sp_uniprot.add_argument("--taxon-id", type=str, help="NCBI Taxon ID, e.g. 9606")
    sp_uniprot.add_argument("--reviewed", action="store_true", help="Only reviewed entries")
    sp_uniprot.add_argument(
        "--unreviewed", action="store_true", help="Only unreviewed entries (ignored if --reviewed is set)"
    )
    sp_uniprot.add_argument(
        "--subcellular-location-uniprot",
        type=str,
        help="Subcellular location label as used by UniProt (e.g. nucleus)",
    )
    sp_uniprot.add_argument(
        "--subcellular-location-go",
        dest="subcellular_location_go",
        action="append",
        help="GO term(s) for subcellular location (e.g. GO:0005634). Can be given multiple times.",
    )
    sp_uniprot.add_argument(
        "--molecular-function-go",
        dest="molecular_function_go",
        action="append",
        help="GO term(s) for molecular function (e.g. GO:0003677). Can be given multiple times.",
    )
    sp_uniprot.add_argument("--limit", type=int, default=10_000, help="SPARQL result limit")
    sp_uniprot.add_argument("--timeout", type=int, default=1_800, help="SPARQL timeout (seconds)")


def _add_search_pdbe_parser(search_subparsers: argparse._SubParsersAction):
    """Add search pdbe subcommand parser."""
    sp_pdbe = search_subparsers.add_parser("pdbe", help="Search PDBe structures for UniProt accessions")
    sp_pdbe.add_argument("uniprot_accs", type=Path, help="Text file with UniProt accessions (one per line)")
    sp_pdbe.add_argument("output_csv", type=Path, help="Output CSV with PDB info per UniProt accession")
    sp_pdbe.add_argument("--limit", type=int, default=10_000, help="SPARQL result limit")
    sp_pdbe.add_argument("--timeout", type=int, default=1_800, help="SPARQL timeout (seconds)")
    sp_pdbe.add_argument("--batch-size", type=int, default=10_000, help="Batch size for SPARQL queries")


def _add_search_alphafold_parser(search_subparsers: argparse._SubParsersAction):
    """Add search alphafold subcommand parser."""
    sp_af = search_subparsers.add_parser("alphafold", help="Search AlphaFold structures for UniProt accessions")
    sp_af.add_argument("uniprot_accs", type=Path, help="Text file with UniProt accessions (one per line)")
    sp_af.add_argument("output_csv", type=Path, help="Output CSV with AlphaFold IDs per UniProt accession")
    sp_af.add_argument("--limit", type=int, default=10_000, help="SPARQL result limit")
    sp_af.add_argument("--timeout", type=int, default=1_800, help="SPARQL timeout (seconds)")
    sp_af.add_argument("--batch-size", type=int, default=10_000, help="Batch size for SPARQL queries")


def _add_retrieve_pdbe_parser(retrieve_subparsers: argparse._SubParsersAction):
    """Add retrieve pdbe subcommand parser."""
    r_pdbe = retrieve_subparsers.add_parser("pdbe", help="Retrieve PDBe gzipped mmCIF files for PDB IDs in CSV")
    r_pdbe.add_argument("pdbe_csv", type=Path, help="CSV from 'search pdbe'")
    r_pdbe.add_argument("output_dir", type=Path, help="Directory to store downloaded PDBe mmCIF files")
    r_pdbe.add_argument(
        "--max-parallel-downloads",
        type=int,
        default=5,
        help="Max parallel downloads (PDBe)",
    )


def _add_retrieve_alphafold_parser(retrieve_subparsers: argparse._SubParsersAction):
    """Add retrieve alphafold subcommand parser."""
    r_af = retrieve_subparsers.add_parser("alphafold", help="Retrieve AlphaFold files for IDs in CSV")
    r_af.add_argument("alphafold_csv", type=Path, help="CSV from 'search alphafold'")
    r_af.add_argument("output_dir", type=Path, help="Directory to store downloaded AlphaFold files")
    r_af.add_argument(
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


def _add_filter_confidence_parser(filter_subparsers: argparse._SubParsersAction):
    """Add filter confidence subcommand parser."""
    f_conf = filter_subparsers.add_parser("confidence", help="Filter AlphaFold PDBs by confidence")
    f_conf.add_argument("input_dir", type=Path, help="Directory with AlphaFold PDB files")
    f_conf.add_argument("output_dir", type=Path, help="Directory to write filtered PDB files")
    f_conf.add_argument("--confidence-threshold", type=float, default=50, help="pLDDT threshold")
    f_conf.add_argument("--min-residues", type=int, default=0, help="Min high-confidence residues to keep")
    f_conf.add_argument("--max-residues", type=int, default=10_000_000, help="Max high-confidence residues to keep")


def _add_filter_chain_parser(filter_subparsers: argparse._SubParsersAction):
    """Add filter chain subcommand parser."""
    f_chain = filter_subparsers.add_parser("chain", help="Keep first UniProt chain from PDBe mmCIF")
    f_chain.add_argument("pdbe_csv", type=Path, help="CSV from 'search pdbe'")
    f_chain.add_argument("input_dir", type=Path, help="Directory with downloaded PDBe mmCIF files")
    f_chain.add_argument("output_dir", type=Path, help="Directory to write single-chain PDB files")
    f_chain.add_argument("--min-residues", type=int, default=0, help="Min residues in chain")
    f_chain.add_argument("--max-residues", type=int, default=10_000_000, help="Max residues in chain")


def _add_filter_residue_parser(filter_subparsers: argparse._SubParsersAction):
    """Add filter residue subcommand parser."""
    f_res = filter_subparsers.add_parser("residue", help="Filter PDBs by number of residues in chain A")
    f_res.add_argument("input_dir", type=Path, help="Directory with PDB files (e.g., from 'filter chain')")
    f_res.add_argument("output_dir", type=Path, help="Directory to write filtered PDB files")
    f_res.add_argument("--min-residues", type=int, default=0, help="Min residues in chain A")
    f_res.add_argument("--max-residues", type=int, default=10_000_000, help="Max residues in chain A")


def _add_search_subcommands(subparsers: argparse._SubParsersAction):
    """Add search command and its subcommands."""
    search = subparsers.add_parser("search", help="Search data sources")
    search_sp = search.add_subparsers(dest="search_cmd", required=True)

    _add_search_uniprot_parser(search_sp)
    _add_search_pdbe_parser(search_sp)
    _add_search_alphafold_parser(search_sp)


def _add_retrieve_subcommands(subparsers: argparse._SubParsersAction):
    """Add retrieve command and its subcommands."""
    retrieve = subparsers.add_parser("retrieve", help="Retrieve structure files")
    retrieve_sp = retrieve.add_subparsers(dest="retrieve_cmd", required=True)

    _add_retrieve_pdbe_parser(retrieve_sp)
    _add_retrieve_alphafold_parser(retrieve_sp)


def _add_filter_subcommands(subparsers: argparse._SubParsersAction):
    """Add filter command and its subcommands."""
    filter_cmd = subparsers.add_parser("filter", help="Filter files")
    filter_sp = filter_cmd.add_subparsers(dest="filter_cmd", required=True)

    _add_filter_confidence_parser(filter_sp)
    _add_filter_chain_parser(filter_sp)
    _add_filter_residue_parser(filter_sp)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Protein Quest CLI", prog="protein-quest")
    parser.add_argument("--log-level", default="WARNING", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_search_subcommands(subparsers)
    _add_retrieve_subcommands(subparsers)
    _add_filter_subcommands(subparsers)

    return parser


def main():
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
    reviewed: bool | None
    if args.reviewed:
        reviewed = True
    elif args.unreviewed:
        reviewed = False
    else:
        reviewed = None

    query = Query(
        taxon_id=args.taxon_id,
        reviewed=reviewed,
        subcellular_location_uniprot=args.subcellular_location_uniprot,
        subcellular_location_go=_as_scalar_or_list(args.subcellular_location_go),
        molecular_function_go=_as_scalar_or_list(args.molecular_function_go),
    )
    accs = search4uniprot(query=query, limit=args.limit, timeout=args.timeout)
    _write_lines(args.output, sorted(accs))


def _handle_search_pdbe(args):
    accs = _read_lines(args.uniprot_accs)
    results = search4pdb(accs, limit=args.limit, timeout=args.timeout, batch_size=args.batch_size)
    _write_pdbe_csv(args.output_csv, results)


def _handle_search_alphafold(args):
    accs = _read_lines(args.uniprot_accs)
    results = search4af(accs, limit=args.limit, timeout=args.timeout, batch_size=args.batch_size)
    _write_alphafold_csv(args.output_csv, results)


def _write_alphafold_summary(summary: EntrySummary, download_dir: Path):
    data = asdict(summary)
    fn = download_dir / f"{summary.entryId}.json"
    fn.write_text(json.dumps(data, indent=2))


def _handle_retrieve_pdbe(args):
    pdb_ids = _read_pdbe_ids_from_csv(args.pdbe_csv)
    pdbe_fetch.fetch(pdb_ids, args.output_dir, max_parallel_downloads=args.max_parallel_downloads)


def _handle_retrieve_alphafold(args):
    download_dir = args.output_dir
    what_af_formats = args.what_af_formats
    if what_af_formats is None:
        what_af_formats = {"pdb"}

    af_ids = [r["af_id"] for r in _read_alphafold_csv(args.alphafold_csv)]
    validated_what = structure(what_af_formats, set[DownloadableFormat])
    rprint(f"Retrieving {len(af_ids)} AlphaFold entries with formats {validated_what}")
    afs = af_fetch(af_ids, download_dir, what=validated_what, max_parallel_downloads=args.max_parallel_downloads)

    for af in tqdm(afs, unit="entry", desc="Writing summaries to disk"):
        # TODO move writing of summary.json to af_fetch function and do concurrently
        if af.summary is None:
            continue
        _write_alphafold_summary(af.summary, download_dir)

    # TODO print nr of files downloaded


def _handle_filter_confidence(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdb_files = sorted(args.input_dir.glob("*.pdb"))
    query = DensityFilterQuery(
        confidence=args.confidence_threshold,
        min_threshold=args.min_residues,
        max_threshold=args.max_residues,
    )
    for _ in tqdm(list(filter_on_density(pdb_files, query, args.output_dir)), unit="pdb"):
        # Side effects are file writes
        pass


def _handle_filter_chain(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for row in tqdm(list(_iter_pdbe_csv(args.pdbe_csv)), unit="pdb"):
        pdb_id = row["pdb_id"]
        mmcif_file = args.input_dir / f"{pdb_id.lower()}.cif.gz"
        # TODO do not filter on min/max residues just write file with single chain
        proteinpdb = ProteinPdbRow(
            id=row["pdb_id"],
            uniprot_chains=row["uniprot_chains"],
            uniprot_acc=row["uniprot_acc"],
            mmcif_file=mmcif_file if mmcif_file.exists() else None,
        )
        query = SingleChainQuery(min_residues=args.min_residues, max_residues=args.max_residues)
        write_single_chain_pdb_file(
            proteinpdb=proteinpdb,
            session_dir=args.output_dir,
            single_chain_dir=args.output_dir,
            query=query,
        )


def _handle_filter_residue(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_files = sorted(args.input_dir.glob("*.pdb"))
    kept_count = 0
    input_files = sorted(args.input_dir.glob("*.pdb"))
    for pdb in tqdm(input_files, unit="pdb"):
        n = nr_residues_in_chain(pdb, chain="A")
        if args.min_residues <= n <= args.max_residues:
            _copy_file(pdb, args.output_dir / pdb.name)
            kept_count += 1

    discarded_count = len(input_files) - kept_count
    rprint(f"Completed filtering on residues. Filtered {kept_count} and discarded {discarded_count} pdb files.")


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


def _read_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def _write_lines(path: Path, lines: Iterable[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _write_pdbe_csv(path: Path, data: dict[str, set]):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["uniprot_acc", "pdb_id", "method", "resolution", "uniprot_chains"]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
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


def _write_alphafold_csv(path: Path, data: dict[str, set[str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["uniprot_acc", "af_id"]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for uniprot_acc, af_ids in sorted(data.items()):
            for af_id in sorted(af_ids):
                writer.writerow({"uniprot_acc": uniprot_acc, "af_id": af_id})


def _read_alphafold_csv(path: Path):
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        yield from reader


def _iter_pdbe_csv(path: Path):
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        # Expect columns: uniprot_acc, pdb_id, method, resolution, uniprot_chains
        yield from reader


def _read_pdbe_ids_from_csv(path: Path) -> set[str]:
    return {row["pdb_id"] for row in _iter_pdbe_csv(path)}


def _copy_file(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    # Use buffered copy without importing shutil for minimal deps
    with src.open("rb") as fsrc, dst.open("wb") as fdst:
        while True:
            buf = fsrc.read(1024 * 1024)
            if not buf:
                break
            fdst.write(buf)
