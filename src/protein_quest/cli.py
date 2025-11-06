"""Module for cli parsers and handlers."""

import argparse
import asyncio
import csv
import logging
import os
import sys
from collections.abc import Callable, Generator, Iterable, Sequence
from contextlib import suppress
from importlib.util import find_spec
from io import BytesIO, TextIOWrapper
from pathlib import Path
from textwrap import dedent

import shtab
from cattrs import structure
from rich.console import Console
from rich.logging import RichHandler
from rich.markdown import Markdown
from rich.panel import Panel
from rich_argparse import ArgumentDefaultsRichHelpFormatter
from tqdm.rich import tqdm

from protein_quest.__version__ import __version__
from protein_quest.alphafold.confidence import ConfidenceFilterQuery, filter_files_on_confidence
from protein_quest.alphafold.fetch import DownloadableFormat, downloadable_formats
from protein_quest.alphafold.fetch import fetch_many as af_fetch
from protein_quest.converter import PositiveInt, converter
from protein_quest.emdb import fetch as emdb_fetch
from protein_quest.filters import filter_files_on_chain, filter_files_on_residues
from protein_quest.go import Aspect, allowed_aspects, search_gene_ontology_term, write_go_terms_to_csv
from protein_quest.io import (
    convert_to_cif_files,
    glob_structure_files,
    locate_structure_file,
    read_structure,
    valid_structure_file_extensions,
)
from protein_quest.pdbe import fetch as pdbe_fetch
from protein_quest.ss import SecondaryStructureFilterQuery, filter_files_on_secondary_structure
from protein_quest.structure import structure2uniprot_accessions
from protein_quest.taxonomy import SearchField, _write_taxonomy_csv, search_fields, search_taxon
from protein_quest.uniprot import (
    ComplexPortalEntry,
    PdbResults,
    Query,
    UniprotDetails,
    filter_pdb_results_on_chain_length,
    map_uniprot_accessions2uniprot_details,
    search4af,
    search4emdb,
    search4interaction_partners,
    search4macromolecular_complexes,
    search4pdb,
    search4uniprot,
)
from protein_quest.utils import (
    Cacher,
    CopyMethod,
    DirectoryCacher,
    PassthroughCacher,
    copy_methods,
    copyfile,
    user_cache_root_dir,
)

console = Console(stderr=True)
rprint = console.print
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
    ).complete = shtab.FILE
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
    parser.add_argument("--min-sequence-length", type=int, help="Minimum length of the canonical sequence.")
    parser.add_argument("--max-sequence-length", type=int, help="Maximum length of the canonical sequence.")
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
        "uniprot_accessions",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line). Use `-` for stdin.",
    ).complete = shtab.FILE
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help=dedent("""\
            Output CSV with following columns:
            `uniprot_accession`, `pdb_id`, `method`, `resolution`, `uniprot_chains`, `chain`, `chain_length`.
            Where `uniprot_chains` is the raw UniProt chain string, for example `A=1-100`.
            and where `chain` is the first chain from `uniprot_chains`, for example `A`
            and `chain_length` is the length of the chain, for example `100`.
            Use `-` for stdout.
        """),
    ).complete = shtab.FILE
    parser.add_argument(
        "--limit", type=int, default=10_000, help="Maximum number of PDB uniprot accessions combinations to return"
    )
    parser.add_argument(
        "--min-residues",
        type=int,
        help="Minimum number of residues required in the chain mapped to the UniProt accession.",
    )
    parser.add_argument(
        "--max-residues",
        type=int,
        help="Maximum number of residues allowed in chain mapped to the UniProt accession.",
    )
    parser.add_argument(
        "--keep-invalid",
        action="store_true",
        help=dedent("""\
            Keep PDB results when chain length could not be determined.
            If not given, such results are dropped.
            Only applies if min/max residues arguments are set.
        """),
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
        "uniprot_accessions",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line). Use `-` for stdin.",
    ).complete = shtab.FILE
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with AlphaFold IDs per UniProt accession. Use `-` for stdout.",
    ).complete = shtab.FILE
    parser.add_argument("--min-sequence-length", type=int, help="Minimum length of the canonical sequence.")
    parser.add_argument("--max-sequence-length", type=int, help="Maximum length of the canonical sequence.")
    parser.add_argument(
        "--limit", type=int, default=10_000, help="Maximum number of Alphafold entry identifiers to return"
    )
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")


def _add_search_emdb_parser(subparsers: argparse._SubParsersAction):
    """Add search emdb subcommand parser."""
    parser = subparsers.add_parser(
        "emdb",
        help="Search Electron Microscopy Data Bank (EMDB) identifiers of given UniProt accessions",
        description=dedent("""\
            Search for Electron Microscopy Data Bank (EMDB) identifiers of given UniProt accessions
            in the Uniprot SPARQL endpoint.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "uniprot_accs",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line). Use `-` for stdin.",
    ).complete = shtab.FILE
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with EMDB IDs per UniProt accession. Use `-` for stdout.",
    ).complete = shtab.FILE
    parser.add_argument("--limit", type=int, default=10_000, help="Maximum number of EMDB entry identifiers to return")
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
    ).complete = shtab.FILE
    parser.add_argument("--limit", type=int, default=100, help="Maximum number of GO term results to return")


def _add_search_taxonomy_parser(subparser: argparse._SubParsersAction):
    """Add search taxonomy subcommand parser."""
    parser = subparser.add_parser(
        "taxonomy",
        help="Search for taxon information in UniProt",
        description=dedent("""\
            Search for taxon information in UniProt.
            Uses https://www.uniprot.org/taxonomy?query=*.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "query", type=str, help="Search query for the taxon. Surround multiple words with quotes (' or \")."
    )
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with taxonomy results. Use `-` for stdout.",
    ).complete = shtab.FILE
    parser.add_argument(
        "--field",
        type=str,
        choices=search_fields,
        help=dedent("""\
            Field to search in. If not given then searches all fields.
            If "tax_id" then searches by taxon ID.
            If "parent" then given a parent taxon ID returns all its children.
            For example, if the parent taxon ID is 9606 (Human), it will return Neanderthal and Denisovan.
        """),
    )
    parser.add_argument("--limit", type=int, default=100, help="Maximum number of results to return")


def _add_search_interaction_partners_parser(subparsers: argparse._SubParsersAction):
    """Add search interaction partners subcommand parser."""
    parser = subparsers.add_parser(
        "interaction-partners",
        help="Search for interaction partners of given UniProt accession",
        description=dedent("""\
            Search for interaction partners of given UniProt accession
            in the Uniprot SPARQL endpoint and Complex Portal.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "uniprot_accession",
        type=str,
        help="UniProt accession (for example P12345).",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        action="append",
        help="UniProt accessions to exclude from the results. For example already known interaction partners.",
    )
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with interaction partners per UniProt accession. Use `-` for stdout.",
    ).complete = shtab.FILE
    parser.add_argument(
        "--limit", type=int, default=10_000, help="Maximum number of interaction partner uniprot accessions to return"
    )
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")


def _add_search_complexes_parser(subparsers: argparse._SubParsersAction):
    """Add search complexes subcommand parser."""
    description = dedent("""\
        Search for complexes in the Complex Portal.
        https://www.ebi.ac.uk/complexportal/

        The output CSV file has the following columns:

        - query_protein: UniProt accession used as query
        - complex_id: Complex Portal identifier
        - complex_url: URL to the Complex Portal entry
        - complex_title: Title of the complex
        - members: Semicolon-separated list of UniProt accessions of complex members
    """)
    parser = subparsers.add_parser(
        "complexes",
        help="Search for complexes in the Complex Portal",
        description=Markdown(description, style="argparse.text"),  # type: ignore using rich formatter makes this OK
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "uniprot_accessions",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line) as query for searching complexes. Use `-` for stdin.",
    ).complete = shtab.FILE
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV file with complex results. Use `-` for stdout.",
    ).complete = shtab.FILE
    parser.add_argument("--limit", type=int, default=100, help="Maximum number of complex results to return")
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")


def _add_search_uniprot_details_parser(subparsers: argparse._SubParsersAction):
    """Add search uniprot details subcommand parser."""
    description = dedent("""\
        Retrieve UniProt details for given UniProt accessions
        from the Uniprot SPARQL endpoint.

        The output CSV file has the following columns:

        - uniprot_accession: UniProt accession.
        - uniprot_id: UniProt ID (mnemonic).
        - sequence_length: Length of the canonical sequence.
        - reviewed: Whether the entry is reviewed (Swiss-Prot) or unreviewed (TrEMBL).
        - protein_name: Recommended protein name.
        - taxon_id: NCBI Taxonomy ID of the organism.
        - taxon_name: Scientific name of the organism.

        The order of the output CSV can be different from the input order.
    """)
    parser = subparsers.add_parser(
        "uniprot-details",
        help="Retrieve UniProt details for given UniProt accessions",
        description=Markdown(description, style="argparse.text"),  # type: ignore using rich formatter makes this OK
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "uniprot_accessions",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="Text file with UniProt accessions (one per line). Use `-` for stdin.",
    ).complete = shtab.FILE
    parser.add_argument(
        "output_csv",
        type=argparse.FileType("w", encoding="UTF-8"),
        help="Output CSV with UniProt details. Use `-` for stdout.",
    ).complete = shtab.FILE
    parser.add_argument("--timeout", type=int, default=1_800, help="Maximum seconds to wait for query to complete")
    parser.add_argument("--batch-size", type=int, default=1_000, help="Number of accessions to query per batch")


def _add_copy_method_arguments(parser):
    parser.add_argument(
        "--copy-method",
        type=str,
        choices=copy_methods,
        default="hardlink",
        help=dedent("""\
            How to make target file be same file as source file.
            By default uses hardlinks to save disk space.
            Note that hardlinks only work within the same filesystem and are harder to track.
            If you want to track cached files easily then use 'symlink'.
            On Windows you need developer mode or admin privileges to create symlinks.
        """),
    )


def _add_cacher_arguments(parser: argparse.ArgumentParser):
    """Add cacher arguments to parser."""
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable caching of files to central location.",
    )
    cache_dir_action = parser.add_argument(
        "--cache-dir",
        type=Path,
        default=user_cache_root_dir(),
        help="Directory to use as cache for files.",
    )
    cache_dir_action.complete = shtab.DIRECTORY  # type: ignore[missing-attribute]
    _add_copy_method_arguments(parser)


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
    ).complete = shtab.FILE
    parser.add_argument(
        "output_dir", type=Path, help="Directory to store downloaded PDBe mmCIF files"
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "--max-parallel-downloads",
        type=int,
        default=5,
        help="Maximum number of parallel downloads",
    )
    _add_cacher_arguments(parser)


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
    ).complete = shtab.FILE
    parser.add_argument(
        "output_dir", type=Path, help="Directory to store downloaded AlphaFold files"
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "--format",
        type=str,
        action="append",
        choices=sorted(downloadable_formats),
        help=dedent("""AlphaFold formats to retrieve. Can be specified multiple times.
            Default is 'cif'."""),
    )
    parser.add_argument(
        "--db-version",
        type=str,
        help="AlphaFold database version to use. If not given, the latest version is used. For example '6'.",
    )
    parser.add_argument(
        "--gzip-files",
        action="store_true",
        help="Whether to gzip the downloaded files. Excludes summary files, they are always uncompressed.",
    )
    parser.add_argument(
        "--all-isoforms",
        action="store_true",
        help=(
            "Whether to return all isoforms of each uniprot entry. "
            "If not given then only the Alphafold entry for the canonical sequence is returned."
        ),
    )
    parser.add_argument(
        "--max-parallel-downloads",
        type=int,
        default=5,
        help="Maximum number of parallel downloads",
    )
    _add_cacher_arguments(parser)


def _add_retrieve_emdb_parser(subparsers: argparse._SubParsersAction):
    """Add retrieve emdb subcommand parser."""
    parser = subparsers.add_parser(
        "emdb",
        help="Retrieve Electron Microscopy Data Bank (EMDB) gzipped 3D volume files for EMDB IDs in CSV.",
        description=dedent("""\
            Retrieve volume files from Electron Microscopy Data Bank (EMDB) website
            for unique EMDB IDs listed in a CSV file.
        """),
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "emdb_csv",
        type=argparse.FileType("r", encoding="UTF-8"),
        help="CSV file with `emdb_id` column. Other columns are ignored. Use `-` for stdin.",
    ).complete = shtab.FILE
    parser.add_argument(
        "output_dir", type=Path, help="Directory to store downloaded EMDB volume files"
    ).complete = shtab.DIRECTORY
    _add_cacher_arguments(parser)


def _add_scheduler_address_argument(parser):
    parser.add_argument(
        "--scheduler-address",
        help=dedent("""Address of the Dask scheduler to connect to.
            If not provided, will create a local cluster.
            If set to `sequential` will run tasks sequentially."""),
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
    parser.add_argument(
        "input_dir", type=Path, help="Directory with AlphaFold mmcif/PDB files"
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "output_dir", type=Path, help="Directory to write filtered mmcif/PDB files"
    ).complete = shtab.DIRECTORY
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
    ).complete = shtab.FILE
    _add_scheduler_address_argument(parser)
    _add_copy_method_arguments(parser)


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
    ).complete = shtab.FILE
    parser.add_argument(
        "input_dir",
        type=Path,
        help=dedent("""\
        Directory with PDB/mmCIF files.
        Expected filenames are `{pdb_id}.cif.gz`, `{pdb_id}.cif`, `{pdb_id}.pdb.gz` or `{pdb_id}.pdb`.
    """),
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "output_dir",
        type=Path,
        help=dedent("""\
        Directory to write the single-chain PDB/mmCIF files. Output files are in same format as input files."""),
    ).complete = shtab.DIRECTORY
    _add_scheduler_address_argument(parser)
    _add_copy_method_arguments(parser)


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
    parser.add_argument(
        "input_dir", type=Path, help="Directory with PDB/mmCIF files (e.g., from 'filter chain')"
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "output_dir",
        type=Path,
        help=dedent("""\
        Directory to write filtered PDB/mmCIF files. Files are copied without modification.
    """),
    ).complete = shtab.DIRECTORY
    parser.add_argument("--min-residues", type=int, default=0, help="Min residues in chain A")
    parser.add_argument("--max-residues", type=int, default=10_000_000, help="Max residues in chain A")
    parser.add_argument(
        "--write-stats",
        type=argparse.FileType("w", encoding="UTF-8"),
        help=dedent("""\
            Write filter statistics to file.
            In CSV format with `<input_file>,<residue_count>,<passed>,<output_file>` columns.
            Use `-` for stdout."""),
    ).complete = shtab.FILE
    _add_copy_method_arguments(parser)


def _add_filter_ss_parser(subparsers: argparse._SubParsersAction):
    """Add filter secondary structure subcommand parser."""
    parser = subparsers.add_parser(
        "secondary-structure",
        help="Filter PDB/mmCIF files by secondary structure",
        description="Filter PDB/mmCIF files by secondary structure",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "input_dir", type=Path, help="Directory with PDB/mmCIF files (e.g., from 'filter chain')"
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "output_dir",
        type=Path,
        help=dedent("""\
            Directory to write filtered PDB/mmCIF files. Files are copied without modification.
        """),
    ).complete = shtab.DIRECTORY
    parser.add_argument("--abs-min-helix-residues", type=int, help="Min residues in helices")
    parser.add_argument("--abs-max-helix-residues", type=int, help="Max residues in helices")
    parser.add_argument("--abs-min-sheet-residues", type=int, help="Min residues in sheets")
    parser.add_argument("--abs-max-sheet-residues", type=int, help="Max residues in sheets")
    parser.add_argument("--ratio-min-helix-residues", type=float, help="Min residues in helices (relative)")
    parser.add_argument("--ratio-max-helix-residues", type=float, help="Max residues in helices (relative)")
    parser.add_argument("--ratio-min-sheet-residues", type=float, help="Min residues in sheets (relative)")
    parser.add_argument("--ratio-max-sheet-residues", type=float, help="Max residues in sheets (relative)")
    parser.add_argument(
        "--write-stats",
        type=argparse.FileType("w", encoding="UTF-8"),
        help=dedent("""
            Write filter statistics to file. In CSV format with columns:
            `<input_file>,<nr_residues>,<nr_helix_residues>,<nr_sheet_residues>,
            <helix_ratio>,<sheet_ratio>,<passed>,<output_file>`.
            Use `-` for stdout.
        """),
    ).complete = shtab.FILE
    _add_copy_method_arguments(parser)


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
    _add_search_emdb_parser(subsubparsers)
    _add_search_go_parser(subsubparsers)
    _add_search_taxonomy_parser(subsubparsers)
    _add_search_interaction_partners_parser(subsubparsers)
    _add_search_complexes_parser(subsubparsers)
    _add_search_uniprot_details_parser(subsubparsers)


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
    _add_retrieve_emdb_parser(subsubparsers)


def _add_filter_subcommands(subparsers: argparse._SubParsersAction):
    """Add filter command and its subcommands."""
    parser = subparsers.add_parser("filter", help="Filter files", formatter_class=ArgumentDefaultsRichHelpFormatter)
    subsubparsers = parser.add_subparsers(dest="filter_cmd", required=True)

    _add_filter_confidence_parser(subsubparsers)
    _add_filter_chain_parser(subsubparsers)
    _add_filter_residue_parser(subsubparsers)
    _add_filter_ss_parser(subsubparsers)


def _add_convert_uniprot_parser(subparsers: argparse._SubParsersAction):
    """Add convert uniprot subcommand parser."""
    parser = subparsers.add_parser(
        "uniprot",
        help="Convert structure files to list of UniProt accessions.",
        description="Convert structure files to list of UniProt accessions. "
        "Uniprot accessions are read from database reference of each structure.",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help=f"Directory with structure files. Supported extensions are {valid_structure_file_extensions}",
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "output",
        type=argparse.FileType("wt", encoding="UTF-8"),
        help="Output text file with UniProt accessions (one per line). Use '-' for stdout.",
    ).complete = shtab.FILE
    parser.add_argument(
        "--grouped",
        action="store_true",
        help="Whether to group accessions by structure file. "
        "If set output changes to `<structure_file1>,<acc1>\\n<structure_file1>,<acc2>` format.",
    )


def _add_convert_structures_parser(subparsers: argparse._SubParsersAction):
    """Add convert structures subcommand parser."""
    parser = subparsers.add_parser(
        "structures",
        help="Convert structure files between formats",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help=f"Directory with structure files. Supported extensions are {valid_structure_file_extensions}",
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=dedent("""\
            Directory to write converted structure files. If not given, files are written to `input_dir`.
        """),
    ).complete = shtab.DIRECTORY
    parser.add_argument(
        "--format",
        type=str,
        choices=("cif",),
        default="cif",
        help="Output format to convert to.",
    )
    _add_copy_method_arguments(parser)


def _add_convert_subcommands(subparsers: argparse._SubParsersAction):
    """Add convert command and its subcommands."""
    parser = subparsers.add_parser(
        "convert",
        help="Convert files between formats",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    subsubparsers = parser.add_subparsers(dest="convert_cmd", required=True)

    _add_convert_structures_parser(subsubparsers)
    _add_convert_uniprot_parser(subsubparsers)


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
    shtab.add_argument_to(parser, ["--print-completion"])

    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_search_subcommands(subparsers)
    _add_retrieve_subcommands(subparsers)
    _add_filter_subcommands(subparsers)
    _add_convert_subcommands(subparsers)
    _add_mcp_command(subparsers)

    return parser


def _name_of(file: TextIOWrapper | BytesIO) -> str:
    try:
        return file.name
    except AttributeError:
        # In pytest BytesIO is used stdout which has no 'name' attribute
        return "<stdout>"


def _handle_search_uniprot(args):
    taxon_id = args.taxon_id
    reviewed = args.reviewed
    subcellular_location_uniprot = args.subcellular_location_uniprot
    subcellular_location_go = args.subcellular_location_go
    molecular_function_go = args.molecular_function_go
    min_sequence_length = args.min_sequence_length
    max_sequence_length = args.max_sequence_length
    limit = args.limit
    timeout = args.timeout
    output_file = args.output

    query = structure(
        {
            "taxon_id": taxon_id,
            "reviewed": reviewed,
            "subcellular_location_uniprot": subcellular_location_uniprot,
            "subcellular_location_go": subcellular_location_go,
            "molecular_function_go": molecular_function_go,
            "min_sequence_length": min_sequence_length,
            "max_sequence_length": max_sequence_length,
        },
        Query,
    )
    rprint("Searching for UniProt accessions")
    accs = search4uniprot(query=query, limit=limit, timeout=timeout)
    rprint(f"Found {len(accs)} UniProt accessions, written to {_name_of(output_file)}")
    _write_lines(output_file, sorted(accs))


def _handle_search_pdbe(args):
    uniprot_accessions = args.uniprot_accessions
    limit = args.limit
    timeout = args.timeout
    output_csv = args.output_csv
    min_residues = converter.structure(args.min_residues, PositiveInt | None)  # pyright: ignore[reportArgumentType]
    max_residues = converter.structure(args.max_residues, PositiveInt | None)  # pyright: ignore[reportArgumentType]
    keep_invalid = args.keep_invalid

    accs = set(_read_lines(uniprot_accessions))
    rprint(f"Finding PDB entries for {len(accs)} uniprot accessions")
    results = search4pdb(accs, limit=limit, timeout=timeout)

    raw_nr_results = len(results)
    raw_total_pdbs = sum([len(v) for v in results.values()])
    if min_residues or max_residues:
        results = filter_pdb_results_on_chain_length(results, min_residues, max_residues, keep_invalid=keep_invalid)
        total_pdbs = sum([len(v) for v in results.values()])
        rprint(f"Before filtering found {raw_total_pdbs} PDB entries for {raw_nr_results} uniprot accessions.")
        rprint(
            f"After filtering on chain length ({min_residues}, {max_residues}) "
            f"remained {total_pdbs} PDB entries for {len(results)} uniprot accessions."
        )
    else:
        rprint(f"Found {raw_total_pdbs} PDB entries for {raw_nr_results} uniprot accessions")

    _write_pdbe_csv(output_csv, results)
    rprint(f"Written to {_name_of(output_csv)}")


def _handle_search_alphafold(args):
    uniprot_accessions = args.uniprot_accessions
    min_sequence_length = converter.structure(args.min_sequence_length, PositiveInt | None)  # pyright: ignore[reportArgumentType]
    max_sequence_length = converter.structure(args.max_sequence_length, PositiveInt | None)  # pyright: ignore[reportArgumentType]
    limit = args.limit
    timeout = args.timeout
    output_csv = args.output_csv

    accs = _read_lines(uniprot_accessions)
    rprint(f"Finding AlphaFold entries for {len(accs)} uniprot accessions")
    results = search4af(
        accs,
        min_sequence_length=min_sequence_length,
        max_sequence_length=max_sequence_length,
        limit=limit,
        timeout=timeout,
    )
    rprint(f"Found {len(results)} AlphaFold entries, written to {_name_of(output_csv)}")
    _write_dict_of_sets2csv(output_csv, results, "af_id")


def _handle_search_emdb(args):
    uniprot_accessions = args.uniprot_accessions
    limit = args.limit
    timeout = args.timeout
    output_csv = args.output_csv

    accs = _read_lines(uniprot_accessions)
    rprint(f"Finding EMDB entries for {len(accs)} uniprot accessions")
    results = search4emdb(accs, limit=limit, timeout=timeout)
    total_emdbs = sum([len(v) for v in results.values()])
    rprint(f"Found {total_emdbs} EMDB entries, written to {_name_of(output_csv)}")
    _write_dict_of_sets2csv(output_csv, results, "emdb_id")


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
    rprint(f"Found {len(results)} GO terms, written to {_name_of(output_csv)}")
    write_go_terms_to_csv(results, output_csv)


def _handle_search_taxonomy(args):
    query: str = args.query
    field: SearchField | None = args.field
    limit: int = args.limit
    output_csv: TextIOWrapper = args.output_csv

    if field:
        rprint(f"Searching for taxon information matching '{query}' in field '{field}'")
    else:
        rprint(f"Searching for taxon information matching '{query}'")
    results = asyncio.run(search_taxon(query=query, field=field, limit=limit))
    rprint(f"Found {len(results)} taxons, written to {_name_of(output_csv)}")
    _write_taxonomy_csv(results, output_csv)


def _handle_search_interaction_partners(args: argparse.Namespace):
    uniprot_accession: str = args.uniprot_accession
    excludes: set[str] = set(args.exclude) if args.exclude else set()
    limit: int = args.limit
    timeout: int = args.timeout
    output_csv: TextIOWrapper = args.output_csv

    rprint(f"Searching for interaction partners of '{uniprot_accession}'")
    results = search4interaction_partners(uniprot_accession, excludes=excludes, limit=limit, timeout=timeout)
    rprint(f"Found {len(results)} interaction partners, written to {_name_of(output_csv)}")
    _write_lines(output_csv, results.keys())


def _handle_search_complexes(args: argparse.Namespace):
    uniprot_accessions = args.uniprot_accessions
    limit = args.limit
    timeout = args.timeout
    output_csv = args.output_csv

    accs = _read_lines(uniprot_accessions)
    rprint(f"Finding complexes for {len(accs)} uniprot accessions")
    results = search4macromolecular_complexes(accs, limit=limit, timeout=timeout)
    rprint(f"Found {len(results)} complexes, written to {_name_of(output_csv)}")
    _write_complexes_csv(results, output_csv)


def _handle_search_uniprot_details(args: argparse.Namespace):
    uniprot_accessions = args.uniprot_accessions
    timeout = args.timeout
    batch_size = args.batch_size
    output_csv: TextIOWrapper = args.output_csv

    accs = _read_lines(uniprot_accessions)
    rprint(f"Retrieving UniProt entry details for {len(accs)} uniprot accessions")
    results = list(map_uniprot_accessions2uniprot_details(accs, timeout=timeout, batch_size=batch_size))
    _write_uniprot_details_csv(output_csv, results)
    rprint(f"Retrieved details for {len(results)} UniProt entries, written to {_name_of(output_csv)}")


def _initialize_cacher(args: argparse.Namespace) -> Cacher:
    if args.no_cache:
        return PassthroughCacher()
    return DirectoryCacher(
        cache_dir=args.cache_dir,
        copy_method=args.copy_method,
    )


def _handle_retrieve_pdbe(args: argparse.Namespace):
    pdbe_csv = args.pdbe_csv
    output_dir = args.output_dir
    max_parallel_downloads = args.max_parallel_downloads
    cacher = _initialize_cacher(args)

    pdb_ids = _read_column_from_csv(pdbe_csv, "pdb_id")
    rprint(f"Retrieving {len(pdb_ids)} PDBe entries")
    result = asyncio.run(
        pdbe_fetch.fetch(pdb_ids, output_dir, max_parallel_downloads=max_parallel_downloads, cacher=cacher)
    )
    rprint(f"Retrieved {len(result)} PDBe entries")


def _handle_retrieve_alphafold(args):
    download_dir = args.output_dir
    raw_formats = args.format
    alphafold_csv = args.alphafold_csv
    max_parallel_downloads = args.max_parallel_downloads
    cacher = _initialize_cacher(args)
    gzip_files = args.gzip_files
    all_isoforms = args.all_isoforms
    db_version = args.db_version

    if raw_formats is None:
        raw_formats = {"cif"}

    # TODO besides `uniprot_accession,af_id\n` csv also allow headless single column format
    af_ids = _read_column_from_csv(alphafold_csv, "af_id")
    formats: set[DownloadableFormat] = structure(raw_formats, set[DownloadableFormat])
    rprint(f"Retrieving {len(af_ids)} AlphaFold entries with formats {formats}")
    afs = af_fetch(
        af_ids,
        download_dir,
        formats=formats,
        db_version=db_version,
        max_parallel_downloads=max_parallel_downloads,
        cacher=cacher,
        gzip_files=gzip_files,
        all_isoforms=all_isoforms,
    )
    total_nr_files = sum(af.nr_of_files() for af in afs)
    rprint(f"Retrieved {total_nr_files} AlphaFold files and {len(afs)} summaries, written to {download_dir}")


def _handle_retrieve_emdb(args):
    emdb_csv = args.emdb_csv
    output_dir = args.output_dir
    cacher = _initialize_cacher(args)

    emdb_ids = _read_column_from_csv(emdb_csv, "emdb_id")
    rprint(f"Retrieving {len(emdb_ids)} EMDB entries")
    result = asyncio.run(emdb_fetch(emdb_ids, output_dir, cacher=cacher))
    rprint(f"Retrieved {len(result)} EMDB entries")


def _handle_filter_confidence(args: argparse.Namespace):
    # we are repeating types here and in add_argument call
    # TODO replace argparse with modern alternative like cyclopts
    # to get rid of duplication
    input_dir = structure(args.input_dir, Path)
    output_dir = structure(args.output_dir, Path)

    confidence_threshold = args.confidence_threshold
    min_residues = args.min_residues
    max_residues = args.max_residues
    stats_file: TextIOWrapper | None = args.write_stats
    copy_method: CopyMethod = structure(args.copy_method, CopyMethod)  # pyright: ignore[reportArgumentType]
    scheduler_address = structure(args.scheduler_address, str | None)  # pyright: ignore[reportArgumentType]

    output_dir.mkdir(parents=True, exist_ok=True)
    input_files = sorted(glob_structure_files(input_dir))
    nr_input_files = len(input_files)
    rprint(f"Starting confidence filtering of {nr_input_files} mmcif/PDB files in {input_dir} directory.")
    query = converter.structure(
        {
            "confidence": confidence_threshold,
            "min_residues": min_residues,
            "max_residues": max_residues,
        },
        ConfidenceFilterQuery,
    )
    if stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["input_file", "residue_count", "passed", "output_file"])

    passed_count = 0
    results = filter_files_on_confidence(
        input_files, query, output_dir, copy_method=copy_method, scheduler_address=scheduler_address
    )
    for r in results:
        if r.filtered_file:
            passed_count += 1
        if stats_file:
            writer.writerow([r.input_file, r.count, r.filtered_file is not None, r.filtered_file])  # pyright: ignore[reportPossiblyUnboundVariable]
    rprint(f"Filtered {passed_count} mmcif/PDB files by confidence, written to {output_dir} directory")
    if stats_file:
        rprint(f"Statistics written to {_name_of(stats_file)}")


def _handle_filter_chain(args):
    input_dir = args.input_dir
    output_dir = structure(args.output_dir, Path)
    pdb_id2chain_mapping_file = args.chains
    scheduler_address = structure(args.scheduler_address, str | None)  # pyright: ignore[reportArgumentType]
    copy_method: CopyMethod = structure(args.copy_method, CopyMethod)  # pyright: ignore[reportArgumentType]

    # make sure files in input dir with entries in mapping file are the same
    # complain when files from mapping file are missing on disk
    rows = list(_iter_csv_rows(pdb_id2chain_mapping_file))
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
        rprint("[red]No valid structure files found. Exiting.")
        sys.exit(1)

    results = filter_files_on_chain(
        file2chain, output_dir, scheduler_address=scheduler_address, copy_method=copy_method
    )

    nr_written = len([r for r in results if r.passed])

    rprint(f"Wrote {nr_written} single-chain PDB/mmCIF files to {output_dir}.")

    for result in results:
        if result.discard_reason:
            rprint(f"[red]Discarding {result.input_file} ({result.discard_reason})[/red]")


def _handle_filter_residue(args):
    input_dir = structure(args.input_dir, Path)
    output_dir = structure(args.output_dir, Path)
    min_residues = structure(args.min_residues, int)
    max_residues = structure(args.max_residues, int)
    copy_method: CopyMethod = structure(args.copy_method, CopyMethod)  # pyright: ignore[reportArgumentType]
    stats_file: TextIOWrapper | None = args.write_stats

    if stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["input_file", "residue_count", "passed", "output_file"])

    nr_passed = 0
    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)
    rprint(f"Filtering {nr_total} files in {input_dir} directory by number of residues in chain A.")
    for r in filter_files_on_residues(
        input_files, output_dir, min_residues=min_residues, max_residues=max_residues, copy_method=copy_method
    ):
        if stats_file:
            writer.writerow([r.input_file, r.residue_count, r.passed, r.output_file])  # pyright: ignore[reportPossiblyUnboundVariable]
        if r.passed:
            nr_passed += 1

    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if stats_file:
        rprint(f"Statistics written to {_name_of(stats_file)}")


def _handle_filter_ss(args):
    input_dir = structure(args.input_dir, Path)
    output_dir = structure(args.output_dir, Path)
    copy_method: CopyMethod = structure(args.copy_method, CopyMethod)  # pyright: ignore[reportArgumentType]
    stats_file: TextIOWrapper | None = args.write_stats

    raw_query = {
        "abs_min_helix_residues": args.abs_min_helix_residues,
        "abs_max_helix_residues": args.abs_max_helix_residues,
        "abs_min_sheet_residues": args.abs_min_sheet_residues,
        "abs_max_sheet_residues": args.abs_max_sheet_residues,
        "ratio_min_helix_residues": args.ratio_min_helix_residues,
        "ratio_max_helix_residues": args.ratio_max_helix_residues,
        "ratio_min_sheet_residues": args.ratio_min_sheet_residues,
        "ratio_max_sheet_residues": args.ratio_max_sheet_residues,
    }
    query = converter.structure(raw_query, SecondaryStructureFilterQuery)
    input_files = sorted(glob_structure_files(input_dir))
    nr_total = len(input_files)
    output_dir.mkdir(parents=True, exist_ok=True)

    if stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(
            [
                "input_file",
                "nr_residues",
                "nr_helix_residues",
                "nr_sheet_residues",
                "helix_ratio",
                "sheet_ratio",
                "passed",
                "output_file",
            ]
        )

    rprint(f"Filtering {nr_total} files in {input_dir} directory by secondary structure.")
    nr_passed = 0
    for input_file, result in filter_files_on_secondary_structure(input_files, query=query):
        output_file: Path | None = None
        if result.passed:
            output_file = output_dir / input_file.name
            copyfile(input_file, output_file, copy_method)
            nr_passed += 1
        if stats_file:
            writer.writerow(  # pyright: ignore[reportPossiblyUnboundVariable]
                [
                    input_file,
                    result.stats.nr_residues,
                    result.stats.nr_helix_residues,
                    result.stats.nr_sheet_residues,
                    round(result.stats.helix_ratio, 3),
                    round(result.stats.sheet_ratio, 3),
                    result.passed,
                    output_file,
                ]
            )
    rprint(f"Wrote {nr_passed} files to {output_dir} directory.")
    if stats_file:
        rprint(f"Statistics written to {_name_of(stats_file)}")


def _handle_mcp(args):
    if find_spec("fastmcp") is None:
        msg = "Unable to start MCP server, please install `protein-quest[mcp]`."
        raise ImportError(msg)

    from protein_quest.mcp_server import mcp  # noqa: PLC0415

    if args.transport == "stdio":
        mcp.run(transport=args.transport)
    else:
        mcp.run(transport=args.transport, host=args.host, port=args.port)


def _handle_convert_uniprot(args):
    input_dir = structure(args.input_dir, Path)
    output_file: TextIOWrapper = args.output
    grouped: bool = args.grouped
    input_files = sorted(glob_structure_files(input_dir))
    if grouped:
        for input_file in tqdm(input_files, unit="file"):
            s = read_structure(input_file)
            uniprot_accessions = structure2uniprot_accessions(s)
            _write_lines(
                output_file, [f"{input_file},{uniprot_accession}" for uniprot_accession in sorted(uniprot_accessions)]
            )
    else:
        uniprot_accessions: set[str] = set()
        for input_file in tqdm(input_files, unit="file"):
            s = read_structure(input_file)
            uniprot_accessions.update(structure2uniprot_accessions(s))
        _write_lines(output_file, sorted(uniprot_accessions))


def _handle_convert_structures(args):
    input_dir = structure(args.input_dir, Path)
    output_dir = input_dir if args.output_dir is None else structure(args.output_dir, Path)
    output_dir.mkdir(parents=True, exist_ok=True)
    copy_method: CopyMethod = structure(args.copy_method, CopyMethod)  # pyright: ignore[reportArgumentType]

    input_files = sorted(glob_structure_files(input_dir))
    rprint(f"Converting {len(input_files)} files in {input_dir} directory to cif format.")
    for _ in tqdm(
        convert_to_cif_files(
            input_files,
            output_dir,
            copy_method=copy_method,
        ),
        total=len(input_files),
        unit="file",
    ):
        pass
    rprint(f"Converted {len(input_files)} files into {output_dir}.")


def _read_lines(file: TextIOWrapper) -> list[str]:
    return [line.strip() for line in file]


def _make_sure_parent_exists(file: TextIOWrapper):
    # Can not create dir for stdout
    with suppress(AttributeError):
        Path(file.name).parent.mkdir(parents=True, exist_ok=True)


def _write_lines(file: TextIOWrapper, lines: Iterable[str]):
    _make_sure_parent_exists(file)
    file.writelines(line + os.linesep for line in lines)


def _write_pdbe_csv(path: TextIOWrapper, data: PdbResults):
    _make_sure_parent_exists(path)
    fieldnames = ["uniprot_accession", "pdb_id", "method", "resolution", "uniprot_chains", "chain", "chain_length"]
    writer = csv.DictWriter(path, fieldnames=fieldnames)
    writer.writeheader()
    for uniprot_accession, entries in sorted(data.items()):
        for e in sorted(entries, key=lambda x: (x.id, x.method)):
            writer.writerow(
                {
                    "uniprot_accession": uniprot_accession,
                    "pdb_id": e.id,
                    "method": e.method,
                    "resolution": e.resolution or "",
                    "uniprot_chains": e.uniprot_chains,
                    "chain": e.chain,
                    "chain_length": e.chain_length,
                }
            )


def _write_dict_of_sets2csv(file: TextIOWrapper, data: dict[str, set[str]], ref_id_field: str):
    _make_sure_parent_exists(file)
    fieldnames = ["uniprot_accession", ref_id_field]

    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    for uniprot_accession, ref_ids in sorted(data.items()):
        for ref_id in sorted(ref_ids):
            writer.writerow({"uniprot_accession": uniprot_accession, ref_id_field: ref_id})


def _iter_csv_rows(file: TextIOWrapper) -> Generator[dict[str, str]]:
    reader = csv.DictReader(file)
    yield from reader


def _read_column_from_csv(file: TextIOWrapper, column: str) -> set[str]:
    return {row[column] for row in _iter_csv_rows(file)}


def _write_complexes_csv(complexes: list[ComplexPortalEntry], output_csv: TextIOWrapper) -> None:
    """Write ComplexPortal information to a CSV file.

    Args:
        complexes: List of ComplexPortalEntry objects.
        output_csv: TextIOWrapper to write the CSV data to.
    """
    writer = csv.writer(output_csv)
    writer.writerow(
        [
            "query_protein",
            "complex_id",
            "complex_url",
            "complex_title",
            "members",
        ]
    )
    for entry in complexes:
        members_str = ";".join(sorted(entry.members))
        writer.writerow(
            [
                entry.query_protein,
                entry.complex_id,
                entry.complex_url,
                entry.complex_title,
                members_str,
            ]
        )


def _write_uniprot_details_csv(
    output_csv: TextIOWrapper,
    uniprot_details_list: Iterable[UniprotDetails],
) -> None:
    if not uniprot_details_list:
        msg = "No UniProt entries found for given accessions"
        raise ValueError(msg)
    # As all props of UniprotDetails are scalar, we can directly unstructure to dicts
    rows = converter.unstructure(uniprot_details_list)
    fieldnames = rows[0].keys()
    writer = csv.DictWriter(output_csv, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


HANDLERS: dict[tuple[str, str | None], Callable] = {
    ("search", "uniprot"): _handle_search_uniprot,
    ("search", "pdbe"): _handle_search_pdbe,
    ("search", "alphafold"): _handle_search_alphafold,
    ("search", "emdb"): _handle_search_emdb,
    ("search", "go"): _handle_search_go,
    ("search", "taxonomy"): _handle_search_taxonomy,
    ("search", "interaction-partners"): _handle_search_interaction_partners,
    ("search", "complexes"): _handle_search_complexes,
    ("search", "uniprot-details"): _handle_search_uniprot_details,
    ("retrieve", "pdbe"): _handle_retrieve_pdbe,
    ("retrieve", "alphafold"): _handle_retrieve_alphafold,
    ("retrieve", "emdb"): _handle_retrieve_emdb,
    ("filter", "confidence"): _handle_filter_confidence,
    ("filter", "chain"): _handle_filter_chain,
    ("filter", "residue"): _handle_filter_residue,
    ("filter", "secondary-structure"): _handle_filter_ss,
    ("mcp", None): _handle_mcp,
    ("convert", "structures"): _handle_convert_structures,
    ("convert", "uniprot"): _handle_convert_uniprot,
}


def main(argv: Sequence[str] | None = None):
    """Main entry point for the CLI.

    Args:
        argv: List of command line arguments. If None, uses sys.argv.
    """
    parser = make_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=args.log_level, handlers=[RichHandler(show_level=False, console=console)])

    # Dispatch table to reduce complexity
    cmd = args.command
    sub = getattr(args, f"{cmd}_cmd", None)
    handler = HANDLERS.get((cmd, sub))
    if handler is None:
        msg = f"Unknown command: {cmd} {sub}"
        raise SystemExit(msg)
    handler(args)
