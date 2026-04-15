"""Search subcommands for protein-quest."""

import asyncio
import contextlib
import csv
import logging
import os
from collections.abc import Iterable, Mapping, Sequence
from typing import Annotated, Any, Literal, cast

from cyclopts import App, Parameter
from cyclopts.types import PositiveInt, StdioPath

from protein_quest.cli.common import (
    BatchSize,
    Common,
    InputFile,
    Limit,
    MaxResidues,
    MaxSequenceLength,
    MinResidues,
    MinSequenceLength,
    OutputFile,
    Timeout,
    console,
)
from protein_quest.converter import converter
from protein_quest.go import Aspect, GoTerm, search_gene_ontology_term
from protein_quest.pdbe_3dbeacons.model import Provider, search_structure_provider_choices
from protein_quest.pdbe_3dbeacons.search import PruneOptions, flatten_structure_summaries, uniprots2structures
from protein_quest.taxonomy import SearchField, Taxon, search_taxon
from protein_quest.uniprot import (
    ComplexPortalEntry,
    PdbChainLengthError,
    PdbResults,
    Query,
    UniprotDetails,
    filter_pdb_results_on_chain_length,
    filter_pdb_results_on_resolution,
    map_uniprot_accessions2uniprot_details,
    search4af,
    search4emdb,
    search4interaction_partners,
    search4macromolecular_complexes,
    search4pdb,
    search4uniprot,
)

rprint = console.print
logger = logging.getLogger(__name__)

search_app = App(name="search", help="Search data sources")


def _name_of_path(file: StdioPath) -> str:
    """Return a display name for a StdioPath."""
    if str(file) == "-":
        return "<stdout>"
    return str(file)


def _write_pdbe_csv(path: StdioPath, data: PdbResults):
    """Write PDBe results to CSV."""
    fieldnames = ("uniprot_accession", "pdb_id", "method", "resolution", "uniprot_chains", "chain", "chain_length")
    with path.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for uniprot_accession, entries in sorted(data.items()):
            for e in sorted(entries, key=lambda x: (x.id, x.method)):
                try:
                    chain_length = e.chain_length
                except PdbChainLengthError:
                    msg = f"Could not determine chain length for {uniprot_accession} / {e.id} chain {e.chain}"
                    msg += f" from '{e.uniprot_chains}'. No chain length for this entry."
                    logger.warning(msg)
                    chain_length = None
                writer.writerow(
                    {
                        "uniprot_accession": uniprot_accession,
                        "pdb_id": e.id,
                        "method": e.method,
                        "resolution": e.resolution or "",
                        "uniprot_chains": e.uniprot_chains,
                        "chain": e.chain,
                        "chain_length": chain_length,
                    }
                )


def _write_dict_of_sets2csv(file: StdioPath, data: dict[str, set[str]], ref_id_field: str):
    """Write a dict of sets to CSV."""
    if str(file) != "-":
        file.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["uniprot_accession", ref_id_field]

    with file.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for uniprot_accession, ref_ids in sorted(data.items()):
            for ref_id in sorted(ref_ids):
                writer.writerow({"uniprot_accession": uniprot_accession, ref_id_field: ref_id})


def _write_lines(file: StdioPath, lines: Iterable[str]):
    """Write lines to a file or stdout."""
    if str(file) != "-":
        file.parent.mkdir(parents=True, exist_ok=True)
    with file.open("w", encoding="utf-8") as f:
        for line in lines:
            f.write(line + os.linesep)


def _write_list_of_dicts_to_csv(file: StdioPath, rows: Sequence[Mapping[str, Any]]):
    """Write a list of dicts to CSV."""
    if str(file) != "-":
        file.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = rows[0].keys()

    with file.open("w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


@contextlib.contextmanager
def write_csv(path: StdioPath):
    """Context manager for writing CSV files.

    Creates parent directories if they do not exist.

    Yields:
        CSV writer object to write rows to.
    """
    if str(path) != "-":
        path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        yield csv.writer(f)


def _write_complexes_csv(complexes: list[ComplexPortalEntry], output_csv: StdioPath):
    """Write ComplexPortal information to a CSV file."""
    with write_csv(output_csv) as writer:
        writer.writerow(("query_protein", "complex_id", "complex_url", "complex_title", "members"))
        for entry in complexes:
            members_str = ";".join(sorted(entry.members))
            writer.writerow(
                (entry.query_protein, entry.complex_id, entry.complex_url, entry.complex_title, members_str)
            )


def _write_uniprot_details_csv(output_csv: StdioPath, uniprot_details_list: list[UniprotDetails]):
    """Write UniProt details to CSV."""
    if not uniprot_details_list:
        msg = "No UniProt entries found for given accessions"
        raise ValueError(msg)
    _write_list_of_dicts_to_csv(output_csv, uniprot_details_list)


def _read_lines(file: StdioPath) -> list[str]:
    """Read lines from a file or stdin."""
    with file.open("r", encoding="utf-8") as f:
        return [line.strip() for line in f]


def _write_go_terms_to_csv(terms: list[GoTerm], csv_file: StdioPath) -> None:
    """Write a list of GO terms to a CSV file.

    Args:
        terms: The list of GO terms to write.
        csv_file: The CSV file to write to. Use `-` for stdout.
    """
    with write_csv(csv_file) as writer:
        writer.writerow(["id", "name", "aspect", "definition"])
        for term in terms:
            writer.writerow([term.id, term.name, term.aspect, term.definition])


def write_taxonomy_csv(taxons: list[Taxon], output_csv: StdioPath) -> None:
    """Write taxon information to a CSV file.

    Args:
        taxons: List of Taxon objects to write to the CSV file.
        output_csv: File object for the output CSV file. Can be a file path or '-' for stdout.
    """
    with write_csv(output_csv) as writer:
        writer.writerow(["taxon_id", "scientific_name", "common_name", "rank", "other_names"])
        for taxon in taxons:
            writer.writerow(
                [
                    taxon.taxon_id,
                    taxon.scientific_name,
                    taxon.common_name,
                    taxon.rank,
                    ";".join(taxon.other_names) if taxon.other_names else "",
                ]
            )


@search_app.command
def uniprot(
    output: OutputFile,
    /,
    *,
    query: Query | None = None,
    limit: Limit = 10_000,
    timeout: Timeout = 1_800,
    _: Common | None = None,
) -> None:
    """Search for UniProt accessions.

    Search for UniProt accessions based on various criteria in the Uniprot SPARQL endpoint.

    Args:
        output: Output text file for UniProt accessions (one per line). Use `-` for stdout.
        query: Search query for UniProtKB. Can be given as a JSON string or as a path to a JSON file.
        limit: Maximum number of uniprot accessions to return.
        timeout: Maximum seconds to wait for query to complete.
        _: Common CLI options.
    """
    if query is None:
        query = Query()

    rprint("Searching for UniProt accessions")
    accs = search4uniprot(query=query, limit=limit, timeout=timeout)
    rprint(f"Found {len(accs)} UniProt accessions, written to {_name_of_path(output)}")
    _write_lines(output, sorted(accs))


@search_app.command
def pdbe(
    uniprot_accessions: InputFile,
    output_csv: OutputFile,
    /,
    *,
    limit: Limit = 10_000,
    timeout: Timeout = 1_800,
    min_residues: MinResidues | None = None,
    max_residues: MaxResidues | None = None,
    keep_invalid: Annotated[bool, Parameter(negative="")] = False,
    top_resolution_per_uniprot_accession: PositiveInt | None = None,
    _: Common | None = None,
) -> None:
    """Search for PDB structures of given UniProt accessions.

    Search for PDB structures of given UniProt accessions in the Uniprot SPARQL endpoint.

    Args:
        uniprot_accessions: Text file with UniProt accessions (one per line). Use `-` for stdin.
        output_csv: Output CSV with following columns:
            `uniprot_accession`, `pdb_id`, `method`, `resolution`, `uniprot_chains`, `chain`, `chain_length`.
            Where `uniprot_chains` is the raw UniProt chain string, for example `A=1-100`.
            And where `chain` is the first chain from `uniprot_chains`, for example `A`.
            And `chain_length` is the length of the chain, for example `100`
            or '' if it could not be determined.
            Use `-` for stdout.
        limit: Maximum number of PDB uniprot accessions combinations to return.
        timeout: Maximum seconds to wait for query to complete.
        min_residues: Minimum number of residues required in the chain mapped to the UniProt accession.
        max_residues: Maximum number of residues allowed in chain mapped to the UniProt accession.
        keep_invalid: Keep PDB results when chain length could not be determined.
        top_resolution_per_uniprot_accession: Retain the top N PDB entries per UniProt accession,
            ranked by best (lowest) resolution first, then by highest residue count.
        _: Common CLI options.
    """
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

    if top_resolution_per_uniprot_accession is not None:
        results = filter_pdb_results_on_resolution(results, top=top_resolution_per_uniprot_accession)
        total_pdbs = sum([len(v) for v in results.values()])
        rprint(
            f"After filtering on top resolution ({top_resolution_per_uniprot_accession}) "
            f"remained {total_pdbs} PDB entries for {len(results)} uniprot accessions."
        )

    _write_pdbe_csv(output_csv, results)
    rprint(f"Written to {output_csv}")


@search_app.command
def alphafold(
    uniprot_accessions: InputFile,
    output_csv: OutputFile,
    /,
    *,
    min_sequence_length: MinSequenceLength | None = None,
    max_sequence_length: MaxSequenceLength | None = None,
    limit: Limit = 10_000,
    timeout: Timeout = 1_800,
    _: Common | None = None,
) -> None:
    """Search for AlphaFold structures of given UniProt accessions.

    Search for AlphaFold structures of given UniProt accessions in the Uniprot SPARQL endpoint.

    Args:
        uniprot_accessions: Text file with UniProt accessions (one per line). Use `-` for stdin.
        output_csv: Output CSV with AlphaFold IDs per UniProt accession.
            CSV has columns: `uniprot_accession`, `af_id`.
            Use `-` for stdout.
        min_sequence_length: Minimum length of the canonical sequence.
        max_sequence_length: Maximum length of the canonical sequence.
        limit: Maximum number of Alphafold entry identifiers to return.
        timeout: Maximum seconds to wait for query to complete.
        _: Common CLI options.
    """
    accs = _read_lines(uniprot_accessions)
    rprint(f"Finding AlphaFold entries for {len(accs)} uniprot accessions")
    results = search4af(
        accs,
        min_sequence_length=min_sequence_length,
        max_sequence_length=max_sequence_length,
        limit=limit,
        timeout=timeout,
    )
    rprint(f"Found {len(results)} AlphaFold entries, written to {output_csv}")
    _write_dict_of_sets2csv(output_csv, results, "af_id")


@search_app.command
def structure(
    uniprot_accessions: InputFile,
    output_csv: OutputFile,
    /,
    *,
    source: Annotated[
        set[Provider | Literal["all"]] | None,
        Parameter(negative=""),
    ] = None,
    min_residues: MinResidues | None = None,
    max_residues: MaxResidues | None = None,
    limit: Limit = 10_000,
    timeout: Timeout = 1_800,
    raw: OutputFile | None = None,
    _: Common | None = None,
) -> None:
    """Search for experimentally determined and predicted structures.

    Search for experimentally determined and predicted structures of given UniProt accessions
    in the 3D Beacons Network API.

    Args:
        uniprot_accessions: Text file with UniProt accessions (one per line). Use `-` for stdin.
        output_csv: Output CSV with following columns:
            `uniprot_accession`, `provider`, `model_identifier`, `model_url`, `model_format`, `chain`, `residue_count`.
            Use `-` for stdout.
        source: Source of the structures to search for. Default `pdbe` and `alphafold`.
            Multiple sources can be given by repeating the `--source` parameter.
            Use 'all' to search all sources.
        min_residues: Minimum number of residues required in the chain mapped to the UniProt accession.
        max_residues: Maximum number of residues allowed in the chain mapped to the UniProt accession.
        limit: Maximum number of structures per uniprot accession per source to return.
        timeout: Maximum seconds to wait for query to complete.
        raw: Path to write raw 3D beacon summaries as JSON.
        _: Common CLI options.
    """
    nsource: set[Provider]
    if source is None:
        nsource = {"pdbe", "alphafold"}
    elif "all" in source:
        nsource = set(search_structure_provider_choices)
    else:
        nsource = cast("set[Provider]", source)

    accs = set(_read_lines(uniprot_accessions))
    rprint(f"Finding structures for {len(accs)} uniprot accessions")

    prune_options = PruneOptions(
        providers=nsource,
        limit=limit,
        min_residues=min_residues,
        max_residues=max_residues,
    )

    results = asyncio.run(
        uniprots2structures(
            accs,
            prune_options,
            timeout=timeout,
        )
    )

    if raw:
        raw.write_bytes(converter.dumps(results))
        rprint(f"Written raw results to {raw}")

    rows = flatten_structure_summaries(results)
    if not rows:
        rprint("No structures found")
        return

    _write_list_of_dicts_to_csv(output_csv, rows)
    rprint(f"Found {len(rows)} structures, written to {output_csv}")


@search_app.command
def emdb(
    uniprot_accessions: InputFile,
    output_csv: OutputFile,
    /,
    *,
    limit: Limit = 10_000,
    timeout: Timeout = 1_800,
    _: Common | None = None,
) -> None:
    """Search for EMDB identifiers of given UniProt accessions.

    Search for Electron Microscopy Data Bank (EMDB) identifiers of given UniProt accessions
    in the Uniprot SPARQL endpoint.

    Args:
        uniprot_accessions: Text file with UniProt accessions (one per line). Use `-` for stdin.
        output_csv: Output CSV with EMDB IDs per UniProt accession.
            CSV has columns: `uniprot_accession`, `emdb_id`. Use `-` for stdout.
        limit: Maximum number of EMDB entry identifiers to return.
        timeout: Maximum seconds to wait for query to complete.
        _: Common CLI options.
    """
    accs = _read_lines(uniprot_accessions)
    rprint(f"Finding EMDB entries for {len(accs)} uniprot accessions")
    results = search4emdb(accs, limit=limit, timeout=timeout)
    total_emdbs = sum([len(v) for v in results.values()])
    rprint(f"Found {total_emdbs} EMDB entries, written to {output_csv}")
    _write_dict_of_sets2csv(output_csv, results, "emdb_id")


@search_app.command
def go(
    term: str,
    output_csv: OutputFile,
    /,
    *,
    aspect: Aspect | None = None,
    limit: Limit = 100,
    _: Common | None = None,
) -> None:
    """Search for Gene Ontology (GO) terms.

    Search for Gene Ontology (GO) terms in the EBI QuickGO API.

    Args:
        term: GO term to search for. For example `apoptosome`.
        output_csv: Output CSV with GO term results.
            CSV has columns: `term`, `id`, `name`, `aspect`, `definition`.
            Use `-` for stdout.
        aspect: Filter on aspect.
        limit: Maximum number of GO term results to return.
        _: Common CLI options.
    """
    if aspect:
        rprint(f"Searching for GO terms matching '{term}' with aspect '{aspect}'")
    else:
        rprint(f"Searching for GO terms matching '{term}'")
    results = asyncio.run(search_gene_ontology_term(term, aspect=aspect, limit=limit))
    rprint(f"Found {len(results)} GO terms, written to {output_csv}")
    _write_go_terms_to_csv(results, output_csv)


@search_app.command
def taxonomy(
    query: str,
    output_csv: OutputFile,
    /,
    *,
    field: SearchField | None = None,
    limit: Limit = 100,
    _: Common | None = None,
) -> None:
    """Search for taxon information in UniProt.

    Search for taxon information in UniProt. Uses https://www.uniprot.org/taxonomy?query=*.

    Args:
        query: Search query for the taxon. Surround multiple words with quotes.
        output_csv: Output CSV with taxonomy results.
            CSV has columns: `tax_id`, `name`, `rank`, `parent_tax_id`, `parent_tax_name`.
            Use `-` for stdout.
        field: Field to search in. If not given then searches all fields.
            If "tax_id" then searches by taxon ID.
            If "parent" then given a parent taxon ID returns all its children.
        limit: Maximum number of results to return.
        _: Common CLI options.
    """
    if field:
        rprint(f"Searching for taxon information matching '{query}' in field '{field}'")
    else:
        rprint(f"Searching for taxon information matching '{query}'")
    results = asyncio.run(search_taxon(query=query, field=field, limit=limit))
    rprint(f"Found {len(results)} taxons, written to {output_csv}")
    write_taxonomy_csv(results, output_csv)


@search_app.command
def interaction_partners(
    uniprot_accession: str,
    output_csv: OutputFile,
    /,
    *,
    exclude: Annotated[list[str] | None, Parameter(negative="")] = None,
    limit: Limit = 10_000,
    timeout: Timeout = 1_800,
    _: Common | None = None,
) -> None:
    """Search for interaction partners of given UniProt accession.

    Search for interaction partners of given UniProt accession in the Uniprot SPARQL endpoint
    and Complex Portal.

    Args:
        uniprot_accession: UniProt accession (for example P12345).
        output_csv: Output CSV with interaction partners per UniProt accession.
            CSV has columns: `uniprot_accession`. Use `-` for stdout.
        exclude: UniProt accessions to exclude from the results.
            Multiple accessions can be given by repeating the `--exclude` option.
        limit: Maximum number of interaction partner uniprot accessions to return.
        timeout: Maximum seconds to wait for query to complete.
        _: Common CLI options.
    """
    excludes = set(exclude) if exclude else set()
    rprint(f"Searching for interaction partners of '{uniprot_accession}'")
    results = search4interaction_partners(uniprot_accession, excludes=excludes, limit=limit, timeout=timeout)
    rprint(f"Found {len(results)} interaction partners, written to {output_csv}")
    _write_lines(output_csv, results.keys())


@search_app.command
def complexes(
    uniprot_accessions: InputFile,
    output_csv: OutputFile,
    /,
    *,
    limit: Limit = 100,
    timeout: Timeout = 1_800,
    _: Common | None = None,
) -> None:
    """Search for complexes in the Complex Portal.

    Search for complexes in the Complex Portal (https://www.ebi.ac.uk/complexportal/).

    The output CSV file has the following columns:

    - query_protein: UniProt accession used as query
    - complex_id: Complex Portal identifier
    - complex_url: URL to the Complex Portal entry
    - complex_title: Title of the complex
    - members: Semicolon-separated list of UniProt accessions of complex members

    Args:
        uniprot_accessions: Text file with UniProt accessions (one per line) as query. Use `-` for stdin.
        output_csv: Output CSV file with complex results. Use `-` for stdout.
        limit: Maximum number of complex results to return.
        timeout: Maximum seconds to wait for query to complete.
        _: Common CLI options.
    """
    accs = _read_lines(uniprot_accessions)
    rprint(f"Finding complexes for {len(accs)} uniprot accessions")
    results = search4macromolecular_complexes(accs, limit=limit, timeout=timeout)
    rprint(f"Found {len(results)} complexes, written to {output_csv}")
    _write_complexes_csv(results, output_csv)


@search_app.command
def uniprot_details(
    uniprot_accessions: InputFile,
    output_csv: OutputFile,
    /,
    *,
    timeout: Timeout = 1_800,
    batch_size: BatchSize = 1_000,
    _: Common | None = None,
) -> None:
    """Search for UniProt details for given UniProt accessions from the UniProt SPARQL endpoint.

    The output CSV file has the following columns:

    - uniprot_accession: UniProt accession.
    - uniprot_id: UniProt ID (mnemonic).
    - sequence_length: Length of the canonical sequence.
    - reviewed: Whether the entry is reviewed (Swiss-Prot) or unreviewed (TrEMBL).
    - protein_name: Recommended protein name.
    - taxon_id: NCBI Taxonomy ID of the organism.
    - taxon_name: Scientific name of the organism.

    The order of the output CSV can be different from the input order.

    Args:
        uniprot_accessions: Text file with UniProt accessions (one per line). Use `-` for stdin.
        output_csv: Output CSV with UniProt details.
            CSV has columns: `uniprot_accession`, `uniprot_id`, `sequence_length`, `reviewed`,
            `protein_name`, `taxon_id`, `taxon_name`.
            Use `-` for stdout.
        timeout: Maximum seconds to wait for query to complete.
        batch_size: Number of accessions to query per batch.
        _: Common CLI options.
    """
    accs = _read_lines(uniprot_accessions)
    rprint(f"Retrieving UniProt entry details for {len(accs)} uniprot accessions")
    results = list(map_uniprot_accessions2uniprot_details(accs, timeout=timeout, batch_size=batch_size))
    _write_uniprot_details_csv(output_csv, results)
    rprint(f"Retrieved details for {len(results)} UniProt entries, written to {output_csv}")
