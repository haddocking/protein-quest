"""Module for fetch Alphafold data."""

import logging
from asyncio import Semaphore
from collections.abc import AsyncGenerator, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast, get_args

import aiofiles
from aiofiles.ospath import exists
from aiohttp_retry import RetryClient
from tqdm.asyncio import tqdm
from yarl import URL

from protein_quest.alphafold.entry_summary import EntrySummary
from protein_quest.converter import converter
from protein_quest.utils import Cacher, PassthroughCacher, friendly_session, retrieve_files, run_async

logger = logging.getLogger(__name__)


DownloadableFormat = Literal[
    "summary",
    "bcif",
    "cif",
    "pdb",
    "paeDoc",
    "amAnnotations",
    "amAnnotationsHg19",
    "amAnnotationsHg38",
    "msaUrl",
    "plddtDocUrl",
]
"""Types of formats that can be downloaded from the AlphaFold web service."""

downloadable_formats: set[DownloadableFormat] = set(get_args(DownloadableFormat))
"""Set of formats that can be downloaded from the AlphaFold web service."""


def _camel_to_snake_case(name: str) -> str:
    """Convert a camelCase string to snake_case."""
    return "".join(["_" + c.lower() if c.isupper() else c for c in name]).lstrip("_")


@dataclass
class AlphaFoldEntry:
    """AlphaFold entry with summary object and optionally local files.

    See https://alphafold.ebi.ac.uk/api-docs for more details on the summary data structure.
    """

    uniprot_accession: str
    summary: EntrySummary
    summary_file: Path | None = None
    bcif_file: Path | None = None
    cif_file: Path | None = None
    pdb_file: Path | None = None
    pae_doc_file: Path | None = None
    am_annotations_file: Path | None = None
    am_annotations_hg19_file: Path | None = None
    am_annotations_hg38_file: Path | None = None
    msa_file: Path | None = None
    plddt_doc_file: Path | None = None

    @classmethod
    def format2attr(cls, dl_format: DownloadableFormat) -> str:
        """Get the attribute name for a specific download format.

        Args:
            dl_format: The format for which to get the attribute name.

        Returns:
            The attribute name corresponding to the download format.

        Raises:
            ValueError: If the format is not valid.
        """
        if dl_format not in downloadable_formats:
            msg = f"Invalid format: {dl_format}. Valid formats are: {downloadable_formats}"
            raise ValueError(msg)
        return _camel_to_snake_case(dl_format) + "_file"

    def by_format(self, dl_format: DownloadableFormat) -> Path | None:
        """Get the file path for a specific format.

        Args:
            dl_format: The format for which to get the file path.

        Returns:
            The file path corresponding to the download format.
            Or None if the file is not set.

        Raises:
            ValueError: If the format is not valid.
        """
        attr = self.format2attr(dl_format)
        return getattr(self, attr, None)

    def nr_of_files(self) -> int:
        """Nr of _file properties that are set

        Returns:
            The number of _file properties that are set.
        """
        return sum(1 for attr in vars(self) if attr.endswith("_file") and getattr(self, attr) is not None)


async def fetch_summary(
    qualifier: str, session: RetryClient, semaphore: Semaphore, save_dir: Path | None, cacher: Cacher
) -> list[EntrySummary]:
    """Fetches a summary from the AlphaFold database for a given qualifier.

    Args:
        qualifier: The uniprot accession for the protein or entry to fetch.
            For example `Q5VSL9`.
        session: An asynchronous HTTP client session with retry capabilities.
        semaphore: A semaphore to limit the number of concurrent requests.
        save_dir: An optional directory to save the fetched summary as a JSON file.
            If set and summary exists then summary will be loaded from disk instead of being fetched from the API.
            If not set then the summary will not be saved to disk and will always be fetched from the API.
        cacher: A cacher to use for caching the fetched summary. Only used if save_dir is not None.

    Returns:
        A list of EntrySummary objects representing the fetched summary.
        When qualifier has multiple isoforms then multiple summaries are returned,
        otherwise a list of a single summary is returned.
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{qualifier}"
    fn: Path | None = None
    if save_dir is not None:
        fn = save_dir / f"{qualifier}.json"
        if await exists(fn):
            logger.debug(f"File {fn} already exists. Skipping download from {url}.")
            async with aiofiles.open(fn, "rb") as f:
                raw_data = await f.read()
            return converter.loads(raw_data, list[EntrySummary])
        cached_file = await cacher.copy_from_cache(Path(fn))
        if cached_file is not None:
            logger.debug(f"Using cached file {cached_file} for summary of {qualifier}.")
            async with aiofiles.open(cached_file, "rb") as f:
                raw_data = await f.read()
            return converter.loads(raw_data, list[EntrySummary])
    async with semaphore, session.get(url) as response:
        response.raise_for_status()
        raw_data = await response.content.read()
        if fn is not None:
            await cacher.write_bytes(Path(fn), raw_data)
        return converter.loads(raw_data, list[EntrySummary])


async def fetch_summaries(
    qualifiers: Iterable[str],
    save_dir: Path | None = None,
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
) -> AsyncGenerator[tuple[str, EntrySummary]]:
    semaphore = Semaphore(max_parallel_downloads)
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
    if cacher is None:
        cacher = PassthroughCacher()
    async with friendly_session() as session:
        tasks = [fetch_summary(qualifier, session, semaphore, save_dir, cacher) for qualifier in qualifiers]
        summaries_per_qualifier: list[list[EntrySummary]] = await tqdm.gather(
            *tasks, desc="Fetching Alphafold summaries"
        )
        for qualifier, summaries in zip(qualifiers, summaries_per_qualifier, strict=True):
            for summary in summaries:
                yield qualifier, summary


async def fetch_many_async(
    uniprot_accessions: Iterable[str],
    save_dir: Path,
    what: set[DownloadableFormat],
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
    gzip_files: bool = False,
    all_isoforms: bool = False,
) -> AsyncGenerator[AlphaFoldEntry]:
    """Asynchronously fetches summaries and files from
    [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).

    Args:
        uniprot_accessions: A set of Uniprot accessions to fetch.
        save_dir: The directory to save the fetched files to.
        what: A set of formats to download.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: A cacher to use for caching the fetched files. Only used if summary is in what set.
        gzip_files: Whether to gzip the downloaded files.
        all_isoforms: Whether to yield all isoforms of each uniprot entry.
            When False then yields only the canonical sequence of uniprot entry.

    Yields:
        A dataclass containing the summary, pdb file, and pae file.
    """
    save_dir_for_summaries = save_dir if "summary" in what and save_dir is not None else None

    summaries = [
        s
        async for s in fetch_summaries(
            uniprot_accessions, save_dir_for_summaries, max_parallel_downloads=max_parallel_downloads, cacher=cacher
        )
        # Filter out isoforms if all_isoforms is False
        # O60481 is canonical and O60481-2 is isoform, so we skip the isoform
        if all_isoforms or s[0] == s[1].uniprotAccession
    ]
    files = files_to_download(what, summaries, gzip_files)

    await retrieve_files(
        files,
        save_dir,
        desc="Downloading AlphaFold files",
        max_parallel_downloads=max_parallel_downloads,
        cacher=cacher,
        gzip_files=gzip_files,
    )

    gzext = ".gz" if gzip_files else ""
    for uniprot_accession, summary in summaries:
        yield AlphaFoldEntry(
            uniprot_accession=uniprot_accession,
            summary=summary,
            summary_file=save_dir / f"{uniprot_accession}.json" if save_dir_for_summaries is not None else None,
            bcif_file=save_dir / (summary.bcifUrl.name + gzext) if "bcif" in what else None,
            cif_file=save_dir / (summary.cifUrl.name + gzext) if "cif" in what else None,
            pdb_file=save_dir / (summary.pdbUrl.name + gzext) if "pdb" in what else None,
            pae_doc_file=save_dir / (summary.paeDocUrl.name + gzext) if "paeDoc" in what else None,
            am_annotations_file=(
                save_dir / (summary.amAnnotationsUrl.name + gzext)
                if "amAnnotations" in what and summary.amAnnotationsUrl
                else None
            ),
            am_annotations_hg19_file=(
                save_dir / (summary.amAnnotationsHg19Url.name + gzext)
                if "amAnnotationsHg19" in what and summary.amAnnotationsHg19Url
                else None
            ),
            am_annotations_hg38_file=(
                save_dir / (summary.amAnnotationsHg38Url.name + gzext)
                if "amAnnotationsHg38" in what and summary.amAnnotationsHg38Url
                else None
            ),
            msa_file=(save_dir / (summary.msaUrl.name + gzext) if "msaUrl" in what and summary.msaUrl else None),
            plddt_doc_file=(
                save_dir / (summary.plddtDocUrl.name + gzext) if "plddtDocUrl" in what and summary.plddtDocUrl else None
            ),
        )


def files_to_download(
    what: set[DownloadableFormat], summaries: Iterable[tuple[str, EntrySummary]], gzip_files: bool
) -> set[tuple[URL, str]]:
    if not (set(what) <= downloadable_formats):
        msg = (
            f"Invalid format(s) specified: {set(what) - downloadable_formats}. "
            f"Valid formats are: {downloadable_formats}"
        )
        raise ValueError(msg)

    url_filename_pairs: set[tuple[URL, str]] = set()
    for _, summary in summaries:
        for fmt in what:
            if fmt == "summary":
                # summary is handled already in fetch_summary
                continue
            url = cast("URL | None", getattr(summary, f"{fmt}Url", None))
            if url is None:
                logger.warning(f"Summary {summary.modelEntityId} does not have a URL for format '{fmt}'. Skipping.")
                continue
            fn = url.name + (".gz" if gzip_files else "")
            url_filename_pair = (url, fn)
            url_filename_pairs.add(url_filename_pair)
    return url_filename_pairs


def fetch_many(
    ids: Iterable[str],
    save_dir: Path,
    what: set[DownloadableFormat],
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
    gzip_files: bool = False,
    all_isoforms: bool = False,
) -> list[AlphaFoldEntry]:
    """Synchronously fetches summaries and pdb and pae files from AlphaFold Protein Structure Database.

    Args:
        ids: A set of Uniprot IDs to fetch.
        save_dir: The directory to save the fetched files to.
        what: A set of formats to download.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: A cacher to use for caching the fetched files. Only used if summary is in what set.
        gzip_files: Whether to gzip the downloaded files.
        all_isoforms: Whether to return all isoforms of each uniprot entry.
            When False then returns only the canonical sequence of uniprot entry.

    Returns:
        A list of AlphaFoldEntry dataclasses containing the summary, pdb file, and pae file.
    """

    async def gather_entries():
        return [
            entry
            async for entry in fetch_many_async(
                ids,
                save_dir,
                what,
                max_parallel_downloads=max_parallel_downloads,
                cacher=cacher,
                gzip_files=gzip_files,
                all_isoforms=all_isoforms,
            )
        ]

    return run_async(gather_entries())


def relative_to(entry: AlphaFoldEntry, session_dir: Path) -> AlphaFoldEntry:
    """Convert paths in an AlphaFoldEntry to be relative to the session directory.

    Args:
        entry: An AlphaFoldEntry instance with absolute paths.
        session_dir: The session directory to which the paths should be made relative.

    Returns:
        An AlphaFoldEntry instance with paths relative to the session directory.
    """
    return AlphaFoldEntry(
        uniprot_accession=entry.uniprot_accession,
        summary=entry.summary,
        summary_file=entry.summary_file.relative_to(session_dir) if entry.summary_file else None,
        bcif_file=entry.bcif_file.relative_to(session_dir) if entry.bcif_file else None,
        cif_file=entry.cif_file.relative_to(session_dir) if entry.cif_file else None,
        pdb_file=entry.pdb_file.relative_to(session_dir) if entry.pdb_file else None,
        pae_doc_file=entry.pae_doc_file.relative_to(session_dir) if entry.pae_doc_file else None,
        am_annotations_file=entry.am_annotations_file.relative_to(session_dir) if entry.am_annotations_file else None,
        am_annotations_hg19_file=(
            entry.am_annotations_hg19_file.relative_to(session_dir) if entry.am_annotations_hg19_file else None
        ),
        am_annotations_hg38_file=(
            entry.am_annotations_hg38_file.relative_to(session_dir) if entry.am_annotations_hg38_file else None
        ),
        msa_file=entry.msa_file.relative_to(session_dir) if entry.msa_file else None,
        plddt_doc_file=entry.plddt_doc_file.relative_to(session_dir) if entry.plddt_doc_file else None,
    )
