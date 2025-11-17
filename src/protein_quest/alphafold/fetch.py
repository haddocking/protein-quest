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
    "msa",
    "plddtDoc",
]
"""Types of formats that can be downloaded from the AlphaFold web service."""

downloadable_formats: set[DownloadableFormat] = set(get_args(DownloadableFormat))
"""Set of formats that can be downloaded from the AlphaFold web service."""

UrlFileNamePair = tuple[URL, str]
"""A tuple of a URL and a filename."""
UrlFileNamePairsOfFormats = dict[DownloadableFormat, UrlFileNamePair]
"""A mapping of DownloadableFormat to UrlFileNamePair."""


def _camel_to_snake_case(name: str) -> str:
    """Convert a camelCase string to snake_case."""
    return "".join(["_" + c.lower() if c.isupper() else c for c in name]).lstrip("_")


@dataclass
class AlphaFoldEntry:
    """AlphaFold entry with summary object and optionally local files.

    See https://alphafold.ebi.ac.uk/api-docs for more details on the summary data structure.
    """

    uniprot_accession: str
    summary: EntrySummary | None = None
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

    def relative_to(self, session_dir: Path) -> "AlphaFoldEntry":
        """Convert paths in an AlphaFoldEntry to be relative to the session directory.

        Args:
            session_dir: The session directory to which the paths should be made relative.

        Returns:
            An AlphaFoldEntry instance with paths relative to the session directory.
        """
        return AlphaFoldEntry(
            uniprot_accession=self.uniprot_accession,
            summary=self.summary,
            summary_file=self.summary_file.relative_to(session_dir) if self.summary_file else None,
            bcif_file=self.bcif_file.relative_to(session_dir) if self.bcif_file else None,
            cif_file=self.cif_file.relative_to(session_dir) if self.cif_file else None,
            pdb_file=self.pdb_file.relative_to(session_dir) if self.pdb_file else None,
            pae_doc_file=self.pae_doc_file.relative_to(session_dir) if self.pae_doc_file else None,
            am_annotations_file=self.am_annotations_file.relative_to(session_dir) if self.am_annotations_file else None,
            am_annotations_hg19_file=(
                self.am_annotations_hg19_file.relative_to(session_dir) if self.am_annotations_hg19_file else None
            ),
            am_annotations_hg38_file=(
                self.am_annotations_hg38_file.relative_to(session_dir) if self.am_annotations_hg38_file else None
            ),
            msa_file=self.msa_file.relative_to(session_dir) if self.msa_file else None,
            plddt_doc_file=self.plddt_doc_file.relative_to(session_dir) if self.plddt_doc_file else None,
        )


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


async def _fetch_many_async_with_summary(
    uniprot_accessions: Iterable[str],
    save_dir: Path,
    formats: set[DownloadableFormat],
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
    gzip_files: bool = False,
    all_isoforms: bool = False,
) -> AsyncGenerator[AlphaFoldEntry]:
    save_dir_for_summaries = save_dir if "summary" in formats else None

    summaries = [
        s
        async for s in fetch_summaries(
            uniprot_accessions, save_dir_for_summaries, max_parallel_downloads=max_parallel_downloads, cacher=cacher
        )
        # Filter out isoforms if all_isoforms is False
        # O60481 is canonical and O60481-2 is isoform, so we skip the isoform
        if all_isoforms or s[0] == s[1].uniprotAccession
    ]
    files = files_to_download(formats, summaries, gzip_files)

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
            bcif_file=save_dir / (summary.bcifUrl.name + gzext) if "bcif" in formats else None,
            cif_file=save_dir / (summary.cifUrl.name + gzext) if "cif" in formats else None,
            pdb_file=save_dir / (summary.pdbUrl.name + gzext) if "pdb" in formats else None,
            pae_doc_file=save_dir / (summary.paeDocUrl.name + gzext) if "paeDoc" in formats else None,
            am_annotations_file=(
                save_dir / (summary.amAnnotationsUrl.name + gzext)
                if "amAnnotations" in formats and summary.amAnnotationsUrl
                else None
            ),
            am_annotations_hg19_file=(
                save_dir / (summary.amAnnotationsHg19Url.name + gzext)
                if "amAnnotationsHg19" in formats and summary.amAnnotationsHg19Url
                else None
            ),
            am_annotations_hg38_file=(
                save_dir / (summary.amAnnotationsHg38Url.name + gzext)
                if "amAnnotationsHg38" in formats and summary.amAnnotationsHg38Url
                else None
            ),
            msa_file=(save_dir / (summary.msaUrl.name + gzext) if "msa" in formats and summary.msaUrl else None),
            plddt_doc_file=(
                save_dir / (summary.plddtDocUrl.name + gzext) if "plddtDoc" in formats and summary.plddtDocUrl else None
            ),
        )


def files_to_download(
    formats: set[DownloadableFormat], summaries: Iterable[tuple[str, EntrySummary]], gzip_files: bool
) -> set[UrlFileNamePair]:
    if not (set(formats) <= downloadable_formats):
        msg = (
            f"Invalid format(s) specified: {set(formats) - downloadable_formats}. "
            f"Valid formats are: {downloadable_formats}"
        )
        raise ValueError(msg)

    url_filename_pairs: set[UrlFileNamePair] = set()
    for _, summary in summaries:
        for fmt in formats:
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


async def fetch_alphafold_db_version() -> str:
    """Fetch the current version of the AlphaFold database.

    Returns:
        The current version of the AlphaFold database as a string. For example: "6".
    """
    url = "https://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv"
    headers = {"Range": "bytes=0-200"}
    logger.debug(f"Detecting AlphaFold DB version from head of {url}")
    async with friendly_session() as session, session.get(url, headers=headers) as response:
        response.raise_for_status()
        raw = await response.content.read(200)
        text = raw.decode("utf-8")
        first_line = text.splitlines()[1]
        version = first_line.split(",")[-1]
        logger.debug(f"Found current AlphaFold DB version is '{version}'")
        return version


def _files_for_alphafold_entry(
    uniprot_accession: str,
    formats: set[DownloadableFormat],
    db_version: str,
    gzip_files: bool,
) -> UrlFileNamePairsOfFormats:
    templates: dict[DownloadableFormat, URL] = {
        "bcif": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-model_v{db_version}.bcif"),
        "cif": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-model_v{db_version}.cif"),
        "pdb": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-model_v{db_version}.pdb"),
        "paeDoc": URL(
            f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-predicted_aligned_error_v{db_version}.json"
        ),
        "amAnnotations": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-aa-substitutions.csv"),
        "amAnnotationsHg19": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-hg19.csv"),
        "amAnnotationsHg38": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-hg38.csv"),
        "msa": URL(f"https://alphafold.ebi.ac.uk/files/msa/AF-{uniprot_accession}-F1-msa_v{db_version}.a3m"),
        "plddtDoc": URL(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-confidence_v{db_version}.json"),
    }
    url_filename_pairs = {}
    for fmt in formats:
        if fmt == "summary":
            # Summaries are downloaded separately as its using API instead of static files
            continue
        if fmt not in templates:
            logger.warning(f"No URL template found for format '{fmt}'. Skipping.")
            continue
        url = templates[cast("DownloadableFormat", fmt)]
        fn = url.name
        if gzip_files:
            fn += ".gz"
        url_filename_pair = (url, fn)
        url_filename_pairs[fmt] = url_filename_pair
    return url_filename_pairs


def files_for_alphafold_entries(
    uniprot_accessions: Iterable[str],
    formats: set[DownloadableFormat],
    db_version: str,
    gzip_files: bool,
) -> dict[str, UrlFileNamePairsOfFormats]:
    """Get the files to download for multiple AlphaFold entries.

    Args:
        uniprot_accessions: A set of Uniprot accessions.
        formats: A set of formats to download.
        db_version: The version of the AlphaFold database to use.
        gzip_files: Whether to download gzipped files. Otherwise downloads uncompressed files.

    Returns:
        A mapping of Uniprot accession to a mapping of DownloadableFormat to UrlFileNamePair.
    """
    return {
        uniprot_accession: _files_for_alphafold_entry(
            uniprot_accession, formats=formats, db_version=db_version, gzip_files=gzip_files
        )
        for uniprot_accession in uniprot_accessions
    }


async def _fetch_many_async_without_summary(
    uniprot_accessions: Iterable[str],
    save_dir: Path,
    formats: set[DownloadableFormat],
    db_version: str | None = None,
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
    gzip_files: bool = False,
) -> AsyncGenerator[AlphaFoldEntry]:
    if db_version is None:
        db_version = await fetch_alphafold_db_version()
    nested_files = files_for_alphafold_entries(
        uniprot_accessions, formats=formats, db_version=db_version, gzip_files=gzip_files
    )
    files: set[UrlFileNamePair] = set()
    for uniprot_accession in uniprot_accessions:
        files.update(nested_files[uniprot_accession].values())

    retrieved_files = await retrieve_files(
        files,
        save_dir,
        desc="Downloading AlphaFold files",
        max_parallel_downloads=max_parallel_downloads,
        cacher=cacher,
        gzip_files=gzip_files,
        raise_for_not_found=False,
    )

    retrieved_files_set = set(retrieved_files)
    for uniprot_accession in uniprot_accessions:
        entry = AlphaFoldEntry(
            uniprot_accession=uniprot_accession,
        )

        for af_format, url_filename_pair in nested_files[uniprot_accession].items():
            _, filename = url_filename_pair
            filepath = save_dir / filename
            if filepath in retrieved_files_set:
                attr = AlphaFoldEntry.format2attr(af_format)
                setattr(entry, attr, filepath)
            # else: File was not found (404) during download, so we leave the attribute as None

        yield entry


def fetch_many_async(
    uniprot_accessions: Iterable[str],
    save_dir: Path,
    formats: set[DownloadableFormat],
    db_version: str | None = None,
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
    gzip_files: bool = False,
    all_isoforms: bool = False,
) -> AsyncGenerator[AlphaFoldEntry]:
    """Asynchronously fetches summaries and/or files from
    [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).

    Args:
        uniprot_accessions: A set of Uniprot accessions to fetch.
        save_dir: The directory to save the fetched files to.
        formats: A set of formats to download.
            If `summary` is in the set then summaries will be fetched using the API endpoint.
            and later the other files will be downloaded using static file URLs.
            If `summary` is not in the set then all files will be downloaded using static file
            URLs only.
        db_version: The version of the AlphaFold database to use. If None, the latest version will be used.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: A cacher to use for caching the fetched files.
        gzip_files: Whether to gzip the downloaded files.
            Summaries are never gzipped.
        all_isoforms: Whether to yield all isoforms of each uniprot entry.
            When False then yields only the canonical sequence per uniprot entry.

    Yields:
        A dataclass containing the summary, pdb file, and pae file.

    Raises:
        ValueError: If 'formats' set is empty.
        ValueError: If all_isoforms is True and 'summary' is not in 'formats' set.
    """
    if len(formats) == 0:
        msg = "At least one format must be specified. The 'formats' argument is empty."
        raise ValueError(msg)
    if "summary" in formats:
        if db_version is not None:
            logger.warning("db_version is ignored when 'summary' is in 'formats' set. Always uses latest version.")
        return _fetch_many_async_with_summary(
            uniprot_accessions,
            save_dir,
            formats,
            max_parallel_downloads=max_parallel_downloads,
            cacher=cacher,
            gzip_files=gzip_files,
            all_isoforms=all_isoforms,
        )
    if all_isoforms:
        msg = "Cannot fetch all isoforms when 'summary' is not in 'formats' set."
        raise ValueError(msg)
    return _fetch_many_async_without_summary(
        uniprot_accessions,
        save_dir,
        formats,
        db_version=db_version,
        max_parallel_downloads=max_parallel_downloads,
        cacher=cacher,
        gzip_files=gzip_files,
    )


# jscpd:ignore-start  # noqa: ERA001
def fetch_many(
    uniprot_accessions: Iterable[str],
    save_dir: Path,
    formats: set[DownloadableFormat],
    db_version: str | None = None,
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
    gzip_files: bool = False,
    all_isoforms: bool = False,
    # jscpd:ignore-end  # noqa: ERA001
) -> list[AlphaFoldEntry]:
    """Synchronously fetches summaries and/or files like cif from AlphaFold Protein Structure Database.

    Args:
        uniprot_accessions: A set of Uniprot accessions to fetch.
        save_dir: The directory to save the fetched files to.
        formats: A set of formats to download.
            If `summary` is in the set then summaries will be fetched using the API endpoint.
            and later the other files will be downloaded using static file URLs.
            If `summary` is not in the set then all files will be downloaded using static file
            URLs only.
            Excluding 'summary' is much faster as it avoids slow API calls.
        db_version: The version of the AlphaFold database to use. If None, the latest version will be used.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: A cacher to use for caching the fetched files.
        gzip_files: Whether to gzip the downloaded files.
            Summaries are never gzipped.
        all_isoforms: Whether to yield all isoforms of each uniprot entry.
            When False then yields only the canonical sequence per uniprot entry.

    Returns:
        A list of AlphaFoldEntry dataclasses containing the summary, pdb file, and pae file.
    """

    async def gather_entries():
        return [
            entry
            async for entry in fetch_many_async(
                uniprot_accessions,
                save_dir,
                formats,
                db_version=db_version,
                max_parallel_downloads=max_parallel_downloads,
                cacher=cacher,
                gzip_files=gzip_files,
                all_isoforms=all_isoforms,
            )
        ]

    return run_async(gather_entries())
