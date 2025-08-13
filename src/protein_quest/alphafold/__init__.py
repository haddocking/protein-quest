import asyncio
import logging
from asyncio import Semaphore
from collections.abc import AsyncGenerator, Iterable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from aiohttp_retry import RetryClient
from cattrs import structure
from tqdm.asyncio import tqdm

from protein_quest.alphafold.entry_summary import EntrySummary
from protein_quest.utils import friendly_session, retrieve_files

logger = logging.getLogger(__name__)


DownloadableFormat = Literal[
    "bcif",
    "cif",
    "pdb",
    "paeImage",
    "paeDoc",
    "amAnnotations",
    "amAnnotationsHg19",
    "amAnnotationsHg38",
]
"""Types of formats that can be downloaded from the AlphaFold web service."""

downloadable_formats: set[DownloadableFormat] = {
    "bcif",
    "cif",
    "pdb",
    "paeImage",
    "paeDoc",
    "amAnnotations",
    "amAnnotationsHg19",
    "amAnnotationsHg38",
}
"""Set of formats that can be downloaded from the AlphaFold web service."""


def _camel_to_snake_case(name: str) -> str:
    """Convert a camelCase string to snake_case."""
    return "".join(["_" + c.lower() if c.isupper() else c for c in name]).lstrip("_")


@dataclass
class AlphaFoldEntry:
    """AlphaFoldEntry represents a minimal single entry in the AlphaFold database.

    See https://alphafold.ebi.ac.uk/api-docs for more details on the API and data structure.
    """

    uniprot_acc: str
    summary: EntrySummary | None
    bcif_file: Path | None = None
    cif_file: Path | None = None
    pdb_file: Path | None = None
    pae_image_file: Path | None = None
    pae_doc_file: Path | None = None
    am_annotations_file: Path | None = None
    am_annotations_hg19_file: Path | None = None
    am_annotations_hg38_file: Path | None = None

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


async def fetch_summmary(qualifier: str, session: RetryClient, semaphore: Semaphore) -> list[EntrySummary]:
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{qualifier}"
    async with semaphore, session.get(url) as response:
        response.raise_for_status()
        data = await response.json()
        return structure(data, list[EntrySummary])


async def fetch_summaries(qualifiers: Iterable[str], max_parallel_downloads: int = 5) -> AsyncGenerator[EntrySummary]:
    semaphore = Semaphore(max_parallel_downloads)
    async with friendly_session() as session:
        tasks = [fetch_summmary(qualifier, session, semaphore) for qualifier in qualifiers]
        summaries_per_qualifier: list[list[EntrySummary]] = await tqdm.gather(
            *tasks, desc="Fetching Alphafold summaries"
        )
        for summaries in summaries_per_qualifier:
            for summary in summaries:
                yield summary


def url2name(url: str) -> str:
    """Given a URL, return the final path component as the name of the file."""
    return url.split("/")[-1]


async def fetch_many_async(
    ids: Iterable[str], save_dir: Path, what: set[DownloadableFormat], max_parallel_downloads: int = 5
) -> AsyncGenerator[AlphaFoldEntry]:
    """Asynchronously fetches summaries and pdb and pae (predicted alignment error) files from
    [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).

    Args:
        ids: A set of Uniprot IDs to fetch.
        save_dir: The directory to save the fetched files to.
        what: A set of formats to download.
        max_parallel_downloads: The maximum number of parallel downloads.

    Yields:
        A dataclass containing the summary, pdb file, and pae file.
    """
    summaries = [s async for s in fetch_summaries(ids, max_parallel_downloads=max_parallel_downloads)]

    files = files_to_download(what, summaries)

    await retrieve_files(
        files,
        save_dir,
        desc="Downloading AlphaFold files",
        max_parallel_downloads=max_parallel_downloads,
    )
    for summary in summaries:
        yield AlphaFoldEntry(
            uniprot_acc=summary.uniprotAccession,
            summary=summary,
            bcif_file=save_dir / url2name(summary.bcifUrl) if "bcif" in what else None,
            cif_file=save_dir / url2name(summary.cifUrl) if "cif" in what else None,
            pdb_file=save_dir / url2name(summary.pdbUrl) if "pdb" in what else None,
            pae_image_file=save_dir / url2name(summary.paeImageUrl) if "paeImage" in what else None,
            pae_doc_file=save_dir / url2name(summary.paeDocUrl) if "paeDoc" in what else None,
            am_annotations_file=(
                save_dir / url2name(summary.amAnnotationsUrl)
                if "amAnnotations" in what and summary.amAnnotationsUrl
                else None
            ),
            am_annotations_hg19_file=(
                save_dir / url2name(summary.amAnnotationsHg19Url)
                if "amAnnotationsHg19" in what and summary.amAnnotationsHg19Url
                else None
            ),
            am_annotations_hg38_file=(
                save_dir / url2name(summary.amAnnotationsHg38Url)
                if "amAnnotationsHg38" in what and summary.amAnnotationsHg38Url
                else None
            ),
        )


def files_to_download(what: set[DownloadableFormat], summaries: Iterable[EntrySummary]) -> set[tuple[str, str]]:
    if not (set(what) <= downloadable_formats):
        msg = (
            f"Invalid format(s) specified: {set(what) - downloadable_formats}. "
            f"Valid formats are: {downloadable_formats}"
        )
        raise ValueError(msg)

    files: set[tuple[str, str]] = set()
    for summary in summaries:
        for fmt in what:
            url = getattr(summary, f"{fmt}Url", None)
            if url is None:
                logger.warning(f"Summary {summary.uniprotAccession} does not have a URL for format '{fmt}'. Skipping.")
                continue
            file = (url, url2name(url))
            files.add(file)
    return files


def fetch_many(
    ids: Iterable[str], save_dir: Path, what: set[DownloadableFormat], max_parallel_downloads: int = 5
) -> list[AlphaFoldEntry]:
    """Synchronously fetches summaries and pdb and pae files from AlphaFold Protein Structure Database.

    Args:
        ids: A set of Uniprot IDs to fetch.
        save_dir: The directory to save the fetched files to.
        what: A set of formats to download.
        max_parallel_downloads: The maximum number of parallel downloads.

    Returns:
        A list of AlphaFoldEntry dataclasses containing the summary, pdb file, and pae file.
    """

    async def gather_entries():
        return [
            entry
            async for entry in fetch_many_async(ids, save_dir, what, max_parallel_downloads=max_parallel_downloads)
        ]

    def run_async_task():
        return asyncio.run(gather_entries())

    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(run_async_task)
        return future.result()


def relative_to(entry: AlphaFoldEntry, session_dir: Path) -> AlphaFoldEntry:
    """Convert paths in an AlphaFoldEntry to be relative to the session directory.

    Args:
        entry: An AlphaFoldEntry instance with absolute paths.
        session_dir: The session directory to which the paths should be made relative.

    Returns:
        An AlphaFoldEntry instance with paths relative to the session directory.
    """
    return AlphaFoldEntry(
        uniprot_acc=entry.uniprot_acc,
        summary=entry.summary,
        bcif_file=entry.bcif_file.relative_to(session_dir) if entry.bcif_file else None,
        cif_file=entry.cif_file.relative_to(session_dir) if entry.cif_file else None,
        pdb_file=entry.pdb_file.relative_to(session_dir) if entry.pdb_file else None,
        pae_image_file=entry.pae_image_file.relative_to(session_dir) if entry.pae_image_file else None,
        pae_doc_file=entry.pae_doc_file.relative_to(session_dir) if entry.pae_doc_file else None,
        am_annotations_file=entry.am_annotations_file.relative_to(session_dir) if entry.am_annotations_file else None,
        am_annotations_hg19_file=(
            entry.am_annotations_hg19_file.relative_to(session_dir) if entry.am_annotations_hg19_file else None
        ),
        am_annotations_hg38_file=(
            entry.am_annotations_hg38_file.relative_to(session_dir) if entry.am_annotations_hg38_file else None
        ),
    )
