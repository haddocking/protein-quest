"""Module for fetching structures from PDBe."""

from collections.abc import Iterable, Mapping
from pathlib import Path

from protein_quest.utils import Cacher, retrieve_files, run_async


def _map_id_mmcif(pdb_id: str) -> tuple[str, str]:
    """
    Map PDB id to a download gzipped mmCIF url and file.

    For example for PDB id "8WAS", the url will be
    "https://www.ebi.ac.uk/pdbe/entry-files/download/8was.cif.gz" and the file will be "8was.cif.gz".

    Args:
        pdb_id: The PDB ID to map.

    Returns:
        A tuple containing the URL to download the mmCIF file and the filename.
    """
    fn = f"{pdb_id.lower()}.cif.gz"
    # On PDBe you can sometimes download an updated mmCIF file,
    # Current url is for the archive mmCIF file
    # TODO check if archive is OK, or if we should try to download the updated file
    # this will cause many more requests, so we should only do this if needed
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{fn}"
    return url, fn


async def fetch(
    ids: Iterable[str], save_dir: Path, max_parallel_downloads: int = 5, cacher: Cacher | None = None
) -> Mapping[str, Path]:
    """Fetches mmCIF files from the PDBe database.

    Args:
        ids: A set of PDB IDs to fetch.
        save_dir: The directory to save the fetched mmCIF files to.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: An optional cacher to use for caching downloaded files.

    Returns:
        A dict of id and paths to the downloaded mmCIF files.
    """

    # The future result, is in a different order than the input ids,
    # so we need to map the ids to the urls and filenames.

    id2urls = {pdb_id: _map_id_mmcif(pdb_id) for pdb_id in ids}
    urls = list(id2urls.values())
    id2paths = {pdb_id: save_dir / fn for pdb_id, (_, fn) in id2urls.items()}

    await retrieve_files(urls, save_dir, max_parallel_downloads, desc="Downloading PDBe mmCIF files", cacher=cacher)
    return id2paths


def sync_fetch(ids: Iterable[str], save_dir: Path, max_parallel_downloads: int = 5) -> Mapping[str, Path]:
    """Synchronously fetches mmCIF files from the PDBe database.

    Args:
        ids: A set of PDB IDs to fetch.
        save_dir: The directory to save the fetched mmCIF files to.
        max_parallel_downloads: The maximum number of parallel downloads.

    Returns:
        A dict of id and paths to the downloaded mmCIF files.
    """
    return run_async(fetch(ids, save_dir, max_parallel_downloads))
