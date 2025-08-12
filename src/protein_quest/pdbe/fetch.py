import asyncio
import concurrent.futures
from collections.abc import Iterable, Mapping
from pathlib import Path

from protein_quest.utils import retrieve_files


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


def fetch(ids: Iterable[str], save_dir: Path, max_parallel_downloads: int = 5) -> Mapping[str, Path]:
    """Fetches mmCIF files from the PDBe database.

    Args:
        ids: A set of PDB IDs to fetch.
        save_dir: The directory to save the fetched mmCIF files to.
        max_parallel_downloads: The maximum number of parallel downloads.

    Returns:
        A dict of id and paths to the downloaded mmCIF files.
    """

    # The future result, is in a different order than the input ids,
    # so we need to map the ids to the urls and filenames.

    id2urls = {pdb_id: _map_id_mmcif(pdb_id) for pdb_id in ids}
    urls = list(id2urls.values())
    id2paths = {pdb_id: save_dir / fn for pdb_id, (_, fn) in id2urls.items()}

    # TODO switch to nest_asyncio instead of this ThreadPoolExecutor
    def run_async_task():
        return asyncio.run(retrieve_files(urls, save_dir, max_parallel_downloads, desc="Downloading PDBe mmCIF files"))

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(run_async_task)
        result = future.result()
        if set(result) != set(id2paths.values()):
            msg = "Not all files were downloaded successfully."
            raise ValueError(msg)
        return id2paths
