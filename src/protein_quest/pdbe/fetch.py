"""Module for fetching structures from PDBe."""

from collections.abc import Iterable, Mapping
from pathlib import Path

from protein_quest.utils import Cacher, read_ids_from_csv, retrieve_files, run_async


def _map_id_mmcif(pdb_id: str, archived: bool = False) -> tuple[str, str]:
    """
    Map PDB id to a download gzipped mmCIF url and file.

    For example for PDB id "8WAS", the url will be
    "https://www.ebi.ac.uk/pdbe/entry-files/download/8was_updated.cif.gz" and the file will be "8was_updated.cif.gz".
    If archived is True, the url will be
    "https://www.ebi.ac.uk/pdbe/entry-files/download/8was.cif.gz" and the file will be "8was.cif.gz".

    Args:
        pdb_id: The PDB ID to map.
        archived: Whether to fetch archived versions of the files or default updated version.

    Returns:
        A tuple containing the URL to download the mmCIF file and the filename.
    """
    fn = f"{pdb_id.lower()}_updated.cif.gz"
    if archived:
        fn = f"{pdb_id.lower()}.cif.gz"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{fn}"
    return url, fn


async def fetch(
    ids: Iterable[str],
    save_dir: Path,
    archived: bool = False,
    max_parallel_downloads: int = 5,
    cacher: Cacher | None = None,
) -> Mapping[str, Path]:
    """Fetches mmCIF files from the PDBe database.

    Args:
        ids: A set of PDB IDs to fetch.
        save_dir: The directory to save the fetched mmCIF files to.
        archived: Whether to fetch archived versions of the files or default updated version.
            The updated version is generated with standardisation of vocabularies,
            and addition of connectivity information for every chemical compound present in the PDB entry.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: An optional cacher to use for caching downloaded files.

    Returns:
        A dict of id and paths to the downloaded mmCIF files.
    """

    # The future result, is in a different order than the input ids,
    # so we need to map the ids to the urls and filenames.

    id2urls = {pdb_id: _map_id_mmcif(pdb_id, archived=archived) for pdb_id in ids}
    urls = list(id2urls.values())
    id2paths = {pdb_id: save_dir / fn for pdb_id, (_, fn) in id2urls.items()}

    await retrieve_files(urls, save_dir, max_parallel_downloads, desc="Downloading PDBe mmCIF files", cacher=cacher)
    return id2paths


def sync_fetch(
    ids: Iterable[str], save_dir: Path, archived: bool = False, max_parallel_downloads: int = 5
) -> Mapping[str, Path]:
    """Synchronously fetches mmCIF files from the PDBe database.

    Args:
        ids: A set of PDB IDs to fetch.
        save_dir: The directory to save the fetched mmCIF files to.
        archived: Whether to fetch archived versions of the files or default updated version.
            The updated version is generated with standardisation of vocabularies,
            and addition of connectivity information for every chemical compound present in the PDB entry.
        max_parallel_downloads: The maximum number of parallel downloads.

    Returns:
        A dict of id and paths to the downloaded mmCIF files.
    """
    return run_async(fetch(ids, save_dir, archived=archived, max_parallel_downloads=max_parallel_downloads))


def read_pdb_ids_from_csv(file: Path) -> set[str]:
    """Reads PDB IDs from a CSV file.

    The CSV file can provide PDB IDs in the ``pdb_id`` column. It can also
    provide generic identifiers through the ``model_provider`` and
    ``model_identifier`` columns. In that case, only rows with
    ``model_provider == "pdbe"`` are used. If the CSV contains only one
    column, every value in that column is treated as an ID, including the
    first row.

    Arguments:
        file: A path to a file containing the CSV data.

    Returns:
        A set of PDB IDs extracted from the CSV file.
    """
    return read_ids_from_csv(file, id_column="pdb_id", model_provider="pdbe")
