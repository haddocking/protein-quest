"""Module dealing with Electron Microscopy Data Bank (EMDB)."""

from collections.abc import Iterable, Mapping
from pathlib import Path

from protein_quest.utils import Cacher, retrieve_files


def _map_id2volume_url(emdb_id: str) -> tuple[str, str]:
    # https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-19583/map/emd_19583.map.gz
    fn = emdb_id.lower().replace("emd-", "emd_") + ".map.gz"
    url = f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emdb_id}/map/{fn}"
    return url, fn


async def fetch(
    emdb_ids: Iterable[str], save_dir: Path, max_parallel_downloads: int = 1, cacher: Cacher | None = None
) -> Mapping[str, Path]:
    """Fetches volume files from the EMDB database.

    Args:
        emdb_ids: A list of EMDB IDs to fetch.
        save_dir: The directory to save the downloaded files.
        max_parallel_downloads: The maximum number of parallel downloads.
        cacher: An optional cacher to use for caching downloaded files.

    Returns:
        A mapping of EMDB IDs to their downloaded files.
    """
    id2urls = {emdb_id: _map_id2volume_url(emdb_id) for emdb_id in emdb_ids}
    urls = list(id2urls.values())
    id2paths = {emdb_id: save_dir / fn for emdb_id, (_, fn) in id2urls.items()}

    # TODO show progress of each item
    # TODO handle failed downloads, by skipping them instead of raising an error
    await retrieve_files(urls, save_dir, max_parallel_downloads, desc="Downloading EMDB volume files", cacher=cacher)
    return id2paths
