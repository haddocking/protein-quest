import asyncio
from collections.abc import Iterable
from contextlib import asynccontextmanager
from pathlib import Path

import aiofiles
import aiohttp
from aiohttp_retry import ExponentialRetry, RetryClient
from tqdm.asyncio import tqdm


async def retrieve_files(
    urls: Iterable[tuple[str, str]],
    save_dir: Path,
    max_parallel_downloads: int = 5,
    retries: int = 3,
    total_timeout: int = 300,
    desc: str = "Downloading files",
) -> list[Path]:
    """Retrieve files from a list of URLs and save them to a directory.

    Args:
        urls: A list of tuples, where each tuple contains a URL and a filename.
        save_dir: The directory to save the downloaded files to.
        max_parallel_downloads: The maximum number of files to download in parallel.
        retries: The number of times to retry a failed download.
        total_timeout: The total timeout for a download in seconds.
        desc: Description for the progress bar.

    Returns:
        A list of paths to the downloaded files.
    """
    save_dir.mkdir(parents=True, exist_ok=True)
    semaphore = asyncio.Semaphore(max_parallel_downloads)
    async with friendly_session(retries, total_timeout) as session:
        tasks = [_retrieve_file(session, url, save_dir / filename, semaphore) for url, filename in urls]
        files: list[Path] = await tqdm.gather(*tasks, desc=desc)
        return files


async def _retrieve_file(
    session: RetryClient,
    url: str,
    save_path: Path,
    semaphore: asyncio.Semaphore,
    ovewrite: bool = False,
    chunk_size: int = 131072,  # 128 KiB
) -> Path:
    """Retrieve a single file from a URL and save it to a specified path.

    Args:
        session: The aiohttp session to use for the request.
        url: The URL to download the file from.
        save_path: The path where the file should be saved.
        semaphore: A semaphore to limit the number of concurrent downloads.
        ovewrite: Whether to overwrite the file if it already exists.
        chunk_size: The size of each chunk to read from the response.

    Returns:
        The path to the saved file.
    """
    if save_path.exists():
        if ovewrite:
            save_path.unlink()
        else:
            return save_path
    async with (
        semaphore,
        aiofiles.open(save_path, "xb") as f,
        session.get(url) as resp,
    ):
        resp.raise_for_status()
        async for chunk in resp.content.iter_chunked(chunk_size):
            await f.write(chunk)
    return save_path


@asynccontextmanager
async def friendly_session(retries: int = 3, total_timeout: int = 300):
    """Create an aiohttp session with retry capabilities.

    Examples:
        Use as async context:

        >>> async with friendly_session(retries=5, total_timeout=60) as session:
        >>>     r = await session.get("https://example.com/api/data")
        >>>     print(r)
        <ClientResponse(https://example.com/api/data) [404 Not Found]>
        <CIMultiDictProxy('Accept-Ranges': 'bytes', ...

    Args:
        retries: The number of retry attempts for failed requests.
        total_timeout: The total timeout for a request in seconds.
    """
    retry_options = ExponentialRetry(attempts=retries)
    timeout = aiohttp.ClientTimeout(total=total_timeout)  # pyrefly: ignore false positive
    async with aiohttp.ClientSession(timeout=timeout) as session:
        client = RetryClient(client_session=session, retry_options=retry_options)
        yield client
