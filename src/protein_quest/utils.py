"""Module for functions that are used in multiple places."""

import asyncio
import hashlib
import logging
import shutil
from collections.abc import Coroutine, Iterable
from contextlib import asynccontextmanager
from functools import lru_cache
from pathlib import Path
from textwrap import dedent
from typing import Any, Literal, Protocol, get_args, runtime_checkable

import aiofiles
import aiofiles.os
import aiohttp
from aiohttp.streams import AsyncStreamIterator
from aiohttp_retry import ExponentialRetry, RetryClient
from platformdirs import user_cache_dir
from tqdm.asyncio import tqdm
from yarl import URL

logger = logging.getLogger(__name__)

CopyMethod = Literal["copy", "symlink", "hardlink"]
"""Methods for copying files."""
copy_methods = set(get_args(CopyMethod))
"""Set of valid copy methods."""


@lru_cache
def _cache_sub_dir(root_cache_dir: Path, filename: str) -> Path:
    """Get the cache sub-directory for a given path.

    To not have too many files in a single directory,
    we create sub-directories based on the hash of the filename.

    Args:
        root_cache_dir: The root directory for the cache.
        filename: The filename to be cached.
    Returns:
        The parent path to the cached file.

    """
    cache_sub_dir = hashlib.sha256(filename.encode("utf-8")).hexdigest()[:4]
    cache_sub_dir_path = root_cache_dir / cache_sub_dir
    cache_sub_dir_path.mkdir(parents=True, exist_ok=True)
    return cache_sub_dir_path


@runtime_checkable
class Cacher(Protocol):
    """Protocol for a cacher."""

    def __contains__(self, item: str | Path) -> bool:
        """Check if a file is in the cache.

        Args:
            item: The filename or Path to check.

        Returns:
            True if the file is in the cache, False otherwise.
        """
        ...

    async def copy_from_cache(self, target: Path) -> Path | None:
        """Copy a file from the cache to a target location if it exists in the cache.

        Assumes:

        - target does not exist.
        - the parent directory of target exists.

        Args:
            target: The path to copy the file to.

        Returns:
            The path to the cached file if it was copied, None otherwise.
        """
        ...

    async def write_iter(self, target: Path, content: AsyncStreamIterator[bytes]) -> Path:
        """Write content to a file and cache it.

        Args:
            target: The path to write the content to.
            content: An async iterator that yields bytes to write to the file.

        Returns:
            The path to the cached file.
        """
        ...

    async def write_bytes(self, target: Path, content: bytes) -> Path:
        """Write bytes to a file and cache it.

        Args:
            target: The path to write the content to.
            content: The bytes to write to the file.

        Returns:
            The path to the cached file.
        """
        ...


class NoopCacher(Cacher):
    """A cacher that caches nothing.

    On writes it just writes to the target path.
    """

    def __contains__(self, item: str | Path) -> bool:
        # We don't have anything cached ever
        return False

    async def copy_from_cache(self, target: Path) -> Path | None:  # noqa: ARG002
        # We don't have anything cached ever
        return None

    async def write_iter(self, target: Path, content: AsyncStreamIterator[bytes]) -> Path:
        target.write_bytes(b"".join([chunk async for chunk in content]))
        return target

    async def write_bytes(self, target: Path, content: bytes) -> Path:
        target.write_bytes(content)
        return target


def user_cache_root_dir() -> Path:
    """Get the users root directory for caching files.

    Returns:
        The path to the user's cache directory for protein-quest.
    """
    return Path(user_cache_dir("protein-quest"))


class DirectoryCacher(Cacher):
    """Class to cache files in a directory."""

    def __init__(
        self,
        cache_dir: Path | None = None,
        copy_method: CopyMethod = "hardlink",
    ) -> None:
        """Initialize the cacher.

        Args:
            cache_dir: The directory to use for caching.
                If None, a default cache directory (~/.cache/protein-quest) is used.
            copy_method: The method to use for copying.
        """
        if cache_dir is None:
            cache_dir = user_cache_root_dir()
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if copy_method not in copy_methods:
            msg = f"Unknown copy method: {copy_method}. Must be one of {copy_methods}."
            raise ValueError(msg)
        self.copy_method: CopyMethod = copy_method

    def __contains__(self, item: str | Path) -> bool:
        cached_file = self._location(item)
        return cached_file.exists()

    def _location(self, item: str | Path) -> Path:
        """Get the location of a cached file.

        Args:
            item: The filename or Path to get the location for.

        Returns:
            The path to the cached file.
        """
        file_name = item.name if isinstance(item, Path) else item
        cache_sub_dir = _cache_sub_dir(self.cache_dir, file_name)
        return cache_sub_dir / file_name

    async def copy_from_cache(self, target: Path) -> Path | None:
        cached_file = self._location(target.name)
        exists = await aiofiles.os.path.exists(str(cached_file))
        if exists:
            await async_copyfile(cached_file, target, copy_method=self.copy_method)
            return cached_file
        return None

    async def write_iter(self, target: Path, content: AsyncStreamIterator[bytes]) -> Path:
        cached_file = self._location(target.name)
        async with aiofiles.open(cached_file, "xb") as f:
            async for chunk in content:
                await f.write(chunk)
        await async_copyfile(cached_file, target, copy_method=self.copy_method)
        return cached_file

    async def write_bytes(self, target: Path, content: bytes) -> Path:
        cached_file = self._location(target.name)
        async with aiofiles.open(cached_file, "xb") as f:
            await f.write(content)
        await async_copyfile(cached_file, target, copy_method=self.copy_method)
        return cached_file


async def retrieve_files(
    urls: Iterable[tuple[URL | str, str]],
    save_dir: Path,
    max_parallel_downloads: int = 5,
    retries: int = 3,
    total_timeout: int = 300,
    desc: str = "Downloading files",
    cacher: Cacher | None = None,
    chunk_size: int = 524288,  # 512 KiB
) -> list[Path]:
    """Retrieve files from a list of URLs and save them to a directory.

    Args:
        urls: A list of tuples, where each tuple contains a URL and a filename.
        save_dir: The directory to save the downloaded files to.
        max_parallel_downloads: The maximum number of files to download in parallel.
        retries: The number of times to retry a failed download.
        total_timeout: The total timeout for a download in seconds.
        desc: Description for the progress bar.
        cacher: An optional cacher to use for caching files.
        chunk_size: The size of each chunk to read from the response.

    Returns:
        A list of paths to the downloaded files.
    """
    save_dir.mkdir(parents=True, exist_ok=True)
    semaphore = asyncio.Semaphore(max_parallel_downloads)
    async with friendly_session(retries, total_timeout) as session:
        tasks = [
            _retrieve_file(
                session=session,
                url=url,
                save_path=save_dir / filename,
                semaphore=semaphore,
                cacher=cacher,
                chunk_size=chunk_size,
            )
            for url, filename in urls
        ]
        files: list[Path] = await tqdm.gather(*tasks, desc=desc)
        return files


async def _retrieve_file(
    session: RetryClient,
    url: URL | str,
    save_path: Path,
    semaphore: asyncio.Semaphore,
    cacher: Cacher | None = None,
    chunk_size: int = 524288,  # 512 KiB
) -> Path:
    """Retrieve a single file from a URL and save it to a specified path.

    Args:
        session: The aiohttp session to use for the request.
        url: The URL to download the file from.
        save_path: The path where the file should be saved.
        semaphore: A semaphore to limit the number of concurrent downloads.
        cacher: An optional cacher to use for caching files.
        chunk_size: The size of each chunk to read from the response.

    Returns:
        The path to the saved file.
    """
    if save_path.exists():
        logger.debug(f"File {save_path} already exists. Skipping download from {url}.")
        return save_path

    if cacher is None:
        cacher = NoopCacher()
    if cached_file := await cacher.copy_from_cache(save_path):
        logger.debug(f"File {save_path} was copied from cache {cached_file}. Skipping download from {url}.")
        return save_path

    async with (
        semaphore,
        session.get(url) as resp,
    ):
        resp.raise_for_status()
        await cacher.write_iter(save_path, resp.content.iter_chunked(chunk_size))
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


class NestedAsyncIOLoopError(RuntimeError):
    """Custom error for nested async I/O loops."""

    def __init__(self) -> None:
        msg = dedent("""\
            Can not run async method from an environment where the asyncio event loop is already running.
            Like a Jupyter notebook.

            Please use the async function directly or
            call `import nest_asyncio; nest_asyncio.apply()` and try again.
            """)
        super().__init__(msg)


def run_async[R](coroutine: Coroutine[Any, Any, R]) -> R:
    """Run an async coroutine with nicer error.

    Args:
        coroutine: The async coroutine to run.

    Returns:
        The result of the coroutine.

    Raises:
        NestedAsyncIOLoopError: If called from a nested async I/O loop like in a Jupyter notebook.
    """
    try:
        return asyncio.run(coroutine)
    except RuntimeError as e:
        raise NestedAsyncIOLoopError from e


def copyfile(source: Path, target: Path, copy_method: CopyMethod = "copy"):
    """Make target path be same file as source by either copying or symlinking or hardlinking.

    Note that hardlink copy method only work within the same filesystem and are harder to track.
    If you want to track cached files easily then use 'symlink'.
    On Windows you need developer mode or admin privileges to create symlinks.

    Args:
        source: The source file to copy or link.
        target: The target file to create.
        copy_method: The method to use for copying.

    Raises:
        FileNotFoundError: If the source file or parent of target does not exist.
        ValueError: If an unknown copy method is provided.
    """
    rel_source = source.relative_to(target.parent, walk_up=True)
    if copy_method == "copy":
        shutil.copyfile(source, target)
    elif copy_method == "symlink":
        target.symlink_to(rel_source)
    elif copy_method == "hardlink":
        target.hardlink_to(source)
    else:
        msg = f"Unknown method: {copy_method}"
        raise ValueError(msg)


async def async_copyfile(
    source: Path,
    target: Path,
    copy_method: CopyMethod = "copy",
):
    """Asynchronously copy a file from source to target using aiofiles.

    Note that hardlink copy method only work within the same filesystem and are harder to track.
    If you want to track cached files easily then use 'symlink'.
    On Windows you need developer mode or admin privileges to create symlinks.

    Args:
        source: The source file to copy.
        target: The target file to create.
        copy_method: The method to use for copying.

    Raises:
        FileNotFoundError: If the source file or parent of target does not exist.
        ValueError: If an unknown copy method is provided.
    """
    if copy_method == "copy":
        # Could use loop of chunks with aiofiles,
        # but shutil is ~1.9x faster on my machine
        # due to fastcopy and sendfile optimizations in shutil.
        await asyncio.to_thread(shutil.copyfile, source, target)
    elif copy_method == "symlink":
        rel_source = source.relative_to(target.parent, walk_up=True)
        await aiofiles.os.symlink(str(rel_source), str(target))
    elif copy_method == "hardlink":
        await aiofiles.os.link(str(source), str(target))
    else:
        msg = f"Unknown method: {copy_method}"
        raise ValueError(msg)
