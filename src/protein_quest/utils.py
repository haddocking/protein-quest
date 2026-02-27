"""Module for functions that are used in multiple places."""

import argparse
import asyncio
import hashlib
import logging
import shutil
from collections.abc import Coroutine, Iterable, Sequence
from contextlib import asynccontextmanager
from functools import lru_cache
from pathlib import Path
from textwrap import dedent
from typing import Any, Literal, Protocol, get_args, runtime_checkable

import aiofiles
import aiofiles.os
import aiohttp
import rich
from aiohttp.streams import AsyncStreamIterator
from aiohttp_retry import ExponentialRetry, RetryClient
from platformdirs import user_cache_dir
from rich_argparse import ArgumentDefaultsRichHelpFormatter
from tqdm.asyncio import tqdm
from yarl import URL

logger = logging.getLogger(__name__)

CopyMethod = Literal["copy", "symlink", "hardlink"]
"""Methods for copying files."""
copy_methods = set(get_args(CopyMethod))
"""Set of valid copy methods."""


@lru_cache
def _cache_sub_dir(root_cache_dir: Path, filename: str, hash_length: int = 4) -> Path:
    """Get the cache sub-directory for a given path.

    To not have too many files in a single directory,
    we create sub-directories based on the hash of the filename.

    Args:
        root_cache_dir: The root directory for the cache.
        filename: The filename to be cached.
        hash_length: The length of the hash to use for the sub-directory.

    Returns:
        The parent path to the cached file.
    """
    full_hash = hashlib.blake2b(filename.encode("utf-8")).hexdigest()
    cache_sub_dir = full_hash[:hash_length]
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

        Raises:
            FileExistsError: If the target file already exists.
        """
        ...

    async def write_bytes(self, target: Path, content: bytes) -> Path:
        """Write bytes to a file and cache it.

        Args:
            target: The path to write the content to.
            content: The bytes to write to the file.

        Returns:
            The path to the cached file.

        Raises:
            FileExistsError: If the target file already exists.
        """
        ...


class PassthroughCacher(Cacher):
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
        if target.exists():
            raise FileExistsError(target)
        target.write_bytes(b"".join([chunk async for chunk in content]))
        return target

    async def write_bytes(self, target: Path, content: bytes) -> Path:
        if target.exists():
            raise FileExistsError(target)
        target.write_bytes(content)
        return target


def user_cache_root_dir() -> Path:
    """Get the users root directory for caching files.

    Returns:
        The path to the user's cache directory for protein-quest.
    """
    return Path(user_cache_dir("protein-quest"))


class DirectoryCacher(Cacher):
    """Class to cache files in a directory.

    Caching logic is based on the file name only.
    If file name of paths are the same then the files are considered the same.

    Attributes:
        cache_dir: The directory to use for caching.
        copy_method: The method to use for copying files.
    """

    def __init__(
        self,
        cache_dir: Path | None = None,
        copy_method: CopyMethod = "hardlink",
    ) -> None:
        """Initialize the cacher.

        If file name of paths are the same then the files are considered the same.

        Args:
            cache_dir: The directory to use for caching.
                If None, a default cache directory (~/.cache/protein-quest) is used.
            copy_method: The method to use for copying.
        """
        if cache_dir is None:
            cache_dir = user_cache_root_dir()
        self.cache_dir: Path = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if copy_method == "copy":
            logger.warning(
                "Using copy as copy_method to cache files is not recommended. "
                "This will use more disk space and be slower than symlink or hardlink."
            )
        if copy_method not in copy_methods:
            msg = f"Unknown copy method: {copy_method}. Must be one of {copy_methods}."
            raise ValueError(msg)
        self.copy_method: CopyMethod = copy_method

    def __contains__(self, item: str | Path) -> bool:
        cached_file = self._as_cached_path(item)
        return cached_file.exists()

    def _as_cached_path(self, item: str | Path) -> Path:
        file_name = item.name if isinstance(item, Path) else item
        cache_sub_dir = _cache_sub_dir(self.cache_dir, file_name)
        return cache_sub_dir / file_name

    async def copy_from_cache(self, target: Path) -> Path | None:
        cached_file = self._as_cached_path(target.name)
        exists = await aiofiles.os.path.exists(str(cached_file))
        if exists:
            await async_copyfile(cached_file, target, copy_method=self.copy_method)
            return cached_file
        return None

    async def write_iter(self, target: Path, content: AsyncStreamIterator[bytes]) -> Path:
        cached_file = self._as_cached_path(target.name)
        # Write file to cache dir
        async with aiofiles.open(cached_file, "xb") as f:
            async for chunk in content:
                await f.write(chunk)
        # Copy to target location
        await async_copyfile(cached_file, target, copy_method=self.copy_method)
        return cached_file

    async def write_bytes(self, target: Path, content: bytes) -> Path:
        cached_file = self._as_cached_path(target.name)
        # Write file to cache dir
        async with aiofiles.open(cached_file, "xb") as f:
            await f.write(content)
        # Copy to target location
        await async_copyfile(cached_file, target, copy_method=self.copy_method)
        return cached_file

    def populate_cache(self, source_dir: Path) -> dict[Path, Path]:
        """Populate the cache from an existing directory.

        This will copy all files from the source directory to the cache directory.
        If a file with the same name already exists in the cache, it will be skipped.

        Args:
            source_dir: The directory to populate the cache from.

        Returns:
            A dictionary mapping source file paths to their cached paths.

        Raises:
            NotADirectoryError: If the source_dir is not a directory.
        """
        if not source_dir.is_dir():
            raise NotADirectoryError(source_dir)
        cached = {}
        for file_path in source_dir.iterdir():
            if not file_path.is_file():
                continue
            cached_path = self._as_cached_path(file_path.name)
            if cached_path.exists():
                logger.debug(f"File {file_path.name} already in cache. Skipping.")
                continue
            copyfile(file_path, cached_path, copy_method=self.copy_method)
            cached[file_path] = cached_path
        return cached


async def retrieve_files(
    urls: Iterable[tuple[URL | str, str]],
    save_dir: Path,
    max_parallel_downloads: int = 5,
    retries: int = 3,
    total_timeout: int = 300,
    desc: str = "Downloading files",
    cacher: Cacher | None = None,
    chunk_size: int = 524288,  # 512 KiB
    gzip_files: bool = False,
    raise_for_not_found: bool = True,
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
        gzip_files: Whether to gzip the downloaded files.
            This requires the server can send gzip encoded content.
        raise_for_not_found: Whether to raise an error for HTTP 404 errors.
            If false then function does not returns Path for which url gave HTTP 404 error and logs as debug message.

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
                gzip_files=gzip_files,
                raise_for_not_found=raise_for_not_found,
            )
            for url, filename in urls
        ]
        raw_files: list[Path | None] = await tqdm.gather(*tasks, desc=desc)
        return [f for f in raw_files if f is not None]


class InvalidContentEncodingError(aiohttp.ClientResponseError):
    """Content encoding is invalid."""


async def _retrieve_file(
    session: RetryClient,
    url: URL | str,
    save_path: Path,
    semaphore: asyncio.Semaphore,
    cacher: Cacher | None = None,
    chunk_size: int = 524288,  # 512 KiB
    gzip_files: bool = False,
    raise_for_not_found=True,
) -> Path | None:
    """Retrieve a single file from a URL and save it to a specified path.

    Args:
        session: The aiohttp session to use for the request.
        url: The URL to download the file from.
        save_path: The path where the file should be saved.
        semaphore: A semaphore to limit the number of concurrent downloads.
        cacher: An optional cacher to use for caching files.
        chunk_size: The size of each chunk to read from the response.
        gzip_files: Whether to gzip the downloaded file.
            This requires the server can send gzip encoded content.
        raise_for_not_found: Whether to raise an error for HTTP 404 errors.
            If false then function returns None on HTTP 404 errors and logs as debug message.

    Returns:
        The path to the saved file.
    """
    if save_path.exists():
        logger.debug(f"File {save_path} already exists. Skipping download from {url}.")
        return save_path

    if cacher is None:
        cacher = PassthroughCacher()
    if cached_file := await cacher.copy_from_cache(save_path):
        logger.debug(f"File {save_path} was copied from cache {cached_file}. Skipping download from {url}.")
        return save_path

    # Alphafold server and many other web servers can return gzipped responses,
    # when we want to save as *.gz, we use raw stream
    # otherwise aiohttp will decompress it automatically for us.
    auto_decompress = not gzip_files
    headers = {"Accept-Encoding": "gzip"}
    async with (
        semaphore,
        session.get(url, headers=headers, auto_decompress=auto_decompress) as resp,
    ):
        if not raise_for_not_found and resp.status == 404:
            logger.debug(f"File not found at {url}, skipping download.")
            return None
        resp.raise_for_status()
        if gzip_files and resp.headers.get("Content-Encoding") != "gzip":
            msg = f"Server did not send gzip encoded content for {url}, can not save as gzipped file."
            raise InvalidContentEncodingError(
                request_info=resp.request_info,
                history=resp.history,
                status=415,
                message=msg,
                headers=resp.headers,
            )
        iterator = resp.content.iter_chunked(chunk_size)
        await cacher.write_iter(save_path, iterator)
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

    Note that the hardlink copy method only works within the same filesystem and is harder to track.
    If you want to track cached files easily then use 'symlink'.
    On Windows you need developer mode or admin privileges to create symlinks.

    Args:
        source: The source file to copy or link.
        target: The target file to create.
        copy_method: The method to use for copying.

    Raises:
        FileNotFoundError: If the source file or parent of target does not exist.
        FileExistsError: If the target file already exists.
        ValueError: If an unknown copy method is provided.
    """
    if copy_method == "copy":
        shutil.copyfile(source, target)
    elif copy_method == "symlink":
        rel_source = source.absolute().relative_to(target.parent.absolute(), walk_up=True)
        target.symlink_to(rel_source)
    elif copy_method == "hardlink":
        target.hardlink_to(source)
    else:
        msg = f"Unknown method: {copy_method}. Valid methods are: {copy_methods}"
        raise ValueError(msg)


async def async_copyfile(
    source: Path,
    target: Path,
    copy_method: CopyMethod = "copy",
):
    """Asynchronously make target path be same file as source by either copying or symlinking or hardlinking.

    Note that the hardlink copy method only works within the same filesystem and is harder to track.
    If you want to track cached files easily then use 'symlink'.
    On Windows you need developer mode or admin privileges to create symlinks.

    Args:
        source: The source file to copy.
        target: The target file to create.
        copy_method: The method to use for copying.

    Raises:
        FileNotFoundError: If the source file or parent of target does not exist.
        FileExistsError: If the target file already exists.
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
        msg = f"Unknown method: {copy_method}. Valid methods are: {copy_methods}"
        raise ValueError(msg)


def populate_cache_command(raw_args: Sequence[str] | None = None):
    """Command line interface to populate the cache from an existing directory.

    Can be called from the command line as:

    ```bash
    python3 -m protein_quest.utils populate-cache /path/to/source/dir
    ```

    Args:
        raw_args: The raw command line arguments to parse. If None, uses sys.argv.
    """
    root_parser = argparse.ArgumentParser(formatter_class=ArgumentDefaultsRichHelpFormatter)
    subparsers = root_parser.add_subparsers(dest="command")

    desc = "Populate the cache directory with files from the source directory."
    populate_cache_parser = subparsers.add_parser(
        "populate-cache",
        help=desc,
        description=desc,
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    populate_cache_parser.add_argument("source_dir", type=Path)
    populate_cache_parser.add_argument(
        "--cache-dir",
        type=Path,
        default=user_cache_root_dir(),
        help="Directory to use for caching. If not provided, a default cache directory is used.",
    )
    populate_cache_parser.add_argument(
        "--copy-method",
        type=str,
        default="hardlink",
        choices=copy_methods,
        help="Method to use for copying files to cache.",
    )

    args = root_parser.parse_args(raw_args)
    if args.command == "populate-cache":
        source_dir = args.source_dir
        cacher = DirectoryCacher(cache_dir=args.cache_dir, copy_method=args.copy_method)
        cached_files = cacher.populate_cache(source_dir)
        rich.print(f"Cached {len(cached_files)} files from {source_dir} to {cacher.cache_dir}")
        for src, cached in cached_files.items():
            rich.print(f"- {src.relative_to(source_dir)} -> {cached.relative_to(cacher.cache_dir)}")


if __name__ == "__main__":
    populate_cache_command()
