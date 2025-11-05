import asyncio
import gzip
import logging
from pathlib import Path

import pytest
from aiohttp.client import ClientResponseError
from aiohttp.streams import AsyncStreamIterator
from aiohttp.test_utils import TestServer
from aiohttp.web import Application
from aiohttp_compress import compress_middleware
from pytest_aiohttp import AiohttpServer

from protein_quest.utils import (
    CopyMethod,
    DirectoryCacher,
    InvalidContentEncodingError,
    PassthroughCacher,
    async_copyfile,
    copyfile,
    populate_cache_command,
    retrieve_files,
    user_cache_root_dir,
)


def test_copyfile_copy(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    copyfile(src, dst, "copy")
    assert dst.read_text() == "Hello, World!"
    assert not dst.is_symlink()


def test_copyfile_symlink(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    copyfile(src, dst, "symlink")
    assert dst.read_text() == "Hello, World!"
    assert dst.is_symlink()
    assert dst.resolve() == src


def test_copyfile_hardlink(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    copyfile(src, dst, "hardlink")
    assert dst.read_text() == "Hello, World!"
    assert not dst.is_symlink()
    assert src.stat().st_ino == dst.stat().st_ino
    assert src.stat().st_nlink == 2
    assert dst.stat().st_nlink == 2


def test_copyfile_invalid_method(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    with pytest.raises(ValueError, match="Unknown method"):
        copyfile(src, dst, "invalid")  # type: ignore  # noqa: PGH003


@pytest.mark.asyncio
async def test_async_copyfile_copy(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    await async_copyfile(src, dst, "copy")
    assert dst.read_text() == "Hello, World!"
    assert not dst.is_symlink()


@pytest.mark.asyncio
async def test_async_copyfile_symlink(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    await async_copyfile(src, dst, "symlink")
    assert dst.read_text() == "Hello, World!"
    assert dst.is_symlink()
    assert dst.resolve() == src


@pytest.mark.asyncio
async def test_async_copyfile_hardlink(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("Hello, World!")
    await async_copyfile(src, dst, "hardlink")
    assert dst.read_text() == "Hello, World!"
    assert not dst.is_symlink()
    assert src.stat().st_ino == dst.stat().st_ino
    assert src.stat().st_nlink == 2
    assert dst.stat().st_nlink == 2


@pytest.mark.asyncio
async def test_async_copyfile_invalid_method(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    with pytest.raises(ValueError, match="Unknown method"):
        await async_copyfile(src, dst, "invalid")  # type: ignore  # noqa: PGH003


def test_user_cache_root_dir():
    cache_dir = user_cache_root_dir()
    assert cache_dir.name == "protein-quest"


class ByteGenerator(AsyncStreamIterator[bytes]):
    """Mock for return of aiohttp.ClientResponse.content.iter_chunked()"""

    def __init__(self, data: bytes, chunk_size: int = 5):
        self.data = data
        self.chunk_size = chunk_size

    def __aiter__(self):
        self.index = 0
        return self

    async def __anext__(self):
        if self.index >= len(self.data):
            raise StopAsyncIteration
        chunk = self.data[self.index : self.index + self.chunk_size]
        self.index += self.chunk_size
        await asyncio.sleep(0)  # Yield control to the event loop
        return chunk


def test_DirectoryCacher_init_copy_method_invalid(tmp_path: Path):
    with pytest.raises(ValueError, match="Unknown copy method"):
        DirectoryCacher(tmp_path / "cache", copy_method="invalid")  # type: ignore  # noqa: PGH003


@pytest.mark.parametrize("copy_method", ["copy", "symlink", "hardlink"])
class TestDirectoryCacher:
    @pytest.fixture
    def cacher(self, tmp_path: Path, copy_method: CopyMethod) -> DirectoryCacher:
        return DirectoryCacher(tmp_path / "cache", copy_method=copy_method)

    def test_init(self, tmp_path: Path, cacher: DirectoryCacher):
        assert cacher.cache_dir == tmp_path / "cache"
        assert cacher.cache_dir.exists()
        assert cacher.cache_dir.is_dir()

    @pytest.mark.asyncio
    async def test_write_bytes(self, tmp_path: Path, cacher: DirectoryCacher):
        target = tmp_path / "test.txt"

        cache_file = await cacher.write_bytes(target, b"Hello, World!")

        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

    def test_in_missing(self, cacher: DirectoryCacher):
        assert "any_file.txt" not in cacher

    @pytest.mark.asyncio
    async def test_in_present(self, tmp_path: Path, cacher: DirectoryCacher):
        # Fill cache
        target = tmp_path / "test.txt"
        await cacher.write_bytes(target, b"Hello, World!")

        assert target in cacher

    @pytest.mark.asyncio
    async def test_copy_from_cache_missing(self, tmp_path: Path, cacher: DirectoryCacher):
        target = tmp_path / "test.txt"

        result = await cacher.copy_from_cache(target)
        assert result is None

    @pytest.mark.asyncio
    async def test_copy_from_cache_present(self, tmp_path: Path, cacher: DirectoryCacher):
        other = tmp_path / "other"
        other.mkdir()
        source = other / "test.txt"
        await cacher.write_bytes(source, b"Hello, World!")

        target = tmp_path / "test.txt"
        result = await cacher.copy_from_cache(target)

        assert result is not None
        assert result.exists()
        assert result.read_bytes() == b"Hello, World!"

    @pytest.mark.asyncio
    async def test_write_iter(self, tmp_path: Path, cacher: DirectoryCacher):
        target = tmp_path / "test.txt"

        cache_file = await cacher.write_iter(target, ByteGenerator(b"Hello, World!"))

        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

    @pytest.mark.asyncio
    async def test_write_bytes_twice(self, tmp_path: Path, cacher: DirectoryCacher):
        target = tmp_path / "test.txt"

        # First write
        cache_file = await cacher.write_bytes(target, b"Hello, World!")
        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

        # Second write
        with pytest.raises(FileExistsError, match=str(cache_file)):
            await cacher.write_bytes(target, b"Goodbye, World!")

    @pytest.mark.asyncio
    async def test_write_iter_twice(self, tmp_path: Path, cacher: DirectoryCacher):
        target = tmp_path / "test.txt"

        # First write
        cache_file = await cacher.write_iter(target, ByteGenerator(b"Hello, World!"))
        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

        # Second write
        with pytest.raises(FileExistsError, match=str(cache_file)):
            await cacher.write_iter(target, ByteGenerator(b"Goodbye, World!"))


class TestPassthroughCacher:
    @pytest.mark.asyncio
    async def test_write_bytes(self, tmp_path: Path):
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"

        cache_file = await cacher.write_bytes(target, b"Hello, World!")

        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

    def test_in_missing(self, tmp_path: Path):
        cacher = PassthroughCacher()
        assert "any_file.txt" not in cacher

    @pytest.mark.asyncio
    async def test_in_present(self, tmp_path: Path):
        # Fill cache
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"
        await cacher.write_bytes(target, b"Hello, World!")

        # PassthroughCacher never has anything cached
        assert target not in cacher

    @pytest.mark.asyncio
    async def test_copy_from_cache_missing(self, tmp_path: Path):
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"

        result = await cacher.copy_from_cache(target)
        assert result is None

    @pytest.mark.asyncio
    async def test_copy_from_cache_present(self, tmp_path: Path):
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"
        await cacher.write_bytes(target, b"Hello, World!")

        result = await cacher.copy_from_cache(target)

        # PassthroughCacher never has anything cached
        assert result is None

    @pytest.mark.asyncio
    async def test_write_iter(self, tmp_path: Path):
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"

        cache_file = await cacher.write_iter(target, ByteGenerator(b"Hello, World!"))

        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

    @pytest.mark.asyncio
    async def test_write_bytes_twice(self, tmp_path: Path):
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"

        # First write
        cache_file = await cacher.write_bytes(target, b"Hello, World!")
        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

        # Second write
        with pytest.raises(FileExistsError, match=str(target)):
            await cacher.write_bytes(target, b"Goodbye, World!")

    @pytest.mark.asyncio
    async def test_write_iter_twice(self, tmp_path: Path):
        cacher = PassthroughCacher()
        target = tmp_path / "test.txt"

        # First write
        cache_file = await cacher.write_iter(target, ByteGenerator(b"Hello, World!"))
        assert cache_file.exists()
        assert cache_file.read_bytes() == b"Hello, World!"

        # Second write
        with pytest.raises(FileExistsError, match=str(target)):
            await cacher.write_iter(target, ByteGenerator(b"Goodbye, World!"))


def test_populate_cache_command_with_hardlink(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    src1 = source_dir / "file1.txt"
    src2 = source_dir / "file2.txt"
    src1.write_text("File 1")
    src2.write_text("File 2")

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    populate_cache_command(
        [
            "populate-cache",
            str(source_dir),
            "--cache-dir",
            str(cache_dir),
            "--copy-method",
            "hardlink",
        ]
    )

    assert src1.stat().st_nlink == 2
    assert src2.stat().st_nlink == 2

    captured = capsys.readouterr()
    assert "->" in captured.out


async def static_server(tmp_path: Path, aiohttp_server: AiohttpServer) -> TestServer:
    host_dir = tmp_path / "host"
    host_dir.mkdir()
    (host_dir / "file1.txt").write_text("This is file 1.")
    app = Application()
    app.router.add_static("/", str(host_dir))
    return await aiohttp_server(app)


@pytest.mark.asyncio
async def test_retrieve_files(tmp_path: Path, aiohttp_server: AiohttpServer):
    server = await static_server(tmp_path, aiohttp_server)
    urls = [
        (server.make_url("file1.txt"), "file1.txt"),
    ]
    save_dir = tmp_path / "downloads"

    downloaded_files = await retrieve_files(
        urls,
        save_dir=save_dir,
    )

    assert downloaded_files == [save_dir / "file1.txt"]
    assert (save_dir / "file1.txt").exists()
    assert (save_dir / "file1.txt").read_text() == "This is file 1."


@pytest.mark.asyncio
async def test_retrieve_files_no_raise_on_not_found(
    tmp_path: Path,
    aiohttp_server: AiohttpServer,
    caplog: pytest.LogCaptureFixture,
):
    server = await static_server(tmp_path, aiohttp_server)
    url = server.make_url("file2.txt")
    urls = [
        (url, "file2.txt"),
    ]
    save_dir = tmp_path / "downloads"

    with caplog.at_level(logging.DEBUG, logger="protein_quest.utils"):
        downloaded_files = await retrieve_files(
            urls,
            save_dir=save_dir,
            raise_for_not_found=False,
        )

    assert downloaded_files == []
    assert not (save_dir / "file2.txt").exists()
    assert f"File not found at {url}, skipping download." in caplog.messages


@pytest.mark.asyncio
async def test_retrieve_files_raise_on_not_found(
    tmp_path: Path,
    aiohttp_server: AiohttpServer,
):
    server = await static_server(tmp_path, aiohttp_server)
    url = server.make_url("file2.txt")
    urls = [
        (url, "file2.txt"),
    ]
    save_dir = tmp_path / "downloads"

    with pytest.raises(ClientResponseError) as exc_info:
        await retrieve_files(
            urls,
            save_dir=save_dir,
            raise_for_not_found=True,
        )

    assert exc_info.value.status == 404
    assert not (save_dir / "file2.txt").exists()


@pytest.mark.asyncio
async def test_retrieve_files_already_exists(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
):
    save_dir = tmp_path
    existing_file = save_dir / "file1.txt"
    existing_file.write_text("Existing content.")
    urls = [
        ("https://example.com/file1.txt", "file1.txt"),
    ]

    with caplog.at_level(logging.DEBUG, logger="protein_quest.utils"):
        downloaded_files = await retrieve_files(
            urls,
            save_dir=save_dir,
        )

    assert downloaded_files == [existing_file]
    assert existing_file.read_text() == "Existing content."
    assert "already exists. Skipping download" in caplog.text


@pytest.mark.asyncio
async def test_retrieve_files_with_directory_cacher_hit(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cacher = DirectoryCacher(cache_dir, copy_method="copy")

    # Pre-fill cache
    cached_file = tmp_path / "cached_file.txt"
    await cacher.write_bytes(cached_file, b"Cached content.")

    save_dir = tmp_path / "downloads"
    urls = [
        ("https://example.com/cached_file.txt", "cached_file.txt"),
    ]

    with caplog.at_level(logging.DEBUG, logger="protein_quest.utils"):
        downloaded_files = await retrieve_files(
            urls,
            save_dir=save_dir,
            cacher=cacher,
        )

    expected_file = save_dir / "cached_file.txt"
    assert downloaded_files == [expected_file]
    assert expected_file.exists()
    assert expected_file.read_bytes() == b"Cached content."
    assert "was copied from cache" in caplog.text
    assert "Skipping download from" in caplog.text


@pytest.mark.asyncio
async def test_retrieve_files_gzipped_on_unsupported_server(
    tmp_path: Path,
    aiohttp_server: AiohttpServer,
):
    server = await static_server(tmp_path, aiohttp_server)
    urls = [
        (server.make_url("file1.txt"), "file1.txt.gz"),
    ]
    save_dir = tmp_path / "downloads"
    save_dir.mkdir()

    with pytest.raises(InvalidContentEncodingError) as exc_info:
        await retrieve_files(
            urls,
            save_dir=save_dir,
            gzip_files=True,
        )

    assert not (save_dir / "file1.txt.gz").exists()
    assert "Server did not send gzip encoded content" in str(exc_info.value)


@pytest.mark.asyncio
async def test_retrieve_files_gzipped(
    tmp_path: Path,
    aiohttp_server: AiohttpServer,
):
    host_dir = tmp_path / "host"
    host_dir.mkdir()
    (host_dir / "file1.txt").write_text("This is file 1.")
    app = Application()
    app.router.add_static("/", str(host_dir))
    app.middlewares.append(compress_middleware)  # type: ignore  # noqa: PGH003
    server = await aiohttp_server(app)

    urls = [
        (server.make_url("file1.txt"), "file1.txt.gz"),
    ]
    save_dir = tmp_path / "downloads"
    save_dir.mkdir()

    downloaded_files = await retrieve_files(
        urls,
        save_dir=save_dir,
        gzip_files=True,
    )

    expected_file = save_dir / "file1.txt.gz"
    assert downloaded_files == [expected_file]
    assert expected_file.exists()
    with gzip.open(expected_file, "rt") as f:
        content = f.read()
        assert content == "This is file 1."
