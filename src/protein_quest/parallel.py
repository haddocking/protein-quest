"""Dask helper functions."""

import logging
import os
import sys
import warnings
from collections.abc import Callable, Collection, Iterator
from contextlib import contextmanager, suppress
from typing import Concatenate, ParamSpec, cast

from dask.distributed import Client, LocalCluster
from distributed.deploy.cluster import Cluster
from distributed.diagnostics.progress import format_time
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from psutil import cpu_count
from tornado.ioloop import IOLoop

logger = logging.getLogger(__name__)


@contextmanager
def configure_dask_scheduler(
    scheduler_address: str | Cluster | None,
    name: str,
    nproc: int = 1,
) -> Iterator[str | Cluster]:
    """Context manager that offers a Dask cluster.

    If scheduler_address is None then creates a local Dask cluster
    else returns scheduler_address unchanged and the callee is responsible for cluster cleanup.

    Args:
        scheduler_address: Address of the Dask scheduler to connect to, or None for local cluster.
        name: Name for the Dask cluster.
        nproc: Number of processes to use per worker for CPU support.

    Yields:
        The scheduler address as a string or a cluster.
    """
    if scheduler_address is not None:
        # Pass through existing scheduler address or cluster
        yield scheduler_address
        return
    cluster = _configure_cpu_dask_scheduler(nproc, name)
    logger.info(f"Using local Dask cluster: {cluster}")
    try:
        yield cluster
    finally:
        cluster.close()


def nr_cpus() -> int:
    """Determine the number of CPU cores to use.

    If the environment variables SLURM_CPUS_PER_TASK or OMP_NUM_THREADS are set,
    their value is used. Otherwise, the number of physical CPU cores is returned.

    Returns:
        The number of CPU cores to use.

    Raises:
        ValueError: If the number of physical CPU cores cannot be determined.
    """
    physical_cores = cpu_count(logical=False)
    if physical_cores is None:
        msg = "Cannot determine number of logical CPU cores."
        raise ValueError(msg)
    for var in ["SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS"]:
        value = os.environ.get(var)
        if value is not None:
            logger.warning(
                'Not using all CPU cores (%s) of machine, environment variable "%s" is set to %s.',
                physical_cores,
                var,
                value,
            )
            return int(value)
    return physical_cores


def _configure_cpu_dask_scheduler(nproc: int, name: str) -> LocalCluster:
    total_cpus = nr_cpus()
    n_workers = total_cpus // nproc
    # Use single thread per worker to prevent GIL slowing down the computations
    return LocalCluster(name=name, threads_per_worker=1, n_workers=n_workers)


class MyProgressBar(ProgressBar):
    """Show progress of Dask computations.

    Copy of distributed.diagnostics.progressbar.TextProgressBar that:

    - prints to stderr instead of stdout
    - Can have its interval (in seconds) set with `TQDM_MININTERVAL` environment variable

    """

    __loop: IOLoop | None = None

    def __init__(
        self,
        keys,
        scheduler=None,
        interval="100ms",
        width=40,
        loop=None,
        complete=True,
        start=True,
        **kwargs,  # noqa: ARG002
    ):
        self._loop_runner = loop_runner = LoopRunner(loop=loop)
        if interval == "100ms":
            interval_env = os.getenv("TQDM_MININTERVAL")
            if interval_env is not None:
                interval = interval_env + "s"

        super().__init__(keys, scheduler, interval, complete)
        self.width = width

        if start:
            loop_runner.run_sync(self.listen)

    @property
    def loop(self) -> IOLoop | None:
        loop = self.__loop
        if loop is None:
            # If the loop is not running when this is called, the LoopRunner.loop
            # property will raise a DeprecationWarning
            # However subsequent calls might occur - eg atexit, where a stopped
            # loop is still acceptable - so we cache access to the loop.
            self.__loop = loop = self._loop_runner.loop
        return loop

    @loop.setter
    def loop(self, value: IOLoop) -> None:
        warnings.warn("setting the loop property is deprecated", DeprecationWarning, stacklevel=2)
        self.__loop = value

    def _draw_bar(self, remaining, all, **kwargs):  # noqa: A002, ARG002
        frac = (1 - remaining / all) if all else 1.0
        bar = "#" * int(self.width * frac)
        percent = int(100 * frac)
        elapsed = format_time(self.elapsed)
        msg = "\r[{0:<{1}}] | {2}% Completed | {3}".format(bar, self.width, percent, elapsed)
        with suppress(ValueError):
            sys.stderr.write(msg)
            sys.stderr.flush()

    def _draw_stop(self, **kwargs):  # noqa: ARG002
        sys.stderr.write("\33[2K\r")
        sys.stderr.flush()


# Generic type parameters used across helpers
P = ParamSpec("P")


def dask_map_with_progress[T, R, **P](
    client: Client,
    func: Callable[Concatenate[T, P], R],
    iterable: Collection[T],
    *args: P.args,
    **kwargs: P.kwargs,
) -> list[R]:
    """
    Wrapper for map, progress, and gather of Dask that returns a correctly typed list.

    Environment variables:
        - Set interval (in seconds) of progress updates with `TQDM_MININTERVAL`
        - Disabled by setting `TQDM_DISABLE` to any value

    Args:
        client: Dask client.
        func: Function to map; first parameter comes from ``iterable`` and any
            additional parameters can be provided positionally via ``*args`` or
            as keyword arguments via ``**kwargs``.
        iterable: Collection of arguments to map over.
        *args: Additional positional arguments to pass to client.map().
        **kwargs: Additional keyword arguments to pass to client.map().

    Returns:
        List of results of type returned by `func` function.
    """
    if client.dashboard_link:
        logger.info(f"Follow progress on dask dashboard at: {client.dashboard_link}")
    futures = client.map(func, iterable, *args, **kwargs)
    if not os.getenv("TQDM_DISABLE"):
        MyProgressBar(futures)
    results = client.gather(futures)
    return cast("list[R]", results)
