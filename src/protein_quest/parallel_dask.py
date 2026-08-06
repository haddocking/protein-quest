"""Dask-backed parallel helpers."""

import logging
import os
import sys
import warnings
from collections.abc import Callable, Collection, Generator
from contextlib import contextmanager, suppress
from typing import Any, Concatenate, ParamSpec, cast, override

from dask.distributed import Client, LocalCluster
from distributed.deploy.cluster import Cluster
from distributed.diagnostics.progress import format_time
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from tornado.ioloop import IOLoop

from protein_quest.parallel_common import MapWithProgressOptions

HAS_DISTRIBUTED = True

try:
    from psutil import cpu_count
except ImportError:

    def cpu_count(logical: bool = True) -> int | None:  # noqa: ARG001
        return None


logger = logging.getLogger(__name__)
P = ParamSpec("P")


@contextmanager
def configure_dask_scheduler(
    scheduler_address: str | Cluster | None,
    name: str,
    nproc: int = 1,
) -> Generator[str | Cluster]:
    """Context manager that offers a Dask cluster."""
    if scheduler_address is not None:
        yield scheduler_address
        return
    cluster = _configure_cpu_dask_scheduler(nproc, name)
    logger.info(f"Using local Dask cluster: {cluster}")
    try:
        yield cluster
    finally:
        cluster.close()


def nr_cpus() -> int:
    """Determine the number of CPU cores to use."""
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
    return LocalCluster(name=name, threads_per_worker=1, n_workers=n_workers)


class MyProgressBar(ProgressBar):
    """Show progress of Dask computations."""

    __loop: IOLoop | None = None

    def __init__(
        self,
        keys: Any,
        scheduler: Any | None = None,
        interval: str = "100ms",
        width: int = 40,
        loop: Any | None = None,
        complete: bool = True,
        start: bool = True,
        **kwargs: Any,  # noqa: ARG002
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
            self.__loop = loop = self._loop_runner.loop
        return loop

    @loop.setter
    def loop(self, value: IOLoop) -> None:
        warnings.warn("setting the loop property is deprecated", DeprecationWarning, stacklevel=2)
        self.__loop = value

    def _draw_bar(self, remaining: int, all: int, **kwargs: Any):  # noqa: A002, ARG002
        frac = (1 - remaining / all) if all else 1.0
        bar = "#" * int(self.width * frac)
        percent = int(100 * frac)
        elapsed = format_time(self.elapsed)
        msg = "\r[{0:<{1}}] | {2}% Completed | {3}".format(bar, self.width, percent, elapsed)
        with suppress(ValueError):
            sys.stderr.write(msg)
            sys.stderr.flush()

    @override
    def _draw_stop(self, **_kwargs: Any):
        sys.stderr.write("\33[2K\r")
        sys.stderr.flush()


def _apply_with_shared_arguments[T, R](
    item: T,
    worker_func: Callable[..., R],
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
) -> R:
    return worker_func(item, *args, **kwargs)


def dask_map_with_progress[T, R, **P](
    client: Client,
    func: Callable[Concatenate[T, P], R],
    iterable: Collection[T],
    *args: P.args,
    **kwargs: P.kwargs,
) -> list[R]:
    """Map over an iterable with Dask, progress reporting, and a typed result."""
    if client.dashboard_link:
        logger.info(f"Follow progress on dask dashboard at: {client.dashboard_link}")
    shared_args = client.scatter([tuple(args)], broadcast=True)[0]
    shared_kwargs = client.scatter([dict(kwargs)], broadcast=True)[0]
    futures = client.map(
        _apply_with_shared_arguments,
        iterable,
        worker_func=func,
        args=shared_args,
        kwargs=shared_kwargs,
    )
    if not os.getenv("TQDM_DISABLE"):
        MyProgressBar(futures)
    results = client.gather(futures)
    return cast("list[R]", results)


def map_dask_with_progress[T, R, **P](
    scheduler_address: str | Cluster | None,
    func: Callable[Concatenate[T, P], R],
    iterable: Collection[T],
    map_with_progress_options: MapWithProgressOptions | None = None,
    *args: P.args,
    **kwargs: P.kwargs,
) -> list[R]:
    """Map a function over an iterable using Dask with scheduler/client setup."""
    opts = map_with_progress_options or {}
    scheduler_name = opts.get("dask_scheduler_name", "map_with_progress")
    with configure_dask_scheduler(scheduler_address, name=scheduler_name) as cluster, Client(cluster) as client:
        client.forward_logging()
        return dask_map_with_progress(client, func, iterable, *args, **kwargs)
