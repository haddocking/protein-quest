"""Dask helper functions."""

import logging
import os
from collections.abc import Callable, Collection
from typing import Concatenate, ParamSpec, cast

from dask.distributed import Client, LocalCluster, progress
from distributed.deploy.cluster import Cluster
from psutil import cpu_count

logger = logging.getLogger(__name__)


def configure_dask_scheduler(
    scheduler_address: str | Cluster | None,
    name: str,
    nproc: int = 1,
) -> str | Cluster:
    """Configure the Dask scheduler by reusing existing or creating a new cluster.

    Args:
        scheduler_address: Address of the Dask scheduler to connect to, or None for local cluster.
        name: Name for the Dask cluster.
        nproc: Number of processes to use per worker for CPU support.

    Returns:
        A Dask Cluster instance or a string address for the scheduler.
    """
    if scheduler_address is None:
        scheduler_address = _configure_cpu_dask_scheduler(nproc, name)
        logger.info(f"Using local Dask cluster: {scheduler_address}")

    return scheduler_address


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
    progress(futures)
    results = client.gather(futures)
    return cast("list[R]", results)
