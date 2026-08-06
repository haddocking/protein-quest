"""Parallel execution helpers with sequential fallback."""

import warnings
from collections.abc import Callable, Collection
from typing import Any, Concatenate, Literal, ParamSpec

from tqdm.auto import tqdm

from protein_quest.parallel_common import MapWithProgressOptions

try:
    from protein_quest.parallel_dask import Cluster, map_dask_with_progress
except (AttributeError, ImportError):
    Cluster = Any

    def map_dask_with_progress[T, R, **P](
        _scheduler_address: str | Cluster | None,
        func: Callable[Concatenate[T, P], R],
        iterable: Collection[T],
        map_with_progress_options: MapWithProgressOptions | None = None,
        *args: P.args,
        **kwargs: P.kwargs,
    ) -> list[R]:
        warnings.warn(
            "distributed is not available; falling back to sequential execution.",
            RuntimeWarning,
            stacklevel=2,
        )
        return _map_sequential_with_progress(func, iterable, map_with_progress_options, *args, **kwargs)


SchedulerAddress = str | Cluster | Literal["sequential"] | None
"""Dask scheduler address.

* "sequential": Run sequentially without Dask.
* None: Create a local Dask cluster.
* str: Address of an existing Dask scheduler.
* Cluster: An existing Dask cluster object.
"""

P = ParamSpec("P")


def _map_sequential_with_progress[T, R, **P](
    func: Callable[Concatenate[T, P], R],
    iterable: Collection[T],
    map_with_progress_options: MapWithProgressOptions | None,
    *args: P.args,
    **kwargs: P.kwargs,
) -> list[R]:
    """Map a function over an iterable with a tqdm progress bar.

    Args:
        func: Function to map.
        iterable: Collection of items to map over.
        map_with_progress_options: Options for progress bar.
        *args: Positional arguments to pass to the function.
        **kwargs: Keyword arguments to pass to the function.

    Returns:
        List of results from applying the function to each item in the iterable.
    """
    opts = map_with_progress_options or {}
    desc = opts.get("tqdm_desc", "")
    unit = opts.get("tqdm_unit", "it")
    return [func(x, *args, **kwargs) for x in tqdm(iterable, desc=desc, unit=unit)]


def map_with_progress[T, R, **P](
    scheduler_address: SchedulerAddress,
    func: Callable[Concatenate[T, P], R],
    iterable: Collection[T],
    map_with_progress_options: MapWithProgressOptions | None = None,
    *args: P.args,
    **kwargs: P.kwargs,
) -> list[R]:
    """Map a function over an iterable with optional progress bar.

    Wraps sequential execution with tqdm and parallel execution with Dask.

    Args:
        scheduler_address: "sequential" for local execution, None for local cluster, or existing cluster.
        func: Function to map.
        iterable: Collection of items to map over.
        map_with_progress_options: Options for progress bar and Dask scheduler.
        *args: Positional arguments passed to func.
        **kwargs: Keyword arguments passed to func.

    Returns:
        List of results from func.
    """
    if scheduler_address == "sequential":
        return _map_sequential_with_progress(func, iterable, map_with_progress_options, *args, **kwargs)
    return map_dask_with_progress(scheduler_address, func, iterable, map_with_progress_options, *args, **kwargs)
