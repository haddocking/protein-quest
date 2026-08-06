"""Common types and options for parallel execution."""

from typing import TypedDict


class MapWithProgressOptions(TypedDict, total=False):
    """Options for progress bar and Dask scheduler in [map_with_progress][protein_quest.parallel.map_with_progress].

    Attributes:
        dask_scheduler_name: Name for the Dask scheduler (default: "map_with_progress").
        tqdm_desc: Description for the tqdm progress bar (default: "").
        tqdm_unit: Unit for the tqdm progress bar (default: "it").
    """

    dask_scheduler_name: str
    tqdm_desc: str
    tqdm_unit: str
