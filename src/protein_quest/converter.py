"""Convert json or dict to Python objects."""

from cattrs.preconf.orjson import make_converter
from yarl import URL

type Percentage = float
"""Type alias for percentage values (0.0-100.0)."""
type Ratio = float
"""Type alias for ratio values (0.0-1.0)."""
type PositiveInt = int
"""Type alias for positive integer values (>= 0)."""

converter = make_converter()
"""cattrs converter to read JSON document or dict to Python objects."""
converter.register_structure_hook(URL, lambda v, _: URL(v))
converter.register_unstructure_hook(URL, lambda u: str(u))


@converter.register_structure_hook
def percentage_hook(val, _) -> Percentage:
    value = float(val)
    """Cattrs hook to validate percentage values."""
    if not 0.0 <= value <= 100.0:
        msg = f"Value {value} is not a valid percentage (0.0-100.0)"
        raise ValueError(msg)
    return value


@converter.register_structure_hook
def ratio_hook(val, _) -> Ratio:
    """Cattrs hook to validate ratio values."""
    value = float(val)
    if not 0.0 <= value <= 1.0:
        msg = f"Value {value} is not a valid ratio (0.0-1.0)"
        raise ValueError(msg)
    return value


@converter.register_structure_hook
def positive_int_hook(val, _) -> PositiveInt:
    """Cattrs hook to validate positive integer values."""
    value = int(val)
    if value < 0:
        msg = f"Value {value} is not a valid positive integer (>= 0)"
        raise ValueError(msg)
    return value
