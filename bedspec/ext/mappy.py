"""Functions reading mappy's records as spans, to register with `register_locatable`."""

from bedspec.ext._mappy import as_blocks
from bedspec.ext._mappy import as_interval

__all__ = [
    "as_blocks",
    "as_interval",
]
