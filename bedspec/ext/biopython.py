"""Functions reading Biopython's records as spans, to register with `register_locatable`."""

from bedspec.ext._biopython import as_blocks
from bedspec.ext._biopython import as_interval

__all__ = [
    "as_blocks",
    "as_interval",
]
