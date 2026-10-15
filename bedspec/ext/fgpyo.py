"""Functions reading fgpyo's records as spans, to register with `register_locatable`."""

from bedspec.ext._fgpyo import as_interval
from bedspec.ext._fgpyo import as_paired_blocks
from bedspec.ext._fgpyo import as_paired_intervals

__all__ = [
    "as_interval",
    "as_paired_blocks",
    "as_paired_intervals",
]
