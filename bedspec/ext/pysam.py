"""Functions reading pysam's records as spans, to register with `register_locatable`."""

from bedspec.ext._pysam import as_blocks
from bedspec.ext._pysam import as_interval
from bedspec.ext._pysam import as_paired_blocks
from bedspec.ext._pysam import as_paired_intervals

__all__ = [
    "as_blocks",
    "as_interval",
    "as_paired_blocks",
    "as_paired_intervals",
]
