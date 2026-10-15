"""Functions reading pysam's records as spans, to register with `register_locatable`."""

from bedspec.ext._pysam import pair_blocks
from bedspec.ext._pysam import read_blocks
from bedspec.ext._pysam import read_span
from bedspec.ext._pysam import tabix_span
from bedspec.ext._pysam import variant_span

__all__ = [
    "pair_blocks",
    "read_blocks",
    "read_span",
    "tabix_span",
    "variant_span",
]
