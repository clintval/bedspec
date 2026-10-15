"""Functions reading Biopython's records as spans, to register with `register_for_bedspec`."""

from bedspec.ext._biopython import alignment_blocks

__all__ = [
    "alignment_blocks",
]
