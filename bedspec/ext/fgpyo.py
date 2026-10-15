"""Functions reading fgpyo's records as spans, to register with `register_for_bedspec`."""

from bedspec.ext._fgpyo import supplementary_span
from bedspec.ext._fgpyo import template_blocks

__all__ = [
    "supplementary_span",
    "template_blocks",
]
