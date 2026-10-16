"""Functions reading fgpyo's records as spans, to register with `register_locatable`."""

from bedspec.ext._fgpyo import as_interval
from bedspec.ext._fgpyo import as_template_blocks
from bedspec.ext._fgpyo import as_template_intervals

__all__ = [
    "as_interval",
    "as_template_blocks",
    "as_template_intervals",
]
