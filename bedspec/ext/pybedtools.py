"""Functions reading pybedtools' records as spans, to register with `register_locatable`."""

from bedspec.ext._pybedtools import interval_span

__all__ = [
    "interval_span",
]
