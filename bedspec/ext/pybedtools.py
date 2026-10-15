"""Functions reading pybedtools' records as spans, to register with `register_for_bedspec`."""

from bedspec.ext._pybedtools import interval_span

__all__ = [
    "interval_span",
]
