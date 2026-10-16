"""Functions reading pybedtools' records as spans, to register with `register_locatable`."""

from bedspec.ext._pybedtools import as_interval

__all__ = [
    "as_interval",
]
