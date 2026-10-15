"""Functions reading cyvcf2's records as spans, to register with `register_locatable`."""

from bedspec.ext._cyvcf2 import as_interval

__all__ = [
    "as_interval",
]
