"""Functions reading cyvcf2's records as spans, to register with `register_locatable`."""

from bedspec.ext._cyvcf2 import variant_span

__all__ = [
    "variant_span",
]
