"""Functions reading cyvcf2's records as spans, to register with `register_for_bedspec`."""

try:
    from cyvcf2 import Variant
except ImportError as error:
    raise ImportError("bedspec.ext.cyvcf2 needs cyvcf2: pip install 'bedspec[cyvcf2]'") from error

from bedspec import Bed3


def variant_span(variant: Variant) -> list[Bed3]:
    """Return the reference bases a variant spans."""
    return [Bed3(variant.CHROM, start=variant.start, end=variant.end)]
