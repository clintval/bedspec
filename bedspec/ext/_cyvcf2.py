try:
    from cyvcf2 import Variant
except ImportError as error:
    raise ImportError("bedspec.ext.cyvcf2 needs cyvcf2: pip install 'bedspec[cyvcf2]'") from error

from bedspec import Bed3


def as_interval(variant: Variant) -> list[Bed3]:
    """Return the interval of a variant, its reference allele."""
    return [Bed3(variant.CHROM, start=variant.start, end=variant.end)]
