try:
    from pybwa import AuxHit
except ImportError as error:
    raise ImportError("bedspec.ext.pybwa needs pybwa: pip install 'bedspec[pybwa]'") from error

from bedspec import Bed6
from bedspec.ext._stranded import stranded


def aux_hit_span(hit: AuxHit) -> list[Bed6]:
    """Return the span of an alternative hit from an XA tag, from 1-based, on its strand."""
    return [stranded(hit.refname, hit.start - 1, hit.end, "-" if hit.negative else "+")]
