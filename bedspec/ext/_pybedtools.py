try:
    from pybedtools import Interval
except ImportError as error:
    raise ImportError(
        "bedspec.ext.pybedtools needs pybedtools: pip install 'bedspec[pybedtools]'"
    ) from error

from bedspec import Bed6
from bedspec.ext._stranded import stranded


def as_interval(interval: Interval) -> list[Bed6]:
    """Return the interval of a pybedtools interval, on its strand if it has one."""
    return [stranded(interval.chrom, interval.start, interval.end, interval.strand)]
