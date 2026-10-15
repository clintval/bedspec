from bedspec import Bed6
from bedspec import BedStrand

_STRANDS: dict[object, BedStrand] = {
    "+": BedStrand.Positive,
    "-": BedStrand.Negative,
    1: BedStrand.Positive,
    -1: BedStrand.Negative,
}


def stranded(refname: str, start: int, end: int, strand: object) -> Bed6:
    """Return a span on a strand given as `+`, `-`, `1` or `-1`, or on none for anything else."""
    return Bed6(refname, start=start, end=end, name=None, score=None, strand=_STRANDS.get(strand))
