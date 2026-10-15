"""Functions reading mappy's records as spans, to register with `register_for_bedspec`."""

try:
    from mappy import Alignment
except ImportError as error:
    raise ImportError("bedspec.ext.mappy needs mappy: pip install 'bedspec[mappy]'") from error

from bedspec import Bed6
from bedspec.ext._stranded import stranded

_ALIGNED = (0, 7, 8)
_ON_REFERENCE = (0, 2, 3, 7, 8)


def alignment_blocks(alignment: Alignment) -> list[Bed6]:
    """Return the aligned blocks of an alignment, split at introns and deletions, on its strand."""
    blocks: list[Bed6] = []
    start = alignment.r_st
    for length, operator in alignment.cigar:
        if operator in _ALIGNED:
            blocks.append(stranded(alignment.ctg, start, start + length, alignment.strand))
        if operator in _ON_REFERENCE:
            start += length
    return blocks
