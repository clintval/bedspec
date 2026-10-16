try:
    from Bio.Align import Alignment
    from Bio.SeqRecord import SeqRecord
except ImportError as error:
    raise ImportError(
        "bedspec.ext.biopython needs Biopython: pip install 'bedspec[biopython]'"
    ) from error

from bedspec import Bed6
from bedspec.ext._blocks import joined_blocks
from bedspec.ext._stranded import stranded

_UNKNOWN_ID = "<unknown id>"
"""The id Biopython gives a `SeqRecord` made without one."""


def _refname(alignment: Alignment) -> str:
    """Return the name of an alignment's target, the reference its spans are on.

    Raises:
        ValueError: If the target is not a `SeqRecord`, or has no id but Biopython's placeholder.
    """
    target = alignment.target
    if not isinstance(target, SeqRecord):
        kind = type(target).__name__
        raise ValueError(
            f"An alignment needs a SeqRecord target to name its reference, not {kind}!"
        )
    if not target.id or target.id == _UNKNOWN_ID:
        raise ValueError(
            f"An alignment's SeqRecord target needs an id to name its reference, not {target.id!r}!"
        )
    return target.id


def _strand(alignment: Alignment) -> str:
    """Return the strand of an alignment: minus if its target and query run opposite ways."""
    target, query = alignment.coordinates[0], alignment.coordinates[1]
    return "-" if (target[0] > target[-1]) != (query[0] > query[-1]) else "+"


def as_interval(alignment: Alignment) -> list[Bed6]:
    """Return the interval of an alignment on its target, on the strand of its query.

    Raises:
        ValueError: If the alignment's target is not a `SeqRecord`, so has no name.
    """
    target = alignment.coordinates[0]
    refname, strand = _refname(alignment), _strand(alignment)
    return [stranded(refname, int(target.min()), int(target.max()), strand)]


def as_blocks(alignment: Alignment) -> list[Bed6]:
    """Return the aligned blocks of an alignment on its target, on the strand of its query.

    Raises:
        ValueError: If the alignment's target is not a `SeqRecord`, so has no name.
    """
    refname, strand = _refname(alignment), _strand(alignment)
    ordered = sorted(
        (int(min(start, end)), int(max(start, end))) for start, end in alignment.aligned[0]
    )
    return [stranded(refname, start, end, strand) for start, end in joined_blocks(ordered)]
