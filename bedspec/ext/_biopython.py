try:
    from Bio.Align import Alignment
except ImportError as error:
    raise ImportError(
        "bedspec.ext.biopython needs Biopython: pip install 'bedspec[biopython]'"
    ) from error

from bedspec import Bed6
from bedspec.ext._blocks import joined_blocks
from bedspec.ext._stranded import stranded


def _strand(alignment: Alignment) -> str:
    """Return the strand of an alignment's query on its target."""
    query = alignment.coordinates[1]
    return "-" if query[0] > query[-1] else "+"


def as_interval(alignment: Alignment) -> list[Bed6]:
    """Return the interval of an alignment on its target, on the strand of its query."""
    target = alignment.coordinates[0]
    refname, strand = alignment.target.id, _strand(alignment)
    return [stranded(refname, int(target.min()), int(target.max()), strand)]


def as_blocks(alignment: Alignment) -> list[Bed6]:
    """Return the aligned blocks of an alignment on its target, on the strand of its query."""
    refname, strand = alignment.target.id, _strand(alignment)
    blocks = joined_blocks((int(start), int(end)) for start, end in alignment.aligned[0])
    return [stranded(refname, start, end, strand) for start, end in blocks]
