"""Functions reading Biopython's records as spans, to register with `register_for_bedspec`."""

try:
    from Bio.Align import Alignment
except ImportError as error:
    raise ImportError(
        "bedspec.ext.biopython needs Biopython: pip install 'bedspec[biopython]'"
    ) from error

from bedspec import Bed6
from bedspec.ext._stranded import stranded


def alignment_blocks(alignment: Alignment) -> list[Bed6]:
    """Return the aligned blocks of an alignment on its target, on the strand of its query."""
    query = alignment.coordinates[1]
    strand = "-" if query[0] > query[-1] else "+"
    refname = alignment.target.id
    return [stranded(refname, int(start), int(end), strand) for start, end in alignment.aligned[0]]
