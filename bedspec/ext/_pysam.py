import re

try:
    import pysam
    from pysam.libctabixproxies import BedProxy
    from pysam.libctabixproxies import GFF3Proxy
    from pysam.libctabixproxies import GTFProxy
except ImportError as error:
    raise ImportError("bedspec.ext.pysam needs pysam: pip install 'bedspec[pysam]'") from error

from bedspec import Bed3
from bedspec import Bed6
from bedspec import ReferenceSpan
from bedspec.ext._blocks import joined_blocks
from bedspec.ext._stranded import stranded

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_ALIGNED = "M=X"
_ON_REFERENCE = "M=XDN"
_ON_REFERENCE_CODES = frozenset({0, 2, 3, 7, 8})


def _consumes_reference(read: pysam.AlignedSegment) -> bool:
    """Return whether a read's CIGAR consumes a reference base, which htslib assumes one does."""
    return any(code in _ON_REFERENCE_CODES for code, _ in read.cigartuples or ())


def as_interval(
    record: pysam.AlignedSegment | pysam.VariantRecord | GTFProxy | GFF3Proxy | BedProxy,
) -> list[ReferenceSpan]:
    """Return a record's interval, the one span it is found by.

    A read's interval is its whole alignment, introns and deletions included, on its strand, and
    an unmapped read has none; a variant's is its reference allele; and a GTF, GFF3 or BED line
    read through tabix has its own, on its strand if it has one.
    """
    if isinstance(record, pysam.AlignedSegment):
        refname, end = record.reference_name, record.reference_end
        if record.is_unmapped or refname is None or end is None:
            return []
        if end - record.reference_start == 1 and not _consumes_reference(record):
            return []
        return [stranded(refname, record.reference_start, end, "-" if record.is_reverse else "+")]
    if isinstance(record, pysam.VariantRecord):
        return [Bed3(record.contig, start=record.start, end=record.stop)]
    try:
        strand: object = record.strand
    except KeyError:
        strand = None
    return [stranded(str(record.contig), int(record.start), int(record.end), strand)]


def as_blocks(read: pysam.AlignedSegment) -> list[Bed6]:
    """Return the aligned blocks of a read, split at introns and deletions, on its strand."""
    if read.is_unmapped or read.reference_name is None:
        return []
    strand = "-" if read.is_reverse else "+"
    blocks = joined_blocks(read.get_blocks())
    return [stranded(read.reference_name, start, end, strand) for start, end in blocks]


def _mate(read: pysam.AlignedSegment) -> tuple[str, int, list[tuple[int, str]]] | None:
    """Return the reference, start and CIGAR operations of a read's mate, or None if it has none.

    Raises:
        ValueError: If the read is paired, its mate is mapped, and it has no MC tag.
    """
    refname = read.next_reference_name
    if not read.is_paired or read.mate_is_unmapped or refname is None:
        return None
    if read.next_reference_start < 0:
        return None
    if not read.has_tag("MC"):
        raise ValueError(f"{read.query_name} has no MC tag: run `samtools fixmate -m` first!")
    mate_cigar = str(read.get_tag("MC"))  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
    if mate_cigar == "*":
        return None
    operations = [(int(length), operator) for length, operator in _CIGAR.findall(mate_cigar)]
    return refname, read.next_reference_start, operations


def as_paired_intervals(read: pysam.AlignedSegment) -> list[Bed3]:
    """Return the interval of a read and of its mate, from the MC tag, never joined, unstranded.

    Raises:
        ValueError: If the read is paired, its mate is mapped, and it has no MC tag.
    """
    intervals = [Bed3(span.refname, start=span.start, end=span.end) for span in as_interval(read)]
    mate = _mate(read)
    if mate is None:
        return intervals
    refname, start, operations = mate
    length = sum(length for length, operator in operations if operator in _ON_REFERENCE)
    intervals.append(Bed3(refname, start=start, end=start + length))
    return intervals


def as_paired_blocks(read: pysam.AlignedSegment) -> list[Bed3]:
    """Return the aligned blocks of a read and of its mate, from the MC tag, without strand.

    Raises:
        ValueError: If the read is paired, its mate is mapped, and it has no MC tag.
    """
    blocks = [Bed3(span.refname, start=span.start, end=span.end) for span in as_blocks(read)]
    mate = _mate(read)
    if mate is None:
        return blocks
    refname, start, operations = mate
    mate_blocks: list[tuple[int, int]] = []
    for length, operator in operations:
        if operator in _ALIGNED:
            mate_blocks.append((start, start + length))
        if operator in _ON_REFERENCE:
            start += length
    blocks.extend(Bed3(refname, start=start, end=end) for start, end in joined_blocks(mate_blocks))
    return blocks
