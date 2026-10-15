"""Functions reading pysam's records as spans, to register with `register_for_bedspec`."""

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
from bedspec.ext._stranded import stranded

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def read_blocks(read: pysam.AlignedSegment) -> list[Bed6]:
    """Return the aligned blocks of a read, split at introns and deletions, on its strand."""
    if read.is_unmapped or read.reference_name is None:
        return []
    strand = "-" if read.is_reverse else "+"
    return [stranded(read.reference_name, start, end, strand) for start, end in read.get_blocks()]


def read_span(read: pysam.AlignedSegment) -> list[Bed6]:
    """Return the whole alignment of a read, introns and deletions included, on its strand."""
    if read.is_unmapped or read.reference_name is None or read.reference_end is None:
        return []
    strand = "-" if read.is_reverse else "+"
    return [stranded(read.reference_name, read.reference_start, read.reference_end, strand)]


def pair_blocks(read: pysam.AlignedSegment) -> list[Bed3]:
    """Return the aligned blocks of a read and of its mate, from the MC tag, without strand.

    Raises:
        ValueError: If the read is paired, its mate is mapped, and it has no MC tag.
    """
    spans = [Bed3(span.refname, start=span.start, end=span.end) for span in read_blocks(read)]
    if not read.is_paired or read.mate_is_unmapped or read.next_reference_name is None:
        return spans
    if not read.has_tag("MC"):
        raise ValueError(f"{read.query_name} has no MC tag: run `samtools fixmate -m` first!")
    refname, start = read.next_reference_name, read.next_reference_start
    mate_cigar = str(read.get_tag("MC"))  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
    for length, operator in _CIGAR.findall(mate_cigar):
        if operator in "M=X":
            spans.append(Bed3(refname, start=start, end=start + int(length)))
        if operator in "M=XDN":
            start += int(length)
    return spans


def variant_span(variant: pysam.VariantRecord) -> list[Bed3]:
    """Return the reference bases a variant spans."""
    return [Bed3(variant.contig, start=variant.start, end=variant.stop)]


def tabix_span(record: GTFProxy | GFF3Proxy | BedProxy) -> list[Bed6]:
    """Return the span of a GTF, GFF3 or BED line read through tabix, on its strand if any."""
    try:
        strand: object = record.strand
    except KeyError:
        strand = None
    return [stranded(str(record.contig), int(record.start), int(record.end), strand)]
