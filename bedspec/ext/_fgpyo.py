from itertools import chain

try:
    from fgpyo.sam import SupplementaryAlignment
    from fgpyo.sam import Template
except ImportError as error:
    raise ImportError("bedspec.ext.fgpyo needs fgpyo: pip install 'bedspec[fgpyo]'") from error

from bedspec import Bed3
from bedspec import Bed6
from bedspec.ext._pysam import read_blocks
from bedspec.ext._stranded import stranded


def template_blocks(template: Template) -> list[Bed3]:
    """Return the aligned blocks of a template's primary and supplementary records, unstranded."""
    records = chain(template.primary_recs(), template.r1_supplementals, template.r2_supplementals)
    return [
        Bed3(span.refname, start=span.start, end=span.end)
        for record in records
        for span in read_blocks(record)
    ]


def supplementary_span(alignment: SupplementaryAlignment) -> list[Bed6]:
    """Return the span of a supplementary alignment from an SA tag, on its strand."""
    strand = "+" if alignment.is_forward else "-"
    return [stranded(alignment.reference_name, alignment.start, alignment.end, strand)]
