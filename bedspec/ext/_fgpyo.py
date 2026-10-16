from collections.abc import Iterator
from itertools import chain

try:
    from fgpyo.sam import SupplementaryAlignment
    from fgpyo.sam import Template
    from pysam import AlignedSegment
except ImportError as error:
    raise ImportError("bedspec.ext.fgpyo needs fgpyo: pip install 'bedspec[fgpyo]'") from error

from bedspec import Bed6
from bedspec import ReferenceSpan
from bedspec.ext import _pysam
from bedspec.ext._stranded import stranded


def as_interval(alignment: SupplementaryAlignment) -> list[Bed6]:
    """Return the interval of the alignment an SA tag entry describes, on its strand."""
    strand = "+" if alignment.is_forward else "-"
    return [stranded(alignment.reference_name, alignment.start, alignment.end, strand)]


def _records(template: Template) -> Iterator[AlignedSegment]:
    """Yield the primary records of a template, then its supplementary records."""
    return chain(template.primary_recs(), template.r1_supplementals, template.r2_supplementals)


def as_template_intervals(template: Template) -> list[ReferenceSpan]:
    """Return the interval of each primary and supplementary record of a template, on its strand."""
    return [span for record in _records(template) for span in _pysam.as_interval(record)]


def as_template_blocks(template: Template) -> list[Bed6]:
    """Return the aligned blocks of each primary and supplementary record of a template."""
    return [span for record in _records(template) for span in _pysam.as_blocks(record)]
