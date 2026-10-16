"""Tests of the functions that read fgpyo's records as spans."""

from fgpyo.sam import SupplementaryAlignment
from fgpyo.sam import Template
from fgpyo.sam.builder import SamBuilder

import bedspec.ext.fgpyo
from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import register_locatable
from bedspec.ext.fgpyo import as_interval
from bedspec.ext.fgpyo import as_template_blocks
from bedspec.ext.fgpyo import as_template_intervals
from bedspec.overlap import TreeDetector


def stranded(start: int, end: int, strand: BedStrand) -> Bed6:
    """Return a span on chr1 on a strand."""
    return Bed6("chr1", start=start, end=end, name=None, score=None, strand=strand)


def test_as_template_blocks_are_the_blocks_of_its_mapped_records_on_their_strands() -> None:
    """Test that a template's blocks are those of its records, each on its record's strand."""
    builder = SamBuilder()
    r1, r2 = builder.add_pair(
        chrom="chr1", start1=100, start2=900, cigar1="10M500N10M", cigar2="20M"
    )
    assert as_template_blocks(Template.build([r1, r2])) == [
        stranded(100, 110, BedStrand.Positive),
        stranded(610, 620, BedStrand.Positive),
        stranded(900, 920, BedStrand.Negative),
    ]


def test_the_paired_spans_of_a_template_skip_an_unmapped_mate() -> None:
    """Test that a template with an unmapped mate is found by its mapped one."""
    mapped, unmapped = SamBuilder().add_pair(chrom="chr1", start1=300, cigar1="20M")
    template = Template.build([mapped, unmapped])
    assert as_template_blocks(template) == [stranded(300, 320, BedStrand.Positive)]
    assert as_template_intervals(template) == [stranded(300, 320, BedStrand.Positive)]


def test_as_template_intervals_are_the_intervals_of_its_mapped_records_on_their_strands() -> None:
    """Test that a template's paired intervals are its records', never joined, on their strands."""
    builder = SamBuilder()
    r1, r2 = builder.add_pair(
        chrom="chr1", start1=100, start2=900, cigar1="10M500N10M", cigar2="20M"
    )
    assert as_template_intervals(Template.build([r1, r2])) == [
        stranded(100, 620, BedStrand.Positive),
        stranded(900, 920, BedStrand.Negative),
    ]


def test_as_interval_of_a_supplementary_alignment_is_0_based_on_its_strand() -> None:
    """Test that an SA tag entry, 1-based in the tag, is read 0-based on its strand."""
    alignment = SupplementaryAlignment.parse("chr2,1001,-,30M70S,60,0")
    assert as_interval(alignment) == [
        Bed6("chr2", start=1000, end=1030, name=None, score=None, strand=BedStrand.Negative)
    ]


def test_registered_templates_are_kept_whole_when_either_read_is_on_target() -> None:
    """Test that templates registered by their blocks are found by either read."""
    builder = SamBuilder()
    on = Template.build(builder.add_pair(chrom="chr1", start1=100, start2=900))
    off = Template.build(builder.add_pair(chrom="chr1", start1=100, start2=300))
    targets = TreeDetector([Bed3("chr1", start=950, end=960)])
    with register_locatable(Template, as_template_blocks):
        assert [template for template in (on, off) if targets.overlaps(template)] == [on]


def test_bedspec_ext_fgpyo_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.fgpyo holds its functions and none of the names it imports."""
    module = bedspec.ext.fgpyo
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert (
        public
        == sorted(module.__all__)
        == ["as_interval", "as_template_blocks", "as_template_intervals"]
    )
