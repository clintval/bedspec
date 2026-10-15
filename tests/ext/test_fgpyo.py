"""Tests of the functions that read fgpyo's records as spans."""

from fgpyo.sam import SupplementaryAlignment
from fgpyo.sam import Template
from fgpyo.sam.builder import SamBuilder

import bedspec.ext.fgpyo
from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import register_for_bedspec
from bedspec.ext.fgpyo import supplementary_span
from bedspec.ext.fgpyo import template_blocks
from bedspec.overlap import TreeDetector


def test_template_blocks_are_the_blocks_of_its_mapped_records_unstranded() -> None:
    """Test that a template's blocks are those of its primary and supplementary records."""
    builder = SamBuilder()
    r1, r2 = builder.add_pair(
        chrom="chr1", start1=100, start2=900, cigar1="10M500N10M", cigar2="20M"
    )
    assert template_blocks(Template.build([r1, r2])) == [
        Bed3("chr1", start=100, end=110),
        Bed3("chr1", start=610, end=620),
        Bed3("chr1", start=900, end=920),
    ]


def test_template_blocks_skip_an_unmapped_mate() -> None:
    """Test that a template with an unmapped mate is found by its mapped one."""
    mapped, unmapped = SamBuilder().add_pair(chrom="chr1", start1=300, cigar1="20M")
    assert template_blocks(Template.build([mapped, unmapped])) == [Bed3("chr1", start=300, end=320)]


def test_supplementary_span_is_0_based_on_its_strand() -> None:
    """Test that an SA tag entry, 1-based in the tag, is read 0-based on its strand."""
    alignment = SupplementaryAlignment.parse("chr2,1001,-,30M70S,60,0")
    assert supplementary_span(alignment) == [
        Bed6("chr2", start=1000, end=1030, name=None, score=None, strand=BedStrand.Negative)
    ]


def test_registered_templates_are_kept_whole_when_either_read_is_on_target() -> None:
    """Test that templates registered by their blocks are found by either read."""
    builder = SamBuilder()
    on = Template.build(builder.add_pair(chrom="chr1", start1=100, start2=900))
    off = Template.build(builder.add_pair(chrom="chr1", start1=100, start2=300))
    targets = TreeDetector([Bed3("chr1", start=950, end=960)])
    with register_for_bedspec(Template, template_blocks):
        assert [template for template in (on, off) if targets.overlaps(template)] == [on]


def test_bedspec_ext_fgpyo_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.fgpyo holds its functions and none of the names it imports."""
    module = bedspec.ext.fgpyo
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert public == sorted(module.__all__) == ["supplementary_span", "template_blocks"]
