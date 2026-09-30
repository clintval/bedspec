from dataclasses import replace

from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec.overlap import TreeDetector

PLUS = Bed6(refname="chr1", start=1, end=9, name="plus", score=None, strand=BedStrand.Positive)
MINUS = Bed6(refname="chr1", start=1, end=9, name="minus", score=None, strand=BedStrand.Negative)
UNSTRANDED = Bed6(refname="chr1", start=1, end=9, name="none", score=None, strand=None)
BED3 = Bed3(refname="chr1", start=1, end=9)
QUERY = Bed6(refname="chr1", start=2, end=3, name=None, score=None, strand=BedStrand.Positive)


def detector() -> TreeDetector[Bed3 | Bed6]:
    """Build a detector holding a feature on each strand, one without a strand, and a BED3."""
    return TreeDetector([PLUS, MINUS, UNSTRANDED, BED3])


def test_strands_are_ignored_by_default() -> None:
    """Test that the strand of a feature is ignored unless asked for."""
    assert set(detector().overlapping(QUERY)) == {PLUS, MINUS, UNSTRANDED, BED3}


def test_stranded_queries_find_only_features_on_the_same_strand() -> None:
    """Test that only features on the query's strand are found when asked for."""
    assert set(detector().overlapping(QUERY, stranded=True)) == {PLUS}
    assert detector().overlaps(QUERY, stranded=True)
    assert set(detector().enclosing(QUERY, stranded=True)) == {PLUS}
    assert set(detector().enclosed_by(replace(QUERY, start=1, end=9), stranded=True)) == {PLUS}


def test_the_opposite_strand_is_queried_with_a_flipped_query() -> None:
    """Test that features on the opposite strand are found by flipping the query's strand."""
    flipped = replace(QUERY, strand=BedStrand.Positive.opposite())
    assert set(detector().overlapping(flipped, stranded=True)) == {MINUS}


def test_features_without_a_strand_never_match_a_stranded_query() -> None:
    """Test that a query or feature without a strand never matches a stranded query."""
    assert set(detector().overlapping(BED3, stranded=True)) == set()
    assert not detector().overlaps(UNSTRANDED, stranded=True)
