from dataclasses import dataclass
from dataclasses import replace

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedPE
from bedspec import BedStrand
from bedspec.overlap import TreeDetector

PLUS = Bed6(refname="chr1", start=1, end=9, name="plus", score=None, strand=BedStrand.Positive)
MINUS = Bed6(refname="chr1", start=1, end=9, name="minus", score=None, strand=BedStrand.Negative)
UNSTRANDED = Bed6(refname="chr1", start=1, end=9, name="none", score=None, strand=None)
BED3 = Bed3(refname="chr1", start=1, end=9)
QUERY = Bed6(refname="chr1", start=2, end=3, name=None, score=None, strand=BedStrand.Positive)
PAIR = BedPE(
    refname1="chr1",
    start1=10,
    end1=20,
    refname2="chr1",
    start2=50,
    end2=60,
    name=None,
    score=None,
    strand1=BedStrand.Positive,
    strand2=BedStrand.Negative,
)


@dataclass(frozen=True)
class StrandedPoint(Bed2):
    """A point on a strand, whose 1-length span has no strand of its own."""

    strand: BedStrand | None


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


def test_stranded_queries_compare_the_strand_of_the_pair_end_they_overlap() -> None:
    """Test that a stranded query matches a pair only by the strand of an end it overlaps."""
    detector: TreeDetector[BedPE] = TreeDetector([PAIR])
    first = replace(QUERY, start=12, end=13)
    second = replace(QUERY, start=55, end=56, strand=BedStrand.Negative)
    both = replace(QUERY, start=0, end=100)

    assert list(detector.overlapping(first, stranded=True)) == [PAIR]
    assert list(detector.overlapping(second, stranded=True)) == [PAIR]
    assert list(detector.overlapping(both, stranded=True)) == [PAIR]
    assert not detector.overlaps(replace(first, strand=BedStrand.Negative), stranded=True)
    assert not detector.overlaps(replace(second, strand=BedStrand.Positive), stranded=True)
    assert detector.overlaps(replace(both, strand=BedStrand.Negative), stranded=True)


def test_stranded_enclosing_compares_the_strand_of_the_pair_end_that_encloses() -> None:
    """Test that stranded enclosing uses the strand of the pair end that encloses the query."""
    detector: TreeDetector[BedPE] = TreeDetector([PAIR])
    inside = replace(QUERY, start=12, end=13)

    assert list(detector.enclosing(inside, stranded=True)) == [PAIR]
    assert not list(detector.enclosing(replace(inside, strand=BedStrand.Negative), stranded=True))


def test_a_pair_with_ends_on_different_strands_is_not_enclosed_by_a_stranded_query() -> None:
    """Test that a stranded query encloses a pair only when both ends are on its strand."""
    around = replace(QUERY, start=0, end=100)
    same = replace(PAIR, strand2=BedStrand.Positive)
    detector: TreeDetector[BedPE] = TreeDetector([PAIR, same])

    assert set(detector.enclosed_by(around)) == {PAIR, same}
    assert list(detector.enclosed_by(around, stranded=True)) == [same]
    assert not list(detector.enclosed_by(replace(around, strand=BedStrand.Negative), stranded=True))


def test_a_span_without_a_strand_takes_the_strand_of_its_feature() -> None:
    """Test that a stranded query matches a feature by its own strand when its span has none."""
    point = StrandedPoint(refname="chr1", start=2, strand=BedStrand.Positive)
    detector: TreeDetector[StrandedPoint] = TreeDetector([point])

    assert list(detector.overlapping(QUERY, stranded=True)) == [point]
    assert not detector.overlaps(replace(QUERY, strand=BedStrand.Negative), stranded=True)
