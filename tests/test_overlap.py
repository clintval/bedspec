import re
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from dataclasses import replace
from random import Random
from typing import Any
from unittest.mock import MagicMock

import pytest
from typing_extensions import override

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed6
from bedspec import BedLike
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import ReferenceSpan
from bedspec.overlap import TreeDetector

PAIR = BedPE(
    refname1="chr1",
    start1=10,
    end1=20,
    refname2="chr1",
    start2=50,
    end2=60,
    name=None,
    score=None,
    strand1=None,
    strand2=None,
)


@dataclass(frozen=True)
class Blocked(Bed3):
    """A record whose territory is a block of two bases at each end of its span."""

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield the first two and the last two bases of this record."""
        yield Bed3(refname=self.refname, start=self.start, end=self.start + 2)
        yield Bed3(refname=self.refname, start=self.end - 2, end=self.end)


BLOCKED = Blocked(refname="chr1", start=10, end=30)


@dataclass(frozen=True)
class Flaky(Bed3):
    """A record whose territory raises after its first span."""

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield this record, then raise."""
        yield self
        raise ValueError("a span of this record is invalid")


@dataclass(frozen=True)
class Hollow(Bed3):
    """A record whose territory has no spans."""

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield no spans."""
        yield from ()


@dataclass(frozen=True)
class Multi(BedLike):
    """A record of any spans on any references, with a strand of its own."""

    parts: tuple[Bed3 | Bed6, ...]
    strand: BedStrand | None

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield the parts of this record."""
        yield from self.parts


@dataclass
class Span:
    """A feature on a reference sequence that is not a BED record."""

    refname: str
    start: int
    end: int


@dataclass
class Region:
    """A feature on a reference sequence with a text field named territory, not a BED record."""

    refname: str
    start: int
    end: int
    territory: str = "EMEA"


def test_overlap_detector_as_iterable() -> None:
    """Test we can iterate over all the intervals we put into the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed3(refname="chr2", start=4, end=5)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2])
    assert set(detector) == {bed1, bed2}


def test_we_can_mix_types_in_the_overlap_detector() -> None:
    """Test mix input types when building the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: TreeDetector[Bed3 | Bed4] = TreeDetector([bed1, bed2])
    assert set(detector) == {bed1, bed2}


def test_we_can_add_a_feature_to_the_overlap_detector() -> None:
    """Test we can add a feature to the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: TreeDetector[Bed3 | Bed4] = TreeDetector()
    detector.add(bed1)
    detector.add(bed2)
    assert set(detector) == {bed1, bed2}


def test_we_can_add_all_features_to_the_overlap_detector() -> None:
    """Test we can add all features to the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: TreeDetector[Bed3 | Bed4] = TreeDetector()
    beds: list[Bed3 | Bed4] = [bed1, bed2]
    detector.add(*beds)
    assert set(detector) == {bed1, bed2}


def test_we_can_query_with_different_type_in_the_overlap_detector() -> None:
    """Test we can query with a different type in the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr1", start=1, end=2, name="Clint Valentine")
    detector: TreeDetector[Bed3] = TreeDetector([bed1])
    assert set(detector.overlapping(bed2)) == {bed1}


def test_we_can_those_enclosing_intervals() -> None:
    """Test that we can get intervals enclosing a given query feature."""
    bed1 = Bed3(refname="chr1", start=1, end=5)
    bed2 = Bed3(refname="chr1", start=3, end=9)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2])
    assert set(detector.enclosing(Bed3(refname="chr1", start=2, end=5))) == {bed1}
    assert set(detector.enclosing(Bed3(refname="chr1", start=3, end=8))) == {bed2}
    assert set(detector.enclosing(Bed3(refname="chr1", start=4, end=9))) == {bed2}
    assert set(detector.enclosing(Bed3(refname="chr1", start=3, end=9))) == {bed2}
    assert set(detector.enclosing(Bed3(refname="chr1", start=2, end=10))) == set()
    assert set(detector.enclosing(Bed3(refname="chr1", start=1, end=10))) == set()


def test_we_can_those_enclosed_by_intervals() -> None:
    """Test that we can get intervals enclosed by a given query feature."""
    bed1 = Bed3(refname="chr1", start=1, end=5)
    bed2 = Bed3(refname="chr1", start=3, end=9)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2])
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=2, end=5))) == set()
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=3, end=8))) == set()
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=4, end=9))) == set()
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=3, end=9))) == {bed2}
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=2, end=10))) == {bed2}
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=1, end=10))) == {bed1, bed2}


def test_we_can_query_for_overlapping_features() -> None:
    """Test we can query for features that overlap using the overlap detector."""
    bed1 = Bed3(refname="chr1", start=2, end=5)
    bed2 = Bed3(refname="chr1", start=4, end=10)
    bed3 = Bed3(refname="chr2", start=4, end=5)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert set(detector.overlapping(Bed3("chr1", start=0, end=1))) == set()
    assert set(detector.overlapping(Bed3("chr1", start=2, end=3))) == {bed1}
    assert set(detector.overlapping(Bed3("chr1", start=4, end=5))) == {bed1, bed2}
    assert set(detector.overlapping(Bed3("chr1", start=5, end=6))) == {bed2}
    assert set(detector.overlapping(Bed3("chr2", start=0, end=1))) == set()
    assert set(detector.overlapping(Bed3("chr2", start=4, end=5))) == {bed3}


def test_we_can_query_if_at_least_one_feature_overlaps() -> None:
    """Test we can query if at least one feature overlaps using the overlap detector."""
    bed1 = Bed3(refname="chr1", start=2, end=5)
    bed2 = Bed3(refname="chr1", start=4, end=10)
    bed3 = Bed3(refname="chr2", start=4, end=5)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert not detector.overlaps(Bed3("chr1", start=0, end=1))
    assert detector.overlaps(Bed3("chr1", start=2, end=3))
    assert detector.overlaps(Bed3("chr1", start=4, end=5))
    assert detector.overlaps(Bed3("chr1", start=5, end=6))
    assert not detector.overlaps(Bed3("chr2", start=0, end=1))
    assert detector.overlaps(Bed3("chr2", start=4, end=5))


def test_we_can_add_to_the_overlap_detector_after_and_before_queries() -> None:
    """Test we can add to the overlap detector after and before queries."""
    bed1 = Bed3(refname="chr1", start=2, end=5)
    bed2 = Bed3(refname="chr1", start=6, end=10)
    bed3 = Bed3(refname="chr2", start=4, end=5)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert not detector.overlaps(Bed3("chr1", start=5, end=6))

    detector.add(Bed3("chr1", start=5, end=6))

    assert detector.overlaps(Bed3("chr1", start=5, end=6))


def test_half_open_features_that_abut_do_not_overlap() -> None:
    """Test that half-open features which share only an endpoint do not overlap."""
    bed = Bed3(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3] = TreeDetector([bed])
    assert not detector.overlaps(Bed3(refname="chr1", start=5, end=10))
    assert not detector.overlaps(Bed3(refname="chr1", start=20, end=25))
    assert list(detector.overlapping(Bed3(refname="chr1", start=9, end=11))) == [bed]
    assert list(detector.overlapping(Bed3(refname="chr1", start=19, end=21))) == [bed]


def test_querying_an_unknown_reference_finds_nothing() -> None:
    """Test that querying a reference sequence with no features returns no overlaps."""
    detector: TreeDetector[Bed3] = TreeDetector([Bed3(refname="chr1", start=10, end=20)])
    assert list(detector.overlapping(Bed3(refname="chr2", start=10, end=20))) == []


def features_added_out_of_order() -> list[Bed3]:
    """Return features on two references, added out of order by start."""
    return [
        Bed3(refname="chr1", start=30, end=40),
        Bed3(refname="chr2", start=5, end=9),
        Bed3(refname="chr1", start=20, end=25),
        Bed3(refname="chr1", start=10, end=15),
        Bed3(refname="chr2", start=1, end=3),
    ]


def test_querying_while_iterating_visits_every_feature_once() -> None:
    """Test that a query inside a loop over the detector neither skips nor repeats features."""
    features = features_added_out_of_order()
    detector: TreeDetector[Bed3] = TreeDetector(features)

    seen: list[Bed3] = []
    for feature in detector:
        seen.append(feature)
        _ = list(detector.overlapping(feature))

    assert sorted(seen, key=lambda f: (f.refname, f.start)) == sorted(
        features, key=lambda f: (f.refname, f.start)
    )


def test_iteration_order_does_not_depend_on_queries() -> None:
    """Test that features are iterated in the order added, whether or not a query ran."""
    queried: TreeDetector[Bed3] = TreeDetector(features_added_out_of_order())
    _ = queried.overlaps(Bed3(refname="chr1", start=0, end=1))
    fresh: TreeDetector[Bed3] = TreeDetector(features_added_out_of_order())

    assert list(fresh) == features_added_out_of_order()
    assert list(queried) == features_added_out_of_order()


@pytest.mark.parametrize("refname", ["chr1", "chr3"])
def test_adding_while_iterating_raises(refname: str) -> None:
    """Test that adding a feature while iterating raises rather than skipping or repeating any."""
    detector: TreeDetector[Bed3] = TreeDetector(features_added_out_of_order())

    with pytest.raises(RuntimeError, match="changed during iteration"):
        for _ in detector:
            detector.add(Bed3(refname=refname, start=0, end=1))  # noqa: B909


def test_overlaps_finds_a_query_inside_a_long_feature_after_a_nested_one() -> None:
    """Test that a query inside a long feature, past a shorter one nested in it, overlaps."""
    long = Bed3(refname="chr1", start=0, end=100)
    nested = Bed3(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3] = TreeDetector([long, nested])
    query = Bed3(refname="chr1", start=30, end=40)

    assert list(detector.overlapping(query)) == [long]
    assert detector.overlaps(query)


def test_overlaps_agrees_with_overlapping_on_random_features() -> None:
    """Test that overlaps is True exactly when overlapping finds a feature, on random features."""
    rng = Random(42)
    for _ in range(200):
        starts = [rng.randrange(1000) for _ in range(rng.randrange(40))]
        features = [
            Bed3(refname="chr1", start=start, end=start + rng.choice((1, 10, 100, 500)))
            for start in starts
        ]
        detector: TreeDetector[Bed3] = TreeDetector(features)
        for _ in range(50):
            start = rng.randrange(1100)
            query = Bed3(refname="chr1", start=start, end=start + rng.randrange(50))
            assert detector.overlaps(query) is any(True for _ in detector.overlapping(query))


def random_pair(rng: Random) -> BedPE:
    """Return a pair with ends of random references, positions, lengths, and strands."""
    refnames = ("chr1", "chr2")
    strands = (BedStrand.Positive, BedStrand.Negative, None)
    start1, start2 = rng.randrange(1000), rng.randrange(1000)
    return BedPE(
        refname1=rng.choice(refnames),
        start1=start1,
        end1=start1 + rng.choice((0, 1, 10, 100)),
        refname2=rng.choice(refnames),
        start2=start2,
        end2=start2 + rng.choice((0, 1, 10, 100)),
        name=None,
        score=None,
        strand1=rng.choice(strands),
        strand2=rng.choice(strands),
    )


def test_overlaps_agrees_with_overlapping_on_random_pairs() -> None:
    """Test that overlaps agrees with overlapping, which yields each pair once, on random pairs."""
    rng = Random(42)
    for _ in range(200):
        pairs = [random_pair(rng) for _ in range(rng.randrange(20))]
        detector: TreeDetector[BedPE] = TreeDetector(pairs)
        for _ in range(50):
            start = rng.randrange(1100)
            query = Bed6(
                refname=rng.choice(("chr1", "chr2")),
                start=start,
                end=start + rng.randrange(50),
                name=None,
                score=None,
                strand=rng.choice((BedStrand.Positive, BedStrand.Negative, None)),
            )
            for stranded in (False, True):
                hits = list(detector.overlapping(query, stranded=stranded))
                assert detector.overlaps(query, stranded=stranded) is bool(hits)
                assert len(hits) == len({id(hit) for hit in hits})


def test_a_feature_that_is_not_a_bed_record_is_its_own_span() -> None:
    """Test that a feature which is not a BED record is found by its own start and end."""
    span = Span(refname="chr1", start=10, end=20)
    detector: TreeDetector[Span] = TreeDetector([span])

    assert list(detector) == [span]
    assert list(detector.overlapping(Bed3(refname="chr1", start=15, end=16))) == [span]
    assert list(detector.enclosing(Bed3(refname="chr1", start=15, end=16))) == [span]
    assert not detector.overlaps(Bed3(refname="chr1", start=20, end=21))


def test_a_point_is_found_by_its_single_base() -> None:
    """Test that a point feature is found only by a query holding or flanking its single base."""
    point = Bed2(refname="chr1", start=5)
    detector: TreeDetector[Bed2] = TreeDetector([point])

    assert list(detector.overlapping(Bed3(refname="chr1", start=5, end=6))) == [point]
    assert list(detector.overlapping(Bed3(refname="chr1", start=0, end=10))) == [point]
    assert list(detector.overlapping(Bed3(refname="chr1", start=5, end=5))) == [point]
    assert not detector.overlaps(Bed3(refname="chr1", start=4, end=5))
    assert not detector.overlaps(Bed3(refname="chr1", start=6, end=7))
    assert list(detector.enclosing(Bed3(refname="chr1", start=5, end=6))) == [point]
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=0, end=10))) == [point]


def test_a_pair_on_one_reference_is_found_by_either_end() -> None:
    """Test that a pair with both ends on one reference is found by a query on either end."""
    detector: TreeDetector[BedPE] = TreeDetector([PAIR])

    assert list(detector.overlapping(Bed3(refname="chr1", start=12, end=13))) == [PAIR]
    assert list(detector.overlapping(Bed3(refname="chr1", start=55, end=56))) == [PAIR]
    assert not detector.overlaps(Bed3(refname="chr1", start=20, end=50))
    assert not detector.overlaps(Bed3(refname="chr2", start=12, end=13))


def test_a_pair_across_two_references_is_found_by_either_end() -> None:
    """Test that a pair with ends on two references is found by a query on either end."""
    pair = replace(PAIR, refname2="chr2")
    detector: TreeDetector[BedPE] = TreeDetector([pair])

    assert list(detector.overlapping(Bed3(refname="chr1", start=12, end=13))) == [pair]
    assert list(detector.overlapping(Bed3(refname="chr2", start=55, end=56))) == [pair]
    assert not detector.overlaps(Bed3(refname="chr1", start=55, end=56))
    assert not detector.overlaps(Bed3(refname="chr2", start=12, end=13))


def test_a_pair_is_yielded_once_by_a_query_over_both_ends() -> None:
    """Test that a query over both ends of a pair yields the pair only once."""
    detector: TreeDetector[BedPE] = TreeDetector([PAIR])
    query = Bed3(refname="chr1", start=0, end=100)

    assert detector.overlaps(query)
    assert list(detector.overlapping(query)) == [PAIR]
    assert list(detector.enclosed_by(query)) == [PAIR]


def test_a_query_between_the_blocks_of_a_feature_finds_nothing() -> None:
    """Test that a query in the gap between the blocks of a feature's territory finds nothing."""
    detector: TreeDetector[Blocked] = TreeDetector([BLOCKED])
    gap = Bed3(refname="chr1", start=12, end=28)

    assert not detector.overlaps(gap)
    assert list(detector.overlapping(gap)) == []
    assert list(detector.enclosing(Bed3(refname="chr1", start=15, end=20))) == []


def test_a_query_in_either_block_of_a_feature_finds_it_once() -> None:
    """Test that a query in either or both blocks of a feature's territory finds it once."""
    detector: TreeDetector[Blocked] = TreeDetector([BLOCKED])

    for start, end in ((10, 11), (29, 30), (11, 29), (0, 100)):
        query = Bed3(refname="chr1", start=start, end=end)
        assert list(detector.overlapping(query)) == [BLOCKED]


def test_enclosing_compares_the_query_with_each_span() -> None:
    """Test that a feature encloses a query only when one of its spans encloses the query."""
    detector: TreeDetector[Blocked] = TreeDetector([BLOCKED])

    assert list(detector.enclosing(Bed3(refname="chr1", start=10, end=12))) == [BLOCKED]
    assert list(detector.enclosing(Bed3(refname="chr1", start=28, end=30))) == [BLOCKED]
    assert list(detector.enclosing(Bed3(refname="chr1", start=10, end=13))) == []
    assert list(detector.enclosing(Bed3(refname="chr1", start=11, end=29))) == []


def test_enclosed_by_requires_every_block_inside_the_query() -> None:
    """Test that a feature is enclosed by a query only when all of its blocks are inside it."""
    detector: TreeDetector[Blocked] = TreeDetector([BLOCKED])

    assert list(detector.enclosed_by(Bed3(refname="chr1", start=10, end=30))) == [BLOCKED]
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=0, end=100))) == [BLOCKED]
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=9, end=12))) == []
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=28, end=31))) == []
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=11, end=29))) == []


def test_a_pair_is_enclosed_only_when_both_ends_are_inside_the_query() -> None:
    """Test that a pair on one reference is enclosed by a query only when both ends are inside."""
    detector: TreeDetector[BedPE] = TreeDetector([PAIR])

    assert list(detector.enclosed_by(Bed3(refname="chr1", start=10, end=60))) == [PAIR]
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=0, end=30))) == []
    assert list(detector.enclosed_by(Bed3(refname="chr1", start=40, end=100))) == []


def test_a_pair_across_two_references_is_never_enclosed() -> None:
    """Test that a pair with ends on two references is not enclosed by a query on either one."""
    detector: TreeDetector[BedPE] = TreeDetector([replace(PAIR, refname2="chr2")])

    assert list(detector.enclosed_by(Bed3(refname="chr1", start=0, end=100))) == []
    assert list(detector.enclosed_by(Bed3(refname="chr2", start=0, end=100))) == []


def test_the_same_feature_added_twice_is_found_twice() -> None:
    """Test that a feature added twice is two features, found and iterated twice."""
    bed = Bed3(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3 | BedPE] = TreeDetector([bed, PAIR, bed, PAIR])
    query = Bed3(refname="chr1", start=0, end=100)

    assert list(TreeDetector([bed, bed]).overlapping(query)) == [bed, bed]
    assert Counter(detector.overlapping(query)) == Counter({bed: 2, PAIR: 2})
    assert list(detector) == [bed, PAIR, bed, PAIR]


def test_a_pair_added_after_a_query_is_yielded_once() -> None:
    """Test that a pair added after a query of features with one span each is yielded once."""
    bed = Bed3(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3 | BedPE] = TreeDetector([bed])
    query = Bed3(refname="chr1", start=0, end=100)

    assert list(detector.overlapping(query)) == [bed]
    detector.add(PAIR)
    assert Counter(detector.overlapping(query)) == Counter({bed: 1, PAIR: 1})


def test_iteration_yields_each_feature_once_in_the_order_added() -> None:
    """Test that iteration yields every feature once, in the order added, whatever its spans."""
    features: list[Bed2 | Bed3 | BedPE] = [
        replace(PAIR, refname2="chr2"),
        Bed3(refname="chr2", start=1, end=3),
        BLOCKED,
        Bed2(refname="chr1", start=5),
        PAIR,
    ]
    detector: TreeDetector[Bed2 | Bed3 | BedPE] = TreeDetector(features)

    assert list(detector) == features


def test_a_feature_with_a_territory_field_is_its_own_span() -> None:
    """Test that a feature which is not a BED record is its own span, whatever its fields."""
    region = Region(refname="chr1", start=10, end=20)
    detector: TreeDetector[Region] = TreeDetector([region])

    assert list(detector) == [region]
    assert list(detector.overlapping(Bed3(refname="chr1", start=15, end=16))) == [region]


def test_a_mock_feature_is_its_own_span() -> None:
    """Test that a mock with a reference, start, and end is found by them."""
    mock = MagicMock()
    mock.refname, mock.start, mock.end = "chr1", 10, 20
    detector: TreeDetector[MagicMock] = TreeDetector([mock])

    assert list(detector.overlapping(Bed3(refname="chr1", start=15, end=16))) == [mock]


def test_adding_a_feature_whose_territory_raises_adds_nothing() -> None:
    """Test that adding features adds none of them when any feature's territory raises."""
    bed = Bed3(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3] = TreeDetector([bed])
    query = Bed3(refname="chr1", start=0, end=100)

    with pytest.raises(ValueError, match="invalid"):
        detector.add(
            Bed3(refname="chr1", start=30, end=40), Flaky(refname="chr1", start=50, end=60)
        )

    assert list(detector) == [bed]
    assert list(detector.overlapping(query)) == [bed]
    assert list(detector.enclosed_by(query)) == [bed]


def test_adding_a_feature_without_spans_is_refused() -> None:
    """Test that adding a feature with no spans raises, naming it, and adds nothing."""
    hollow = Hollow(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3] = TreeDetector()

    with pytest.raises(ValueError, match=re.escape(repr(hollow))):
        detector.add(Bed3(refname="chr1", start=0, end=5), hollow)

    assert list(detector) == []
    assert not detector.overlaps(Bed3(refname="chr1", start=0, end=100))


STRANDS = (BedStrand.Positive, BedStrand.Negative, None)


def random_span(rng: Random) -> Bed3 | Bed6:
    """Return a short, possibly zero-length, span on a random strand, or on none."""
    refname, start = rng.choice(("chr1", "chr2")), rng.randrange(60)
    end = start + rng.choice((0, 0, 1, 2, 5, 15))
    if rng.random() < 0.4:
        return Bed3(refname=refname, start=start, end=end)
    strand = rng.choice(STRANDS)
    return Bed6(refname=refname, start=start, end=end, name=None, score=None, strand=strand)


def random_feature(rng: Random) -> Bed2 | Bed3 | Bed6 | BedPE | Multi:
    """Return a random span, point, pair, blocked record, or record of one to three spans."""
    kind = rng.randrange(5)
    start = rng.randrange(60)
    if kind == 0:
        return random_span(rng)
    if kind == 1:
        return Bed2(refname=rng.choice(("chr1", "chr2")), start=start)
    if kind == 2:
        return Blocked(refname=rng.choice(("chr1", "chr2")), start=start, end=start + 10)
    if kind == 3:
        return Multi(
            tuple(random_span(rng) for _ in range(rng.randrange(1, 4))), rng.choice(STRANDS)
        )
    first, second = random_span(rng), random_span(rng)
    return BedPE(
        refname1=first.refname,
        start1=first.start,
        end1=first.end,
        refname2=second.refname,
        start2=second.start,
        end2=second.end,
        name=None,
        score=None,
        strand1=rng.choice(STRANDS),
        strand2=rng.choice(STRANDS),
    )


def random_query(rng: Random) -> Bed6:
    """Return a possibly zero-length query on one of three references, on a strand or none."""
    refname, start = rng.choice(("chr1", "chr2", "chr3")), rng.randrange(70)
    end = start + rng.choice((0, 0, 1, 3, 10, 30, 80))
    strand = rng.choice(STRANDS)
    return Bed6(refname=refname, start=start, end=end, name=None, score=None, strand=strand)


def closed(span: ReferenceSpan) -> tuple[int, int]:
    """Return the closed interval of bases a span covers, or flanks if it is zero-length."""
    if span.start == span.end:
        return max(span.start - 1, 0), span.start
    return span.start, span.end - 1


def span_strand(span: ReferenceSpan, feature: Any) -> BedStrand | None:
    """Return the strand of a span, or of its feature if the span has none."""
    strand: BedStrand | None = getattr(span, "strand", None) or getattr(feature, "strand", None)
    return strand


def expected(features: list[Any], query: Bed6, stranded: bool) -> dict[str, list[Any]]:
    """Return the features each query method should find, by checking every span of each."""
    found: dict[str, list[Any]] = {"overlapping": [], "enclosing": [], "enclosed_by": []}
    start, end = closed(query)
    for feature in features:
        every = list(feature.territory()) if isinstance(feature, BedLike) else [feature]
        spans = [
            span
            for span in every
            if span.refname == query.refname
            and (
                not stranded
                or query.strand is not None
                and span_strand(span, feature) is query.strand
            )
        ]
        if any(closed(span)[0] <= end and start <= closed(span)[1] for span in spans):
            found["overlapping"].append(feature)
        if any(span.start <= query.start and query.end <= span.end for span in spans):
            found["enclosing"].append(feature)
        if len(spans) == len(every) and all(
            query.start <= span.start and span.end <= query.end for span in spans
        ):
            found["enclosed_by"].append(feature)
    return found


def test_queries_agree_with_checking_every_span_of_random_features() -> None:
    """Test that every query agrees with checking every span of random features, in batches."""
    rng = Random(42)
    for _ in range(300):
        detector: TreeDetector[Any] = TreeDetector()
        added: list[Any] = []
        for _ in range(rng.randrange(1, 4)):
            batch = [random_feature(rng) for _ in range(rng.randrange(8))]
            if added and rng.random() < 0.2:
                batch.append(rng.choice(added))
            detector.add(*batch)
            added.extend(batch)
            assert list(map(id, detector)) == list(map(id, added))
            for _ in range(10):
                query = random_query(rng)
                for stranded in (False, True):
                    found = expected(added, query, stranded)
                    for method, features in found.items():
                        result = getattr(detector, method)(query, stranded=stranded)
                        assert Counter(map(id, result)) == Counter(map(id, features))
                    assert detector.overlaps(query, stranded=stranded) is bool(found["overlapping"])
