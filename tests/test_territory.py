from collections.abc import Iterator
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from random import Random
from typing import Any

import pytest
from typing_extensions import override

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedPE
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import ReferenceSpan
from bedspec import Territory


@dataclass(frozen=True)
class Blocked(Bed3):
    """A record whose spans are a block of two bases at each end of it."""

    @override
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield the first two and the last two bases of this record."""
        yield Bed3(refname=self.refname, start=self.start, end=self.start + 2)
        yield Bed3(refname=self.refname, start=self.end - 2, end=self.end)


@dataclass
class Span:
    """A span on a reference sequence that is not a BED record."""

    refname: str
    start: int
    end: int


def pair(refname1: str, start1: int, end1: int, refname2: str, start2: int, end2: int) -> BedPE:
    """Return a pair of spans without a name, score, or strands."""
    return BedPE(
        refname1=refname1,
        start1=start1,
        end1=end1,
        refname2=refname2,
        start2=start2,
        end2=end2,
        name=None,
        score=None,
        strand1=None,
        strand2=None,
    )


def chr1(*spans: tuple[int, int]) -> Territory:
    """Return the territory of spans on chr1, each given by its start and end."""
    return Territory(Bed3("chr1", start=start, end=end) for start, end in spans)


def bounds(territory: Territory) -> list[tuple[str, int, int]]:
    """Return the reference, start and end of each span of a territory, in order."""
    return [(span.refname, span.start, span.end) for span in territory]


def test_overlapping_spans_are_joined() -> None:
    """Test that spans that overlap are joined into one."""
    assert bounds(chr1((0, 10), (5, 15))) == [("chr1", 0, 15)]
    assert bounds(chr1((5, 15), (0, 10))) == [("chr1", 0, 15)]


def test_abutting_spans_are_joined() -> None:
    """Test that a span that starts where another ends is joined to it."""
    assert bounds(chr1((0, 5), (5, 10))) == [("chr1", 0, 10)]
    assert bounds(chr1((5, 10), (0, 5), (10, 12))) == [("chr1", 0, 12)]


def test_nested_spans_are_held_by_the_outer_span() -> None:
    """Test that spans inside another span are held by that span."""
    assert bounds(chr1((0, 20), (5, 10), (0, 20), (19, 20))) == [("chr1", 0, 20)]


def test_disjoint_spans_are_kept_apart_and_sorted_by_start() -> None:
    """Test that spans with a gap between them are kept apart, by start, even by one base."""
    assert bounds(chr1((20, 30), (0, 10), (11, 15))) == [
        ("chr1", 0, 10),
        ("chr1", 11, 15),
        ("chr1", 20, 30),
    ]


def test_spans_on_different_references_are_never_joined() -> None:
    """Test that spans on different references are kept apart, though their positions overlap."""
    territory = Territory([Bed3("chr1", start=0, end=10), Bed3("chr2", start=5, end=15)])
    assert bounds(territory) == [("chr1", 0, 10), ("chr2", 5, 15)]


def test_references_are_in_the_order_they_were_first_added() -> None:
    """Test that references are iterated in the order they were first added, not by name."""
    territory = Territory([
        Bed3("chr2", start=50, end=60),
        Bed3("chr10", start=5, end=6),
        Bed3("chr1", start=0, end=1),
        Bed3("chr2", start=0, end=10),
    ])
    assert bounds(territory) == [
        ("chr2", 0, 10),
        ("chr2", 50, 60),
        ("chr10", 5, 6),
        ("chr1", 0, 1),
    ]


def test_spans_are_yielded_as_bed3_records() -> None:
    """Test that a territory yields BED3 records, whatever kind of records it was built from."""
    territory = Territory([Bed6("chr1", start=0, end=10, name="a", score=1, strand=None)])
    assert [type(span) for span in territory] == [Bed3]
    assert list(territory) == [Bed3("chr1", start=0, end=10)]


def test_an_empty_territory_holds_no_bases() -> None:
    """Test that a territory built from no spans holds no bases."""
    territory = Territory()
    assert list(territory) == []
    assert territory.length == 0
    assert not territory
    assert territory == Territory([])


def test_a_zero_length_span_holds_no_bases() -> None:
    """Test that a zero-length span, such as an insertion, adds no bases and no reference."""
    territory = Territory([Bed3("chr1", start=5, end=5), Bed3("chr2", start=0, end=0)])
    assert territory == Territory()
    assert territory.length == 0
    assert not territory


def test_a_zero_length_span_neither_joins_nor_splits_spans() -> None:
    """Test that a zero-length span between or inside spans leaves them as they were."""
    assert bounds(chr1((0, 5), (5, 5), (6, 10))) == [("chr1", 0, 5), ("chr1", 6, 10)]
    assert bounds(chr1((0, 5), (6, 6), (7, 10))) == [("chr1", 0, 5), ("chr1", 7, 10)]
    assert bounds(chr1((0, 10), (5, 5))) == [("chr1", 0, 10)]
    assert chr1((0, 10)) - chr1((5, 5)) == chr1((0, 10))
    assert chr1((0, 10)) & chr1((5, 5)) == Territory()


def test_the_intersection_of_partly_overlapping_spans() -> None:
    """Test that the intersection of partly overlapping spans holds only their shared bases."""
    assert bounds(chr1((0, 10)) & chr1((5, 15))) == [("chr1", 5, 10)]
    assert bounds(chr1((5, 15)) & chr1((0, 10))) == [("chr1", 5, 10)]


def test_the_intersection_with_a_span_inside_another() -> None:
    """Test that the intersection with a span inside another is the inner span, in either order."""
    assert bounds(chr1((0, 20)) & chr1((5, 10))) == [("chr1", 5, 10)]
    assert bounds(chr1((5, 10)) & chr1((0, 20))) == [("chr1", 5, 10)]
    assert bounds(chr1((0, 20)) & chr1((0, 20))) == [("chr1", 0, 20)]


def test_the_intersection_of_spans_that_share_no_bases_is_empty() -> None:
    """Test that disjoint, abutting, and other-reference spans intersect in no bases."""
    assert chr1((0, 5)) & chr1((6, 10)) == Territory()
    assert chr1((0, 5)) & chr1((5, 10)) == Territory()
    assert chr1((0, 5)) & Territory([Bed3("chr2", start=0, end=5)]) == Territory()
    assert chr1((0, 5)) & Territory() == Territory()


def test_the_intersection_of_many_spans_with_many() -> None:
    """Test that the intersection of several spans with several holds every shared base."""
    first = chr1((0, 10), (20, 30), (40, 50))
    second = chr1((5, 25), (28, 42), (49, 60))
    assert bounds(first & second) == [
        ("chr1", 5, 10),
        ("chr1", 20, 25),
        ("chr1", 28, 30),
        ("chr1", 40, 42),
        ("chr1", 49, 50),
    ]


def test_the_intersection_keeps_the_order_of_the_left_references() -> None:
    """Test that the intersection yields references in the order of the left territory."""
    first = Territory([Bed3("chr2", start=0, end=10), Bed3("chr1", start=0, end=10)])
    second = Territory([Bed3("chr1", start=5, end=15), Bed3("chr2", start=5, end=15)])
    assert bounds(first & second) == [("chr2", 5, 10), ("chr1", 5, 10)]
    assert bounds(second & first) == [("chr1", 5, 10), ("chr2", 5, 10)]


def test_the_difference_of_partly_overlapping_spans() -> None:
    """Test that subtracting a partly overlapping span trims the bases it shares."""
    assert bounds(chr1((0, 10)) - chr1((5, 15))) == [("chr1", 0, 5)]
    assert bounds(chr1((5, 15)) - chr1((0, 10))) == [("chr1", 10, 15)]


def test_the_difference_with_a_span_inside_another_splits_it() -> None:
    """Test that subtracting a span inside another leaves the bases either side of it."""
    assert bounds(chr1((0, 20)) - chr1((5, 10))) == [("chr1", 0, 5), ("chr1", 10, 20)]
    assert bounds(chr1((0, 20)) - chr1((0, 10))) == [("chr1", 10, 20)]
    assert bounds(chr1((0, 20)) - chr1((10, 20))) == [("chr1", 0, 10)]


def test_the_difference_with_an_enclosing_span_is_empty() -> None:
    """Test that subtracting a span that encloses another leaves none of its bases."""
    assert chr1((5, 10)) - chr1((0, 20)) == Territory()
    assert chr1((5, 10)) - chr1((5, 10)) == Territory()


def test_the_difference_with_spans_that_share_no_bases_changes_nothing() -> None:
    """Test that subtracting disjoint, abutting, other-reference, or no spans changes nothing."""
    assert chr1((0, 5)) - chr1((6, 10)) == chr1((0, 5))
    assert chr1((5, 10)) - chr1((0, 5), (10, 15)) == chr1((5, 10))
    assert chr1((0, 5)) - Territory([Bed3("chr2", start=0, end=5)]) == chr1((0, 5))
    assert chr1((0, 5)) - Territory() == chr1((0, 5))
    assert Territory() - chr1((0, 5)) == Territory()


def test_the_difference_of_many_spans_with_many() -> None:
    """Test that one span can trim several, and several spans can trim one."""
    assert bounds(chr1((0, 10), (20, 30)) - chr1((5, 25))) == [("chr1", 0, 5), ("chr1", 25, 30)]
    assert bounds(chr1((0, 30)) - chr1((2, 4), (10, 12), (28, 30))) == [
        ("chr1", 0, 2),
        ("chr1", 4, 10),
        ("chr1", 12, 28),
    ]


def test_the_union_joins_spans_of_both_territories() -> None:
    """Test that the union joins spans of both territories that overlap or abut."""
    assert bounds(chr1((0, 5), (20, 30)) | chr1((5, 10), (25, 40))) == [
        ("chr1", 0, 10),
        ("chr1", 20, 40),
    ]


def test_the_union_puts_the_left_references_first() -> None:
    """Test that the union yields the left territory's references, then the right's new ones."""
    first = Territory([Bed3("chr2", start=0, end=1), Bed3("chr1", start=0, end=1)])
    second = Territory([Bed3("chr3", start=0, end=1), Bed3("chr1", start=5, end=6)])
    assert [span.refname for span in first | second] == ["chr2", "chr1", "chr1", "chr3"]


def test_the_length_counts_each_base_once() -> None:
    """Test that the length counts each base once, however many spans hold it."""
    territory = Territory([
        Bed3("chr1", start=0, end=10),
        Bed3("chr1", start=5, end=15),
        Bed3("chr1", start=20, end=21),
        Bed3("chr2", start=0, end=3),
    ])
    assert territory.length == 15 + 1 + 3


def test_contains_finds_a_base_by_its_reference_and_position() -> None:
    """Test that contains tests one base, holding a span's start and not its end."""
    territory = chr1((10, 20), (30, 40))
    assert [position for position in range(50) if territory.contains("chr1", position)] == [
        *range(10, 20),
        *range(30, 40),
    ]
    assert not territory.contains("chr2", 15)


def test_in_tests_that_every_base_of_a_span_is_held() -> None:
    """Test that a span is in a territory when the territory holds every base of it."""
    territory = chr1((10, 20), (20, 30), (40, 50))
    assert Bed3("chr1", start=10, end=30) in territory
    assert Bed3("chr1", start=15, end=25) in territory
    assert Bed3("chr1", start=42, end=43) in territory
    assert Bed3("chr1", start=5, end=15) not in territory
    assert Bed3("chr1", start=25, end=45) not in territory
    assert Bed3("chr1", start=0, end=5) not in territory
    assert Bed3("chr2", start=15, end=16) not in territory
    assert Span("chr1", 12, 18) in territory


def test_in_tests_that_either_base_beside_a_zero_length_span_is_held() -> None:
    """Test that a zero-length span is in a territory when either base beside it is held."""
    territory = chr1((10, 20), (21, 25))
    insertions = [Bed3("chr1", start=position, end=position) for position in range(30)]
    assert [insertion.start for insertion in insertions if insertion in territory] == [
        *range(10, 26)
    ]
    assert Bed3("chr2", start=10, end=10) not in territory


def test_in_tests_every_span_of_a_bed_record() -> None:
    """Test that a BED record is in a territory only when every span of its territory is."""
    territory = chr1((10, 20), (50, 60))
    assert Bed2("chr1", 10) in territory
    assert Bed2("chr1", 20) not in territory
    assert pair("chr1", 10, 15, "chr1", 50, 60) in territory
    assert pair("chr1", 10, 15, "chr1", 45, 55) not in territory
    assert pair("chr1", 10, 15, "chr2", 50, 60) not in territory
    assert Blocked("chr1", start=18, end=52) in territory
    assert Blocked("chr1", start=18, end=62) not in territory


def test_in_is_false_for_anything_that_is_not_a_span() -> None:
    """Test that an object without a reference, start, and end is never in a territory."""
    territory = chr1((0, 10))
    assert "chr1" not in territory
    assert ("chr1", 5) not in territory


def test_a_point_adds_its_base() -> None:
    """Test that a BED2 point adds its one base."""
    territory = Territory([Bed2("chr1", 5), Bed2("chr1", 6), Bed2("chr1", 8)])
    assert bounds(territory) == [("chr1", 5, 7), ("chr1", 8, 9)]


def test_a_pair_adds_both_of_its_ends() -> None:
    """Test that a BEDPE record adds both of its ends, on one reference or two."""
    assert bounds(Territory([pair("chr1", 10, 20, "chr1", 50, 60)])) == [
        ("chr1", 10, 20),
        ("chr1", 50, 60),
    ]
    assert bounds(Territory([pair("chr2", 50, 60, "chr1", 10, 20)])) == [
        ("chr2", 50, 60),
        ("chr1", 10, 20),
    ]
    assert bounds(Territory([pair("chr1", 10, 20, "chr1", 15, 30)])) == [("chr1", 10, 30)]


def test_a_record_of_several_spans_adds_only_its_spans() -> None:
    """Test that a record adds the spans of its territory and not the bases between them."""
    territory = Territory([Blocked("chr1", start=10, end=30)])
    assert bounds(territory) == [("chr1", 10, 12), ("chr1", 28, 30)]
    assert territory.length == 4
    assert not territory.contains("chr1", 20)
    assert bounds(territory & chr1((0, 100))) == [("chr1", 10, 12), ("chr1", 28, 30)]
    assert bounds(chr1((0, 40)) - territory) == [
        ("chr1", 0, 10),
        ("chr1", 12, 28),
        ("chr1", 30, 40),
    ]


def bed12(start: int, blocks: tuple[tuple[int, int], ...]) -> Bed12:
    """Return a BED12 record on chr1 with the given blocks, each a start offset and a size."""
    offset, size = blocks[-1]
    return Bed12(
        "chr1",
        start=start,
        end=start + offset + size,
        name=None,
        score=None,
        strand=BedStrand.Positive,
        thick_start=None,
        thick_end=None,
        item_rgb=None,
        block_count=len(blocks),
        block_sizes=tuple(size for _, size in blocks),
        block_starts=tuple(offset for offset, _ in blocks),
    )


def test_a_bed12_adds_its_blocks_and_not_its_introns() -> None:
    """Test that a BED12 record adds the bases of its blocks, the same as its own territory."""
    record = bed12(100, ((0, 10), (40, 20), (90, 10)))
    territory = Territory([record])
    assert bounds(territory) == [("chr1", 100, 110), ("chr1", 140, 160), ("chr1", 190, 200)]
    assert territory == record.territory()
    assert territory.length == 40
    assert record in territory
    assert Bed3("chr1", start=100, end=200) not in territory
    assert bounds(chr1((0, 300)) - territory) == [
        ("chr1", 0, 100),
        ("chr1", 110, 140),
        ("chr1", 160, 190),
        ("chr1", 200, 300),
    ]


def test_a_span_that_is_not_a_bed_record_adds_itself() -> None:
    """Test that any object with a reference, start, and end adds itself."""
    assert Territory([Span("chr1", 0, 10), Span("chr1", 5, 15)]) == chr1((0, 15))


def test_a_territory_can_be_built_from_a_territory() -> None:
    """Test that a territory built from another holds the same bases."""
    territory = chr1((0, 5), (10, 15))
    assert Territory(territory) == territory


@pytest.mark.parametrize(
    ("span", "message"),
    [
        (Span("chr1", 5, 4), "end must be greater than or equal to start!"),
        (Span("chr1", -1, 4), "start must be greater than or equal to 0!"),
        (Span("", 0, 4), "refname must not be empty!"),
    ],
)
def test_an_invalid_span_is_refused(span: Span, message: str) -> None:
    """Test that a span that is not a valid BED span is refused."""
    with pytest.raises(ValueError, match=message):
        _ = Territory([span])


def test_combining_territories_leaves_them_unchanged() -> None:
    """Test that combining territories returns a new territory and changes neither of them."""
    first, second = chr1((0, 10), (20, 30)), chr1((5, 25))
    for combined in (first | second, first & second, first - second):
        assert combined is not first
        assert combined is not second
    assert bounds(first) == [("chr1", 0, 10), ("chr1", 20, 30)]
    assert bounds(second) == [("chr1", 5, 25)]


def test_a_territory_cannot_be_changed() -> None:
    """Test that a territory's length cannot be set and it takes no new attributes."""
    territory = chr1((0, 10))
    with pytest.raises(AttributeError):
        territory.length = 5  # type: ignore[misc]  # pyright: ignore[reportAttributeAccessIssue]
    with pytest.raises(AttributeError):
        territory.spans = []  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]
    assert territory.length == 10


def test_territories_holding_the_same_bases_are_equal() -> None:
    """Test that territories are equal, and hash alike, when they hold the same bases."""
    first = Territory([Bed3("chr1", start=0, end=10), Bed3("chr2", start=0, end=5)])
    second = Territory([
        Bed3("chr2", start=0, end=5),
        Bed3("chr1", start=5, end=10),
        Bed3("chr1", start=0, end=5),
    ])
    assert first == second
    assert hash(first) == hash(second)
    assert len({first, second}) == 1
    assert first != Territory([Bed3("chr1", start=0, end=10)])
    assert first != chr1((0, 10), (11, 12))
    assert first != list(first)


def test_combining_with_something_other_than_a_territory_raises() -> None:
    """Test that combining a territory with anything other than a territory raises a TypeError."""
    territory = chr1((0, 10))
    spans: Any = [Bed3("chr1", start=5, end=15)]
    with pytest.raises(TypeError, match="unsupported operand"):
        _ = territory | spans
    with pytest.raises(TypeError, match="unsupported operand"):
        _ = territory & spans
    with pytest.raises(TypeError, match="unsupported operand"):
        _ = territory - spans


def test_the_representation_lists_the_spans() -> None:
    """Test that the representation of a territory lists its spans as BED3 records."""
    assert repr(chr1((0, 5), (5, 10))) == "Territory([Bed3(refname='chr1', start=0, end=10)])"
    assert repr(Territory()) == "Territory([])"


def test_a_territory_written_and_read_back_is_equal(tmp_path: Path) -> None:
    """Test that a territory written as BED3 and read back holds the same bases."""
    path = tmp_path / "territory.bed"
    territory = Territory([pair("chr2", 50, 60, "chr1", 10, 20), Bed2("chr1", 20)])
    with BedWriter.from_path[Bed3](path) as writer:
        for span in territory:
            writer.write(span)
    with BedReader.from_path[Bed3](path) as reader:
        assert Territory(reader) == territory


def random_feature(rng: Random) -> Bed2 | Bed3 | Bed12 | BedPE | Blocked:
    """Return a short, possibly zero-length, span, a point, a pair, or a record of blocks."""
    refname, start = rng.choice(("chr1", "chr2")), rng.randrange(80)
    kind = rng.randrange(5)
    if kind == 0:
        return Bed2(refname, start)
    if kind == 1:
        return Blocked(refname, start=start, end=start + rng.randrange(4, 12))
    if kind == 4:
        offsets, sizes = [0], [rng.randrange(1, 6) for _ in range(rng.randrange(1, 4))]
        for size in sizes[:-1]:
            offsets.append(offsets[-1] + size + rng.randrange(0, 8))
        return Bed12(
            refname,
            start=start,
            end=start + offsets[-1] + sizes[-1],
            name=None,
            score=None,
            strand=rng.choice((BedStrand.Positive, BedStrand.Negative, None)),
            thick_start=None,
            thick_end=None,
            item_rgb=None,
            block_count=len(sizes),
            block_sizes=tuple(sizes),
            block_starts=tuple(offsets),
        )
    end = start + rng.choice((0, 0, 1, 2, 5, 15, 30))
    if kind == 2:
        return Bed3(refname, start=start, end=end)
    other = rng.randrange(80)
    return BedPE(
        refname1=refname,
        start1=start,
        end1=end,
        refname2=rng.choice(("chr1", "chr2")),
        start2=other,
        end2=other + rng.choice((0, 1, 10)),
        name=None,
        score=None,
        strand1=rng.choice((BedStrand.Positive, BedStrand.Negative, None)),
        strand2=None,
    )


def bases_of(features: list[Bed2 | Bed3 | Bed12 | BedPE | Blocked]) -> set[tuple[str, int]]:
    """Return every base of every span of the given features."""
    return {
        (span.refname, position)
        for feature in features
        for span in feature.spans()
        for position in range(span.start, span.end)
    }


def is_held(feature: Bed2 | Bed3 | Bed12 | BedPE | Blocked, bases: set[tuple[str, int]]) -> bool:
    """Return whether every base of every span is held, or either base beside a zero-length span."""
    for span in feature.spans():
        if span.start == span.end:
            if not {(span.refname, span.start - 1), (span.refname, span.start)} & bases:
                return False
        elif not all((span.refname, position) in bases for position in range(span.start, span.end)):
            return False
    return True


def held(territory: Territory) -> set[tuple[str, int]]:
    """Return every base of a territory, checking its spans are the fewest that hold them."""
    spans = list(territory)
    for previous, span in pairwise(spans):
        assert previous.refname != span.refname or previous.end < span.start
    assert all(span.start < span.end for span in spans)
    return {(span.refname, position) for span in spans for position in range(span.start, span.end)}


def test_territories_agree_with_sets_of_bases_on_random_features() -> None:
    """Test that territories and their combinations agree with sets of bases, on random features."""
    rng = Random(42)
    for _ in range(300):
        first = [random_feature(rng) for _ in range(rng.randrange(12))]
        second = [random_feature(rng) for _ in range(rng.randrange(12))]
        territory, other = Territory(first), Territory(second)
        bases, other_bases = bases_of(first), bases_of(second)

        assert held(territory) == bases
        assert territory.length == len(bases)
        assert bool(territory) is bool(bases)
        assert held(territory | other) == bases | other_bases
        assert held(territory & other) == bases & other_bases
        assert held(territory - other) == bases - other_bases
        assert (territory - other) | (territory & other) == territory
        assert (territory == other) is (bases == other_bases)
        for refname in ("chr1", "chr2", "chr3"):
            for position in range(120):
                assert territory.contains(refname, position) is ((refname, position) in bases)
        for feature in second:
            assert (feature in territory) is is_held(feature, bases)
            assert held(feature.territory()) == bases_of([feature])
