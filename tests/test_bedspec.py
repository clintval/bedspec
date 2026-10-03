from dataclasses import is_dataclass

import pytest
from typing_extensions import override

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedColor
from bedspec import BedGraph
from bedspec import BedLike
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import PairBed
from bedspec import PointBed
from bedspec import ReferenceSpan
from bedspec import SimpleBed
from bedspec import Stranded
from bedspec import Territory
from bedspec._bedspec import DataclassInstance


def test_bed_strand() -> None:
    """Test that BED strands behave as string."""
    assert BedStrand("+") == BedStrand.Positive
    assert BedStrand("-") == BedStrand.Negative
    assert str(BedStrand.Positive) == "+"
    assert str(BedStrand.Negative) == "-"


def test_bed_strand_opposite() -> None:
    """Test that we return an opposite BED strand."""
    assert BedStrand.Positive.opposite() == BedStrand.Negative
    assert BedStrand.Negative.opposite() == BedStrand.Positive


def test_bed_color() -> None:
    """Test the small helper class for BED color."""
    assert str(BedColor(2, 3, 4)) == "2,3,4"


@pytest.mark.parametrize(
    "r,g,b",
    [
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
        (256, 0, 0),
        (0, 256, 0),
        (0, 0, 256),
    ],
)
def test_bed_color_validation(r: int, g: int, b: int) -> None:
    """Test that an invalid BED color cannot be made."""
    with pytest.raises(
        ValueError, match=r"RGB color values must be in the range \[0, 255\] but found"
    ):
        _ = BedColor(r, g, b)


def test_bed_color_from_string() -> None:
    """Test that we can build a BED color from a string."""
    assert BedColor.from_string("2,3,4") == BedColor(2, 3, 4)


def test_bed_color_from_string_raises_when_malformed() -> None:
    """Test that we raise an exception when building a BED color from a malformed string."""
    with pytest.raises(ValueError, match="Invalid string '-1,hi,4'. Expected 'int,int,int'!"):
        _ = BedColor.from_string("-1,hi,4")


@pytest.mark.parametrize("bed_type", (PointBed, SimpleBed, PairBed))
def test_bed_type_class_hierarchy(bed_type: type[BedLike]) -> None:
    """Test that all abstract base classes are subclasses of BedLike."""
    assert issubclass(bed_type, BedLike)


@pytest.mark.parametrize("bed_type", (Bed2, Bed3, Bed4, Bed5, Bed6, BedGraph, BedPE))
def test_all_bed_types_are_dataclasses(bed_type: type[BedLike]) -> None:
    """Test that a simple BED record behaves as expected."""
    assert is_dataclass(bed_type)


def test_locatable_structural_type() -> None:
    """Test that the ReferenceSpan structural type is set correctly."""
    span: ReferenceSpan = Bed6(
        refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive
    )
    assert isinstance(span, ReferenceSpan)


def test_stranded_structural_type() -> None:
    """Test that the Stranded structural type is set correctly."""
    stranded: Stranded = Bed6(
        refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive
    )
    assert isinstance(stranded, Stranded)


def test_dataclass_protocol_structural_type() -> None:
    """Test that the dataclass structural type is set correctly."""
    bed: DataclassInstance = Bed2(refname="chr1", start=1)
    assert isinstance(bed, DataclassInstance)


def test_instantiating_all_bed_types() -> None:
    """Test that we can instantiate all builtin BED types."""
    _ = Bed2(refname="chr1", start=1)
    _ = Bed3(refname="chr1", start=1, end=2)
    _ = Bed4(refname="chr1", start=1, end=2, name="foo")
    _ = Bed5(refname="chr1", start=1, end=2, name="foo", score=3)
    _ = Bed6(refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive)
    _ = BedGraph(refname="chr1", start=1, end=2, value=0.2)
    _ = BedPE(
        refname1="chr1",
        start1=1,
        end1=2,
        refname2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )


def test_paired_bed_has_two_interval_properties() -> None:
    """Test that a paired BED has two BED intervals as properties."""
    record = BedPE(
        refname1="chr1",
        start1=1,
        end1=2,
        refname2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    assert record.bed1 == Bed6(refname="chr1", start=1, end=2, name="foo", score=5, strand=BedStrand.Positive)  # fmt: skip  # noqa: E501
    assert record.bed2 == Bed6(refname="chr2", start=3, end=4, name="foo", score=5, strand=BedStrand.Negative)  # fmt: skip  # noqa: E501


def test_point_bed_types_have_a_span() -> None:
    """Test that a point BED has one span of 1-length."""
    expected = Bed3(refname="chr1", start=1, end=2)
    assert list(Bed2(refname="chr1", start=1).spans()) == [expected]


def test_point_bed_types_are_length_1() -> None:
    """Test that a point BED has a length of 1."""
    assert len(Bed2(refname="chr1", start=1)) == 1


def test_simple_bed_types_have_a_span() -> None:
    """Test that simple BEDs are their own only span."""
    for record in (
        Bed3(refname="chr1", start=1, end=2),
        Bed4(refname="chr1", start=1, end=2, name="foo"),
        Bed5(refname="chr1", start=1, end=2, name="foo", score=3),
        Bed6(refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive),
        BedGraph(refname="chr1", start=1, end=2, value=1.0),
    ):
        assert list(record.spans()) == [record]


def test_simple_bed_types_have_length() -> None:
    """Test that a simple BED has the right length."""
    assert len(Bed3(refname="chr1", start=1, end=2)) == 1
    assert len(Bed3(refname="chr1", start=1, end=3)) == 2
    assert len(Bed3(refname="chr1", start=1, end=4)) == 3


def test_paired_bed_types_have_two_spans() -> None:
    """Test that paired BEDs have both their intervals as spans, each with its own strand."""
    record = BedPE(
        refname1="chr1",
        start1=1,
        end1=2,
        refname2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    expected: list[Bed6] = [
        Bed6(refname="chr1", start=1, end=2, name="foo", score=5, strand=BedStrand.Positive),
        Bed6(refname="chr2", start=3, end=4, name="foo", score=5, strand=BedStrand.Negative),
    ]
    assert list(record.spans()) == expected


def bed12(
    start: int,
    end: int,
    blocks: tuple[tuple[int, int], ...] | None,
    strand: BedStrand | None = BedStrand.Negative,
) -> Bed12:
    """Return a BED12 record on chr1 with the given blocks, each a start offset and a size."""
    return Bed12(
        refname="chr1",
        start=start,
        end=end,
        name="tx",
        score=7,
        strand=strand,
        thick_start=None,
        thick_end=None,
        item_rgb=None,
        block_count=None if blocks is None else len(blocks),
        block_sizes=None if blocks is None else tuple(size for _, size in blocks),
        block_starts=None if blocks is None else tuple(offset for offset, _ in blocks),
    )


def test_bed12_spans_are_its_blocks_on_its_strand() -> None:
    """Test that the spans of a BED12 record are its blocks, with its name, score, and strand."""
    record = bed12(100, 200, ((0, 10), (40, 20), (90, 10)))
    assert list(record.spans()) == [
        Bed6(refname="chr1", start=100, end=110, name="tx", score=7, strand=BedStrand.Negative),
        Bed6(refname="chr1", start=140, end=160, name="tx", score=7, strand=BedStrand.Negative),
        Bed6(refname="chr1", start=190, end=200, name="tx", score=7, strand=BedStrand.Negative),
    ]


def test_bed12_spans_keep_abutting_blocks_apart() -> None:
    """Test that blocks that abut are each a span of their own."""
    record = bed12(0, 10, ((0, 4), (4, 6)), strand=None)
    assert [(span.start, span.end) for span in record.spans()] == [(0, 4), (4, 10)]


def test_a_bed12_without_blocks_is_its_own_span() -> None:
    """Test that a BED12 record without blocks is its own only span."""
    record = bed12(100, 200, None)
    assert list(record.spans()) == [record]


def test_the_territory_of_a_bed12_holds_its_blocks_and_not_its_introns() -> None:
    """Test that the territory of a BED12 record holds the bases of its blocks only."""
    territory = bed12(100, 200, ((0, 10), (40, 20), (90, 10))).territory()
    assert isinstance(territory, Territory)
    assert list(territory) == [
        Bed3("chr1", start=100, end=110),
        Bed3("chr1", start=140, end=160),
        Bed3("chr1", start=190, end=200),
    ]
    assert territory.length == 40
    assert not territory.contains("chr1", 120)


def test_the_territory_of_a_record_is_the_territory_of_its_spans() -> None:
    """Test that a record's territory holds the bases of its spans, without their strands."""
    pair = BedPE(
        refname1="chr2",
        start1=10,
        end1=20,
        refname2="chr1",
        start2=15,
        end2=30,
        name=None,
        score=None,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    assert pair.territory() == Territory(pair.spans())
    assert list(pair.territory()) == [
        Bed3("chr2", start=10, end=20),
        Bed3("chr1", start=15, end=30),
    ]
    assert Bed2("chr1", 5).territory() == Territory([Bed3("chr1", start=5, end=6)])
    assert Bed3("chr1", start=5, end=5).territory() == Territory()
    record = Bed6("chr1", start=1, end=9, name="a", score=1, strand=BedStrand.Positive)
    assert list(record.territory()) == [Bed3("chr1", start=1, end=9)]


def test_a_custom_bed_class_must_override_spans_and_not_territory() -> None:
    """Test that a custom BED class that overrides territory() is refused when it is defined."""
    with pytest.raises(TypeError, match=r"Override spans\(\) in custom BED class definitions"):

        class Bad(Bed3):  # pyright: ignore[reportUnusedClass]
            @override  # type: ignore[misc]
            def territory(self) -> Territory:  # pyright: ignore[reportIncompatibleMethodOverride]
                return Territory()


def test_bed12_validation() -> None:
    """Test that we can validate improper BED12 records."""

    def make_bed12(
        thick_start: int | None = None,
        thick_end: int | None = None,
        block_count: int | None = None,
        block_sizes: tuple[int, ...] | None = None,
        block_starts: tuple[int, ...] | None = None,
    ) -> Bed12:
        return Bed12(
            refname="chr1",
            start=2,
            end=10,
            name="bed12",
            score=2,
            strand=BedStrand.Positive,
            thick_start=thick_start,
            thick_end=thick_end,
            item_rgb=BedColor(101, 2, 32),
            block_count=block_count,
            block_sizes=block_sizes,
            block_starts=block_starts,
        )

    with pytest.raises(ValueError, match="end must be greater than or equal to start!"):
        _ = Bed12(
            refname="chr1",
            start=2,
            end=1,
            name="bed12",
            score=2,
            strand=BedStrand.Positive,
            thick_start=None,
            thick_end=None,
            item_rgb=BedColor(101, 2, 32),
            block_count=None,
            block_sizes=None,
            block_starts=None,
        )

    with pytest.raises(
        ValueError, match="thick_start and thick_end must both be None or both be set!"
    ):
        _ = make_bed12(thick_start=1, thick_end=None)
        _ = make_bed12(thick_start=None, thick_end=2)

    with pytest.raises(
        ValueError, match="block_count, block_sizes, block_starts must all be set or unset!"
    ):
        _ = make_bed12(block_count=1, block_sizes=None, block_starts=None)
        _ = make_bed12(block_count=None, block_sizes=(1,), block_starts=(0,))
        _ = make_bed12(block_count=1, block_sizes=None, block_starts=(0,))
        _ = make_bed12(block_count=1, block_sizes=(1,), block_starts=None)

    with pytest.raises(
        ValueError, match="When set, block_count must be greater than or equal to 1!"
    ):
        _ = make_bed12(block_count=-1, block_sizes=(1,), block_starts=(0,))

    with pytest.raises(
        ValueError, match="Length of block_sizes and block_starts must equal block_count!"
    ):
        _ = make_bed12(block_count=1, block_sizes=(1,), block_starts=(0, 1))
        _ = make_bed12(block_count=1, block_sizes=(1, 2), block_starts=(0,))
        _ = make_bed12(block_count=2, block_sizes=(1,), block_starts=(0,))

    with pytest.raises(ValueError, match="block_starts must start with 0!"):
        _ = make_bed12(block_count=1, block_sizes=(1,), block_starts=(1,))

    with pytest.raises(
        ValueError, match="All sizes in block_size must be greater than or equal to one!"
    ):
        _ = make_bed12(block_count=1, block_sizes=(-1,), block_starts=(0,))

    with pytest.raises(
        ValueError, match="The last defined block's end must be equal to the BED end!"
    ):
        _ = make_bed12(block_count=2, block_sizes=(1, 1), block_starts=(0, 4))


def test_make_bedpe_from_pair_of_bed6() -> None:
    """Test that we can make a BedPE record from two Bed6 records."""
    bed1 = Bed6(
        refname="chr1",
        start=1,
        end=2,
        name="bed1",
        score=1,
        strand=BedStrand.Positive,
    )

    bed2 = Bed6(
        refname="chr2",
        start=3,
        end=4,
        name="bed2",
        score=2,
        strand=BedStrand.Negative,
    )

    expected = BedPE(
        refname1="chr1",
        start1=1,
        end1=2,
        refname2="chr2",
        start2=3,
        end2=4,
        name="bed3",
        score=3,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )

    assert BedPE.from_bed6(bed1, bed2, name="bed3", score=3) == expected
