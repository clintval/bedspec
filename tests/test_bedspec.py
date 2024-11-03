from dataclasses import is_dataclass

import pytest

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
from bedspec import GenomicSpan
from bedspec import PairBed
from bedspec import PointBed
from bedspec import SimpleBed
from bedspec import Stranded
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
        BedColor(r, g, b)


def test_bed_color_from_string() -> None:
    """Test that we can build a BED color from a string."""
    assert BedColor.from_string("2,3,4") == BedColor(2, 3, 4)


def test_bed_color_from_string_raises_when_malformed() -> None:
    """Test that we raise an exception when building a BED color from a malformed string."""
    with pytest.raises(ValueError, match="Invalid string '-1,hi,4'. Expected 'int,int,int'!"):
        BedColor.from_string("-1,hi,4")


@pytest.mark.parametrize("bed_type", (PointBed, SimpleBed, PairBed))
def test_bed_type_class_hierarchy(bed_type: type[BedLike]) -> None:
    """Test that all abstract base classes are subclasses of BedLike."""
    assert issubclass(bed_type, BedLike)


@pytest.mark.parametrize("bed_type", (Bed2, Bed3, Bed4, Bed5, Bed6, BedGraph, BedPE))
def test_all_bed_types_are_dataclasses(bed_type: type[BedLike]) -> None:
    """Test that a simple BED record behaves as expected."""
    assert is_dataclass(bed_type)


def test_locatable_structural_type() -> None:
    """Test that the GenomicSpan structural type is set correctly."""
    span: GenomicSpan = Bed6(
        refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive
    )
    assert isinstance(span, GenomicSpan)


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
    Bed2(refname="chr1", start=1)
    Bed3(refname="chr1", start=1, end=2)
    Bed4(refname="chr1", start=1, end=2, name="foo")
    Bed5(refname="chr1", start=1, end=2, name="foo", score=3)
    Bed6(refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive)
    BedGraph(refname="chr1", start=1, end=2, value=0.2)
    BedPE(
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


def test_point_bed_types_have_a_territory() -> None:
    """Test that a point BED has a territory of 1-length."""
    expected = Bed3(refname="chr1", start=1, end=2)
    assert list(Bed2(refname="chr1", start=1).territory()) == [expected]


def test_point_bed_types_are_length_1() -> None:
    """Test that a point BED has a length of 1."""
    assert Bed2(refname="chr1", start=1).length == 1


def test_simple_bed_types_have_a_territory() -> None:
    """Test that simple BEDs are their own territory."""
    for record in (
        Bed3(refname="chr1", start=1, end=2),
        Bed4(refname="chr1", start=1, end=2, name="foo"),
        Bed5(refname="chr1", start=1, end=2, name="foo", score=3),
        Bed6(refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive),
        BedGraph(refname="chr1", start=1, end=2, value=1.0),
    ):
        assert list(record.territory()) == [record]


def test_simple_bed_types_have_length() -> None:
    """Test that a simple BED has the right length."""
    assert Bed3(refname="chr1", start=1, end=2).length == 1
    assert Bed3(refname="chr1", start=1, end=3).length == 2
    assert Bed3(refname="chr1", start=1, end=4).length == 3


def test_simple_bed_validates_start_and_end() -> None:
    """Test that a simple BED record validates its start and end."""
    with pytest.raises(ValueError, match="start must be greater than 0 and less than end!"):
        Bed3(refname="chr1", start=-1, end=5)
    with pytest.raises(ValueError, match="start must be greater than 0 and less than end!"):
        Bed3(refname="chr1", start=5, end=5)
    with pytest.raises(ValueError, match="start must be greater than 0 and less than end!"):
        Bed3(refname="chr1", start=5, end=0)


def test_paired_bed_validates_start_and_end() -> None:
    """Test a paired BED record validates its start and end for both intervals."""
    # fmt: off
    with pytest.raises(ValueError, match="start1 must be greater than 0 and less than end1!"):
        BedPE(refname1="chr1", start1=-1, end1=5, refname2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError, match="start1 must be greater than 0 and less than end1!"):
        BedPE(refname1="chr1", start1=5, end1=5, refname2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError, match="start1 must be greater than 0 and less than end1!"):
        BedPE(refname1="chr1", start1=5, end1=0, refname2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError, match="start2 must be greater than 0 and less than end2!"):
        BedPE(refname1="chr1", start1=1, end1=2, refname2="chr1", start2=-1, end2=5, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError, match="start2 must be greater than 0 and less than end2!"):
        BedPE(refname1="chr1", start1=1, end1=2, refname2="chr1", start2=5, end2=5, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError, match="start2 must be greater than 0 and less than end2!"):
        BedPE(refname1="chr1", start1=1, end1=2, refname2="chr1", start2=5, end2=0, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    # fmt: on


def test_paired_bed_types_have_a_territory() -> None:
    """Test that paired BEDs use both their intervals as their territory."""
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
    assert list(record.territory()) == expected


def test_bed12_validation() -> None:
    """Test that we can validate improper BED12 records."""

    def make_bed12(
        thick_start: int | None = None,
        thick_end: int | None = None,
        block_count: int | None = None,
        block_sizes: list[int] | None = None,
        block_starts: list[int] | None = None,
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

    with pytest.raises(
        ValueError, match="thick_start and thick_end must both be None or both be set!"
    ):
        make_bed12(thick_start=1, thick_end=None)
        make_bed12(thick_start=None, thick_end=2)

    with pytest.raises(
        ValueError, match="block_count, block_sizes, block_starts must all be set or unset!"
    ):
        make_bed12(block_count=1, block_sizes=None, block_starts=None)
        make_bed12(block_count=None, block_sizes=[1], block_starts=[0])
        make_bed12(block_count=1, block_sizes=None, block_starts=[0])
        make_bed12(block_count=1, block_sizes=[1], block_starts=None)

    with pytest.raises(
        ValueError, match="When set, block_count must be greater than or equal to 1!"
    ):
        make_bed12(block_count=-1, block_sizes=[1], block_starts=[0])

    with pytest.raises(
        ValueError, match="Length of block_sizes and block_starts must equal block_count!"
    ):
        make_bed12(block_count=1, block_sizes=[1], block_starts=[0, 1])
        make_bed12(block_count=1, block_sizes=[1, 2], block_starts=[0])
        make_bed12(block_count=2, block_sizes=[1], block_starts=[0])

    with pytest.raises(ValueError, match="block_starts must start with 0!"):
        make_bed12(block_count=1, block_sizes=[1], block_starts=[1])

    with pytest.raises(
        ValueError, match="All sizes in block_size must be greater than or equal to one!"
    ):
        make_bed12(block_count=1, block_sizes=[-1], block_starts=[0])

    with pytest.raises(
        ValueError, match="The last defined block's end must be equal to the BED end!"
    ):
        make_bed12(block_count=2, block_sizes=[1, 1], block_starts=[0, 4])
