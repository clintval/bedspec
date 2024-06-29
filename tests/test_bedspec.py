from dataclasses import dataclass

import pytest

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedGraph
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import PairBed
from bedspec import PointBed
from bedspec import SimpleBed


def test_instantiating_all_bed_types() -> None:
    """Test that we can instantiate all builtin BED types."""
    Bed2(contig="chr1", start=1)
    Bed3(contig="chr1", start=1, end=2)
    Bed4(contig="chr1", start=1, end=2, name="foo")
    Bed5(contig="chr1", start=1, end=2, name="foo", score=3)
    Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive)
    BedGraph(contig="chr1", start=1, end=2, value=0.2)
    BedPE(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
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
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    assert record.bed1 == Bed6(contig="chr1", start=1, end=2, name="foo", score=5, strand=BedStrand.Positive)  # fmt: skip  # noqa: E501
    assert record.bed2 == Bed6(contig="chr2", start=3, end=4, name="foo", score=5, strand=BedStrand.Negative)  # fmt: skip  # noqa: E501


def test_point_bed_types_have_a_territory() -> None:
    """Test that a point BED has a territory of 1-length."""
    expected = Bed3(contig="chr1", start=1, end=2)
    assert list(Bed2(contig="chr1", start=1).territory()) == [expected]


def test_point_bed_types_are_length_1() -> None:
    """Test that a point BED has a length of 1."""
    assert Bed2(contig="chr1", start=1).length == 1


def test_simple_bed_types_have_a_territory() -> None:
    """Test that simple BEDs are their own territory."""
    for record in (
        Bed3(contig="chr1", start=1, end=2),
        Bed4(contig="chr1", start=1, end=2, name="foo"),
        Bed5(contig="chr1", start=1, end=2, name="foo", score=3),
        Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive),
        BedGraph(contig="chr1", start=1, end=2, value=1.0),
    ):
        assert list(record.territory()) == [record]


def test_simple_bed_types_have_length() -> None:
    """Test that a simple BED has the right length."""
    assert Bed3(contig="chr1", start=1, end=2).length == 1
    assert Bed3(contig="chr1", start=1, end=3).length == 2
    assert Bed3(contig="chr1", start=1, end=4).length == 3


def test_simple_bed_validates_start_and_end() -> None:
    """Test that a simple BED record validates its start and end."""
    with pytest.raises(ValueError):
        Bed3(contig="chr1", start=-1, end=5)
    with pytest.raises(ValueError):
        Bed3(contig="chr1", start=5, end=5)
    with pytest.raises(ValueError):
        Bed3(contig="chr1", start=5, end=0)


def test_paired_bed_validates_start_and_end() -> None:
    """Test a paired BED record validates its start and end for both intervals."""
    # fmt: off
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=-1, end1=5, contig2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=5, end1=5, contig2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=5, end1=0, contig2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=1, end1=2, contig2="chr1", start2=-1, end2=5, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=1, end1=2, contig2="chr1", start2=5, end2=5, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=1, end1=2, contig2="chr1", start2=5, end2=0, name="foo", score=5, strand1=BedStrand.Positive, strand2=BedStrand.Positive)  # noqa: E501
    # fmt: on


def test_paired_bed_types_have_a_territory() -> None:
    """Test that paired BEDs use both their intervals as their territory."""
    record = BedPE(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    expected: list[Bed6] = [
        Bed6(contig="chr1", start=1, end=2, name="foo", score=5, strand=BedStrand.Positive),
        Bed6(contig="chr2", start=3, end=4, name="foo", score=5, strand=BedStrand.Negative),
    ]
    assert list(record.territory()) == expected


def test_that_decoding_splits_on_any_whitespace() -> None:
    """Test that we can decode a BED on arbitrary whitespace (tabs or spaces)."""
    assert Bed3.decode("   chr1 \t 1\t \t2  \n") == Bed3(contig="chr1", start=1, end=2)


def test_that_we_can_decode_all_bed_types_from_strings() -> None:
    """Test that we can decode all builtin BED types from strings."""
    # fmt: off
    assert Bed2.decode("chr1\t1") == Bed2(contig="chr1", start=1)
    assert Bed3.decode("chr1\t1\t2") == Bed3(contig="chr1", start=1, end=2)
    assert Bed4.decode("chr1\t1\t2\tfoo") == Bed4(contig="chr1", start=1, end=2, name="foo")
    assert Bed5.decode("chr1\t1\t2\tfoo\t3") == Bed5(contig="chr1", start=1, end=2, name="foo", score=3)  # noqa: E501
    assert Bed6.decode("chr1\t1\t2\tfoo\t3\t+") == Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive)  # noqa: E501
    assert BedGraph.decode("chr1\t1\t2\t0.2") == BedGraph(contig="chr1", start=1, end=2, value=0.2)
    assert BedPE.decode("chr1\t1\t2\tchr2\t3\t4\tfoo\t5\t+\t-") == BedPE(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    # fmt: on


def test_that_we_can_make_our_own_custom_point_bed() -> None:
    """Test that we can make our own custom point BED type."""

    @dataclass
    class Bed2_2(PointBed):
        """A custom BED2+2 record."""

        contig: str
        start: int
        custom1: float
        custom2: int

    decoded = Bed2_2(contig="chr1", start=1, custom1=3.0, custom2=4)
    assert Bed2_2.decode("chr1\t1\t3.0\t4") == decoded
    expected = Bed3(contig="chr1", start=1, end=2)
    assert list(decoded.territory()) == [expected]


def test_that_we_can_make_our_own_custom_simple_bed() -> None:
    """Test that we can make our own custom simple BED type."""

    @dataclass
    class Bed3_2(SimpleBed):
        """A custom BED3+2 record."""

        contig: str
        start: int
        end: int
        custom1: float
        custom2: int

    decoded = Bed3_2(contig="chr1", start=1, end=2, custom1=3.0, custom2=4)
    assert Bed3_2.decode("chr1\t1\t2\t3.0\t4") == decoded
    assert list(decoded.territory()) == [decoded]


def test_that_we_can_make_our_own_custom_paired_bed() -> None:
    """Test that we can make our own custom paired BED type."""

    @dataclass
    class PairedBed6_1(PairBed):
        """A custom BedPE6+2 record."""

        contig1: str
        start1: int
        end1: int
        contig2: str
        start2: int
        end2: int
        custom1: float

    decoded = PairedBed6_1(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        custom1=4.0,
    )

    assert PairedBed6_1.decode("chr1\t1\t2\tchr2\t3\t4\t4.0") == decoded
    territory = list(decoded.territory())
    assert territory == [Bed3(contig="chr1", start=1, end=2), Bed3(contig="chr2", start=3, end=4)]


def test_that_we_cannot_build_a_custom_bed_without_it_being_a_dataclass() -> None:
    """Test that we cannot build a custom BED type without it being a dataclass."""

    class Bed2_1(PointBed):
        """A custom BED2+1 record that is wrongly not a dataclass."""

        contig: str
        start: int
        custom1: float

    with pytest.raises(
        TypeError, match="You must annotate custom BED class definitions with @dataclass!"
    ):
        Bed2_1(contig="chr1", start=1, custom1=3.0)  # type: ignore[abstract]


def test_that_we_get_a_helpful_error_when_we_dont_have_the_right_number_of_fields() -> None:
    """Test that we get a helpful error when we cannot decode due to a wrong number of fields."""

    with pytest.raises(ValueError, match="Expected 3 fields but found 2 in record: 'chr1 1'"):
        Bed3.decode("chr1\t1")

    with pytest.raises(ValueError, match="Expected 2 fields but found 3 in record: 'chr1 1 2'"):
        Bed2.decode("chr1\t1\t2")


def test_that_we_get_a_helpful_error_when_we_cant_decode_the_types() -> None:
    """Test that we get a helpful error when we cannot decode due to having the wrong types."""

    with pytest.raises(
        TypeError,
        match=(
            r"Tried to build the BED field 'start' \(of type 'int'\) from the value 'chr1'"
            r" but couldn't for record 'chr1 chr1'"
        ),
    ):
        Bed2.decode("chr1\tchr1")
