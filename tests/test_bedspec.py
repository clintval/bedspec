import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pytest

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedColor
from bedspec import BedPE
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedType
from bedspec import BedWriter
from bedspec import PairBed
from bedspec import PointBed
from bedspec import SimpleBed
from bedspec import Stranded


def test_bed_strand() -> None:
    """Test that BED strands behave as string."""
    assert BedStrand("+") == BedStrand.POSITIVE
    assert BedStrand("-") == BedStrand.NEGATIVE
    assert BedStrand(".") == BedStrand.UNKNOWN
    assert str(BedStrand.POSITIVE) == "+"
    assert str(BedStrand.NEGATIVE) == "-"
    assert str(BedStrand.UNKNOWN) == "."


def test_bed_color() -> None:
    """Test the small helper class for BED color."""
    assert str(BedColor(2, 3, 4)) == "2,3,4"


def test_bed_type_class_hierarchy() -> None:
    """Test that all abstract base classes are subclasses of Bedtype."""
    for subclass in (PointBed, SimpleBed, PairBed):
        assert issubclass(subclass, BedType)


@pytest.mark.parametrize("bed_type", (Bed2, Bed3, Bed4, Bed5, Bed6, BedPE))
def test_all_bed_types_are_dataclasses(bed_type: type[BedType]) -> None:
    """Test that a simple BED record behaves as expected."""
    assert dataclasses.is_dataclass(bed_type)


def test_stranded_structural_type() -> None:
    """Test that the Stranded structural type is set correctly."""
    _: Stranded = Bed6(
        contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.POSITIVE
    )


def test_dataclass_protocol_structural_type() -> None:
    """Test that the dataclass structural type is set correctly."""
    from bedspec._bedspec import DataclassProtocol

    _: DataclassProtocol = Bed2(contig="chr1", start=1)


def test_instantiating_all_bed_types() -> None:
    """Test that we can instantiate all BED types."""
    Bed2(contig="chr1", start=1)
    Bed3(contig="chr1", start=1, end=2)
    Bed4(contig="chr1", start=1, end=2, name="foo")
    Bed5(contig="chr1", start=1, end=2, name="foo", score=3)
    Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.POSITIVE)
    BedPE(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.POSITIVE,
        strand2=BedStrand.NEGATIVE,
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
        strand1=BedStrand.POSITIVE,
        strand2=BedStrand.NEGATIVE,
    )
    assert record.bed1 == Bed6(contig="chr1", start=1, end=2, name="foo", score=5, strand=BedStrand.POSITIVE)  # fmt: skip  # noqa: E501
    assert record.bed2 == Bed6(contig="chr2", start=3, end=4, name="foo", score=5, strand=BedStrand.NEGATIVE)  # fmt: skip  # noqa: E501


def test_point_bed_types_have_a_territory() -> None:
    """Test that a point BED has a territory of 1-length."""
    expected = Bed3(contig="chr1", start=1, end=2)
    assert list(Bed2(contig="chr1", start=1).territory()) == [expected]


def test_simple_bed_types_have_a_territory() -> None:
    """Test that simple BEDs are their own territory."""
    for record in (
        Bed3(contig="chr1", start=1, end=2),
        Bed4(contig="chr1", start=1, end=2, name="foo"),
        Bed5(contig="chr1", start=1, end=2, name="foo", score=3),
        Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.POSITIVE),
    ):
        assert list(record.territory()) == [record]


def test_simple_bed_validates_start_and_end() -> None:
    """Test that a simple BED record validates its start and end."""
    with pytest.raises(ValueError):
        Bed3(contig="chr1", start=-1, end=5)
    with pytest.raises(ValueError):
        Bed3(contig="chr1", start=5, end=5)
    with pytest.raises(ValueError):
        Bed3(contig="chr1", start=5, end=0)


def test_paired_bed_validates_start_and_end() -> None:
    """Test that a simple BED record validates its start and end."""
    # fmt: off
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=-1, end1=5, contig2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.POSITIVE, strand2=BedStrand.POSITIVE)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=5, end1=5, contig2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.POSITIVE, strand2=BedStrand.POSITIVE)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=5, end1=0, contig2="chr1", start2=1, end2=2, name="foo", score=5, strand1=BedStrand.POSITIVE, strand2=BedStrand.POSITIVE)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=1, end1=2, contig2="chr1", start2=-1, end2=5, name="foo", score=5, strand1=BedStrand.POSITIVE, strand2=BedStrand.POSITIVE)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=1, end1=2, contig2="chr1", start2=5, end2=5, name="foo", score=5, strand1=BedStrand.POSITIVE, strand2=BedStrand.POSITIVE)  # noqa: E501
    with pytest.raises(ValueError):
        BedPE(contig1="chr1", start1=1, end1=2, contig2="chr1", start2=5, end2=0, name="foo", score=5, strand1=BedStrand.POSITIVE, strand2=BedStrand.POSITIVE)  # noqa: E501
    # fmt: on


def test_paired_bed_types_have_a_territory() -> None:
    """Test that simple BEDs are their own territory."""
    record = BedPE(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.POSITIVE,
        strand2=BedStrand.NEGATIVE,
    )
    expected: list[Bed6] = [
        Bed6(contig="chr1", start=1, end=2, name="foo", score=5, strand=BedStrand.POSITIVE),
        Bed6(contig="chr2", start=3, end=4, name="foo", score=5, strand=BedStrand.NEGATIVE),
    ]
    assert list(record.territory()) == expected


def test_that_decoding_splits_on_any_whitespace() -> None:
    """Test that we can decode a BED on arbitrary whitespace."""
    assert Bed3.decode("   chr1 \t 1\t \t2  \n") == Bed3(contig="chr1", start=1, end=2)


def test_all_bed_types_have_fieldnames() -> None:
    # fmt: off
    assert Bed2.decode("chr1\t1") == Bed2(contig="chr1", start=1)
    assert Bed3.decode("chr1\t1\t2") == Bed3(contig="chr1", start=1, end=2)
    assert Bed4.decode("chr1\t1\t2\tfoo") == Bed4(contig="chr1", start=1, end=2, name="foo")
    assert Bed5.decode("chr1\t1\t2\tfoo\t3") == Bed5(contig="chr1", start=1, end=2, name="foo", score=3)  # noqa: E501
    assert Bed6.decode("chr1\t1\t2\tfoo\t3\t+") == Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.POSITIVE)  # noqa: E501
    assert BedPE.decode("chr1\t1\t2\tchr2\t3\t4\tfoo\t5\t+\t-") == BedPE(
        contig1="chr1",
        start1=1,
        end1=2,
        contig2="chr2",
        start2=3,
        end2=4,
        name="foo",
        score=5,
        strand1=BedStrand.POSITIVE,
        strand2=BedStrand.NEGATIVE,
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

    with pytest.raises(TypeError, match="You must mark custom BED records with @dataclass."):
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


@pytest.mark.parametrize(
    "bed,expected",
    [
        [Bed2(contig="chr1", start=1), "chr1\t1\n"],
        [Bed3(contig="chr1", start=1, end=2), "chr1\t1\t2\n"],
        [Bed4(contig="chr1", start=1, end=2, name="foo"), "chr1\t1\t2\tfoo\n"],
        [Bed5(contig="chr1", start=1, end=2, name="foo", score=3), "chr1\t1\t2\tfoo\t3\n"],
        [
            Bed6(contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.POSITIVE),
            "chr1\t1\t2\tfoo\t3\t+\n",
        ],  # fmt: skip
        [
            BedPE(
                contig1="chr1",
                start1=1,
                end1=2,
                contig2="chr2",
                start2=3,
                end2=4,
                name="foo",
                score=5,
                strand1=BedStrand.POSITIVE,
                strand2=BedStrand.NEGATIVE,
            ),
            "chr1\t1\t2\tchr2\t3\t4\tfoo\t5\t+\t-\n",
        ],
    ],
)
def test_bed_writer_can_write_all_bed_types(bed: BedType, expected: str, tmp_path: Path) -> None:
    """Test that the BED writer can write all BED types."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_can_write_all_at_once(tmp_path: Path) -> None:
    """Test that the BED writer can write multiple BED records at once."""
    expected: str = "chr1\t1\t2\nchr2\t3\t4\n"

    def records() -> Iterator[Bed3]:
        yield Bed3(contig="chr1", start=1, end=2)
        yield Bed3(contig="chr2", start=3, end=4)

    with open(tmp_path / "test.bed", "w") as handle:
        BedWriter[Bed3](handle).write_all(records())

    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_remembers_the_type_it_will_write(tmp_path: Path) -> None:
    """Test that the BED writer remembers the type it can only write."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(Bed2(contig="chr1", start=1))
        assert writer.bed_kind is Bed2
        with pytest.raises(
            TypeError,
            match=(
                "BedWriter can only continue to write features of the same type. Will not write a"
                " Bed3 after a Bed2"
            ),
        ):
            writer.write(Bed3(contig="chr1", start=1, end=2))


def test_bed_writer_remembers_the_type_it_will_write_generic(tmp_path: Path) -> None:
    """Test that the generically parameterized BED writer remembers the type it can only write."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter[Bed2](handle)
        assert writer.bed_kind is Bed2
        with pytest.raises(
            TypeError,
            match=(
                "BedWriter can only continue to write features of the same type. Will not write a"
                " Bed3 after a Bed2"
            ),
        ):
            writer.write(Bed3(contig="chr1", start=1, end=2))  # type: ignore[arg-type]


def test_bed_writer_write_comment_with_prefix_pound_symbol(tmp_path: Path) -> None:
    """Test that we can write comments that have a leading pound symbol."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter[Bed2](handle)
        writer.write_comment("# hello mom!")
        writer.write(Bed2(contig="chr1", start=1))
        writer.write_comment("# hello dad!")
        writer.write(Bed2(contig="chr2", start=2))

    expected = "# hello mom!\nchr1\t1\n# hello dad!\nchr2\t2\n"
    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_write_comment_without_prefix_pound_symbol(tmp_path: Path) -> None:
    """Test that we can write comments that do not have a leading pound symbol."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter[Bed2](handle)
        writer.write_comment("track this-is-fine")
        writer.write_comment("browser is mario's enemy?")
        writer.write_comment("hello mom!")
        writer.write(Bed2(contig="chr1", start=1))
        writer.write_comment("hello dad!")
        writer.write(Bed2(contig="chr2", start=2))

    expected = (
        "track this-is-fine\n"
        "browser is mario's enemy?\n"
        "# hello mom!\n"
        "chr1\t1\n"
        "# hello dad!\n"
        "chr2\t2\n"
    )

    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_reader_can_read_bed_records_if_typed(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records if the reader is typed."""

    bed: Bed3 = Bed3(contig="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == "chr1\t1\t2\n"

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader[Bed3](handle)) == [bed]


def test_bed_reader_can_read_bed_records_from_a_path(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records from a path if it is typed."""

    bed: Bed3 = Bed3(contig="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == "chr1\t1\t2\n"

    reader = BedReader[Bed3].from_path(tmp_path / "test.bed", bed_kind=Bed3)
    assert list(reader) == [bed]

    reader = BedReader[Bed3].from_path(str(tmp_path / "test.bed"), bed_kind=Bed3)
    assert list(reader) == [bed]


def test_bed_reader_can_raises_exception_if_not_typed(tmp_path: Path) -> None:
    """Test that the BED reader raises an exception if it is not typed."""

    bed: Bed3 = Bed3(contig="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == "chr1\t1\t2\n"

    with open(tmp_path / "test.bed", "r") as handle:
        with pytest.raises(
            NotImplementedError,
            match="Untyped reading is not yet supported!",
        ):
            list(BedReader(handle))


def test_bed_reader_can_read_bed_records_with_comments(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with comments."""

    bed: Bed3 = Bed3(contig="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write_comment("track this-is-fine")
        writer.write_comment("browser is mario's enemy?")
        writer.write_comment("hello mom!")
        handle.write("\n")  # empty line
        handle.write(" \t\n")  # empty line
        writer.write(bed)
        writer.write_comment("hello dad!")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader[Bed3](handle)) == [bed]
