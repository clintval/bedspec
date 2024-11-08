from pathlib import Path

import pytest

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedGraph
from bedspec import BedLike
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedWriter


# fmt: off
@pytest.mark.parametrize(
    "bed,expected",
    [
        [Bed2(refname="chr1", start=1), "chr1\t1\n"],
        [Bed3(refname="chr1", start=1, end=2), "chr1\t1\t2\n"],
        [Bed4(refname="chr1", start=1, end=2, name="foo"), "chr1\t1\t2\tfoo\n"],
        [Bed5(refname="chr1", start=1, end=2, name="foo", score=3), "chr1\t1\t2\tfoo\t3\n"],
        [Bed6(refname="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive), "chr1\t1\t2\tfoo\t3\t+\n"],  # noqa: E501
        [BedGraph(refname="chr1", start=1, end=2, value=0.2), "chr1\t1\t2\t0.2\n"],
        [
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
            ),
            "chr1\t1\t2\tchr2\t3\t4\tfoo\t5\t+\t-\n",
        ],
    ],
)
# fmt: on
def test_bed_writer_can_write_all_bed_types(bed: BedLike, expected: str, tmp_path: Path) -> None:
    """Test that the BED writer can write all BED types."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle, type(bed))
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_can_be_closed(tmp_path: Path) -> None:
    """Test that we can close a BED writer."""
    path: Path = tmp_path / "test.bed"
    writer = BedWriter(open(path, "w"), Bed3)
    writer.write(Bed3(refname="chr1", start=1, end=2))
    writer.close()

    with pytest.raises(ValueError, match="I/O operation on closed file"):
        writer.write(Bed3(refname="chr1", start=1, end=2))


def test_bed_writer_can_write_bed_records_from_a_path(tmp_path: Path) -> None:
    """Test that the BED write can write BED records from a path if it is typed."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with BedWriter.from_path(tmp_path / "test1.bed", Bed3) as writer:
        writer.write(bed)

    assert (tmp_path / "test1.bed").read_text() == "chr1\t1\t2\n"

    with BedWriter.from_path(str(tmp_path / "test2.bed"), Bed3) as writer:
        writer.write(bed)

    assert (tmp_path / "test2.bed").read_text() == "chr1\t1\t2\n"


def test_bed_writer_remembers_the_type_it_will_write(tmp_path: Path) -> None:
    """Test that the BED writer remembers the type it can only write."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle, Bed2)
        writer.write(Bed2(refname="chr1", start=1))
        with pytest.raises(
            ValueError,
            match="Expected Bed2 but found Bed3!",
        ):
            writer.write(Bed3(refname="chr1", start=1, end=2))


def test_bed_writer_remembers_the_type_it_will_write_generic(tmp_path: Path) -> None:
    """Test that the generically parameterized BED writer remembers the type it can only write."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter(handle, Bed2)
        writer.write(Bed2("chr1", 1))
        with pytest.raises(
            ValueError,
            match="Expected Bed2 but found Bed3!",
        ):
            writer.write(Bed3(refname="chr1", start=1, end=2))  # type: ignore[arg-type]


def test_bed_writer_write_comment_with_prefix_pound_symbol(tmp_path: Path) -> None:
    """Test that we can write comments that have a leading pound symbol."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter(handle, Bed2)
        writer.write_comment("# hello mom!")
        writer.write(Bed2(refname="chr1", start=1))
        writer.write_comment("# hello\ndad!")
        writer.write(Bed2(refname="chr2", start=2))

    expected = "# hello mom!\nchr1\t1\n# hello\n# dad!\nchr2\t2\n"
    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_write_comment_without_prefix_pound_symbol(tmp_path: Path) -> None:
    """Test that we can write comments that do not have a leading pound symbol."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter(handle, Bed2)
        writer.write_comment("track this-is-fine")
        writer.write_comment("browser is mario's enemy?")
        writer.write_comment("hello\nmom!")
        writer.write(Bed2(refname="chr1", start=1))
        writer.write_comment("hello dad!")
        writer.write(Bed2(refname="chr2", start=2))

    expected = (
        "track this-is-fine\n"
        "browser is mario's enemy?\n"
        "# hello\n"
        "# mom!\n"
        "chr1\t1\n"
        "# hello dad!\n"
        "chr2\t2\n"
    )

    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_can_be_used_as_context_manager(tmp_path: Path) -> None:
    """Test that the BED writer can be used as a context manager."""
    with BedWriter(open(tmp_path / "test.bed", "w"), Bed2) as handle:
        handle.write(Bed2(refname="chr1", start=1))
        handle.write(Bed2(refname="chr2", start=2))

    expected = "chr1\t1\nchr2\t2\n"
    assert Path(tmp_path / "test.bed").read_text() == expected
