from pathlib import Path

from typeline import Comment

from bedspec import MISSING_FIELD
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed12
from bedspec import BedColor
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter


def test_bed_reader_can_read_bed_records_from_a_path(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records from a path if it is typed."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter[Bed3](handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == "chr1\t1\t2\n"

    reader = BedReader.from_path[Bed3](tmp_path / "test.bed")
    assert list(reader) == [bed]

    reader = BedReader.from_path[Bed3](str(tmp_path / "test.bed"))
    assert list(reader) == [bed]


def test_bed_reader_can_read_bed_records_with_comments(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with comments."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter[Bed3](handle)
        writer.write_comment("track\nthis-is-fine")
        writer.write_comment("browser is mario's enemy?")
        writer.write_comment("hello mom!")
        handle.write("\n")  # empty line
        handle.write("  \n")  # empty line
        writer.write(bed)
        writer.write_comment("hello dad!")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader[Bed3](handle)) == [bed]


def test_bed_reader_can_read_optional_string_types(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with optional string types."""
    bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

    (tmp_path / "test.bed").write_text(f"chr1\t1\t2\t{MISSING_FIELD}\n")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader[Bed4](handle)) == [bed]


def test_bed_reader_can_read_optional_other_types(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with optional other types."""
    bed: Bed5 = Bed5(refname="chr1", start=1, end=2, name="foo", score=None)

    (tmp_path / "test.bed").write_text(f"chr1\t1\t2\tfoo\t{MISSING_FIELD}\n")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader[Bed5](handle)) == [bed]


def test_bed_reader_can_be_used_as_context_manager(tmp_path: Path) -> None:
    """Test that the BED reader can be used as a context manager."""
    bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

    (tmp_path / "test.bed").write_text(f"chr1\t1\t2\t{MISSING_FIELD}\n")

    with BedReader[Bed4](open(tmp_path / "test.bed")) as reader:
        assert list(reader) == [bed]


def test_we_can_roundtrip_a_bed_record_with_complex_types(tmp_path: Path) -> None:
    """Test that we can roundtrip a BED record with complex types (e.g. lists)."""
    bed12: Bed12 = Bed12(
        refname="chr1",
        start=2,
        end=10,
        name="bed12",
        score=2,
        strand=BedStrand.Positive,
        thick_start=3,
        thick_end=4,
        item_rgb=BedColor(101, 2, 32),
        block_count=2,
        block_sizes=[1, 2],
        block_starts=[0, 6],
    )

    with BedWriter.from_path[Bed12](tmp_path / "test.bed") as writer:
        writer.write(bed12)

    expected: str = "chr1\t2\t10\tbed12\t2\t+\t3\t4\t101,2,32\t2\t1,2\t0,6\n"
    assert Path(tmp_path / "test.bed").read_text() == expected

    with BedReader.from_path[Bed12](tmp_path / "test.bed") as reader:
        assert list(reader) == [bed12]


def test_we_can_roundtrip_a_bed_record_without_a_color(tmp_path: Path) -> None:
    """Test that a BED record without a color is written with a 0 color and read back as None."""
    bed12: Bed12 = Bed12(
        refname="chr1",
        start=2,
        end=10,
        name="bed12",
        score=2,
        strand=BedStrand.Positive,
        thick_start=3,
        thick_end=4,
        item_rgb=None,
        block_count=2,
        block_sizes=[1, 2],
        block_starts=[0, 6],
    )

    with BedWriter.from_path[Bed12](tmp_path / "test.bed") as writer:
        writer.write(bed12)

    expected: str = "chr1\t2\t10\tbed12\t2\t+\t3\t4\t0\t2\t1,2\t0,6\n"
    assert Path(tmp_path / "test.bed").read_text() == expected

    with BedReader.from_path[Bed12](tmp_path / "test.bed") as reader:
        assert list(reader) == [bed12]


def test_bed_reader_hands_comments_to_on_comment(tmp_path: Path) -> None:
    """Test that the BED reader hands each comment line, with its line number, to on_comment."""
    (tmp_path / "test.bed").write_text(
        "track name=test\n# made by a tool\nchr1\t1\t2\nbrowser hide all\n"
    )

    comments: list[Comment] = []
    with BedReader.from_path[Bed3](tmp_path / "test.bed", on_comment=comments.append) as reader:
        assert list(reader) == [Bed3(refname="chr1", start=1, end=2)]

    assert comments == [
        Comment(line_number=1, text="track name=test"),
        Comment(line_number=2, text="# made by a tool"),
        Comment(line_number=4, text="browser hide all"),
    ]
