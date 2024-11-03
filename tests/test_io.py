from pathlib import Path

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
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec._io import MISSING_FIELD


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
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_can_be_closed(tmp_path: Path) -> None:
    """Test that we can close a BED writer."""
    path: Path = tmp_path / "test.bed"
    writer = BedWriter[Bed3](open(path, "w"))
    writer.write(Bed3(refname="chr1", start=1, end=2))
    writer.close()

    with pytest.raises(ValueError, match="I/O operation on closed file"):
        writer.write(Bed3(refname="chr1", start=1, end=2))


def test_bed_writer_can_write_bed_records_from_a_path(tmp_path: Path) -> None:
    """Test that the BED write can write BED records from a path if it is typed."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with BedWriter.from_path(tmp_path / "test1.bed") as writer:
        writer.write(bed)

    assert (tmp_path / "test1.bed").read_text() == "chr1\t1\t2\n"

    with BedWriter.from_path(str(tmp_path / "test2.bed")) as writer:
        writer.write(bed)

    assert (tmp_path / "test2.bed").read_text() == "chr1\t1\t2\n"


def test_bed_writer_remembers_the_type_it_will_write(tmp_path: Path) -> None:
    """Test that the BED writer remembers the type it can only write."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(Bed2(refname="chr1", start=1))
        assert writer.bed_type is Bed2
        with pytest.raises(
            TypeError,
            match=(
                "BedWriter can only continue to write features of the same type. Will not write a"
                " Bed3 after a Bed2"
            ),
        ):
            writer.write(Bed3(refname="chr1", start=1, end=2))


def test_bed_writer_remembers_the_type_it_will_write_generic(tmp_path: Path) -> None:
    """Test that the generically parameterized BED writer remembers the type it can only write."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter[Bed2](handle)
        writer.write(Bed2("chr1", 1))
        with pytest.raises(
            TypeError,
            match=(
                "BedWriter can only continue to write features of the same type. Will not write a"
                " Bed3 after a Bed2"
            ),
        ):
            writer.write(Bed3(refname="chr1", start=1, end=2))  # type: ignore[arg-type]


def test_bed_writer_write_comment_with_prefix_pound_symbol(tmp_path: Path) -> None:
    """Test that we can write comments that have a leading pound symbol."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter[Bed2](handle)
        writer.write_comment("# hello mom!")
        writer.write(Bed2(refname="chr1", start=1))
        writer.write_comment("# hello\ndad!")
        writer.write(Bed2(refname="chr2", start=2))

    expected = "# hello mom!\nchr1\t1\n# hello\n# dad!\nchr2\t2\n"
    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_writer_write_comment_without_prefix_pound_symbol(tmp_path: Path) -> None:
    """Test that we can write comments that do not have a leading pound symbol."""
    with open(tmp_path / "test.bed", "w") as handle:
        writer = BedWriter[Bed2](handle)
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
    with BedWriter[Bed2](open(tmp_path / "test.bed", "w")) as handle:
        handle.write(Bed2(refname="chr1", start=1))
        handle.write(Bed2(refname="chr2", start=2))

    expected = "chr1\t1\nchr2\t2\n"
    assert Path(tmp_path / "test.bed").read_text() == expected


def test_bed_reader_can_read_bed_records_if_typed(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records if the reader is typed."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == "chr1\t1\t2\n"

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader(handle, Bed3)) == [bed]


def test_bed_reader_can_be_closed(tmp_path: Path) -> None:
    """Test that we can close a BED reader."""
    path: Path = tmp_path / "test.bed"
    path.touch()
    reader = BedReader(open(path), Bed3)
    reader.close()

    with pytest.raises(ValueError, match="I/O operation on closed file"):
        next(iter(reader))


def test_bed_reader_can_read_bed_records_from_a_path(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records from a path if it is typed."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write(bed)

    assert Path(tmp_path / "test.bed").read_text() == "chr1\t1\t2\n"

    reader = BedReader.from_path(tmp_path / "test.bed", Bed3)
    assert list(reader) == [bed]

    reader = BedReader.from_path(str(tmp_path / "test.bed"), Bed3)
    assert list(reader) == [bed]


def test_bed_reader_can_read_bed_records_with_comments(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with comments."""
    bed: Bed3 = Bed3(refname="chr1", start=1, end=2)

    with open(tmp_path / "test.bed", "w") as handle:
        writer: BedWriter = BedWriter(handle)
        writer.write_comment("track\nthis-is-fine")
        writer.write_comment("browser is mario's enemy?")
        writer.write_comment("hello mom!")
        handle.write("\n")  # empty line
        handle.write(" \t\n")  # empty line
        writer.write(bed)
        writer.write_comment("hello dad!")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader(handle, Bed3)) == [bed]


def test_bed_reader_can_read_optional_string_types(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with optional string types."""
    bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

    (tmp_path / "test.bed").write_text(f"chr1\t1\t2\t{MISSING_FIELD}\n")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader(handle, Bed4)) == [bed]


def test_bed_reader_can_read_optional_other_types(tmp_path: Path) -> None:
    """Test that the BED reader can read BED records with optional other types."""
    bed: Bed5 = Bed5(refname="chr1", start=1, end=2, name="foo", score=None)

    (tmp_path / "test.bed").write_text(f"chr1\t1\t2\tfoo\t{MISSING_FIELD}\n")

    with open(tmp_path / "test.bed", "r") as handle:
        assert list(BedReader(handle, Bed5)) == [bed]


def test_bed_reader_can_be_used_as_context_manager(tmp_path: Path) -> None:
    """Test that the BED reader can be used as a context manager."""
    bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

    (tmp_path / "test.bed").write_text(f"chr1\t1\t2\t{MISSING_FIELD}\n")

    with BedReader(open(tmp_path / "test.bed"), Bed4) as reader:
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

    with BedWriter.from_path(tmp_path / "test.bed") as writer:
        writer.write(bed12)

    expected: str = "chr1\t2\t10\tbed12\t2\t+\t3\t4\t101,2,32\t2\t1,2\t0,6\n"
    assert Path(tmp_path / "test.bed").read_text() == expected

    with BedReader.from_path(tmp_path / "test.bed", Bed12) as reader:
        assert list(reader) == [bed12]

# @pytest.mark.parametrize("ext", _ALL_GZIP_COMPATIBLE_EXTENSIONS)
# def test_bed_reader_can_read_gzip_compressed(tmp_path: Path, ext: str) -> None:
#     """Test that the BED reader can read gzip compressed paths."""
#     bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

#     with gzip.open(tmp_path / ("test.bed" + ext), "wt") as handle:
#         handle.write(f"chr1\t1\t2\t{MISSING_FIELD}\n")

#     with BedReader[Bed4](gzip.open(tmp_path / ("test.bed" + ext), "rt")) as reader:
#         assert list(reader) == [bed]


# @pytest.mark.parametrize("ext", _ALL_GZIP_COMPATIBLE_EXTENSIONS)
# def test_bed_reader_can_read_gzip_compressed_generic(tmp_path: Path, ext: str) -> None:
#     """Test that the BED reader can read gzip compressed paths."""
#     bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

#     with gzip.open(tmp_path / ("test.bed" + ext), "wt") as handle:
#         handle.write(f"chr1\t1\t2\t{MISSING_FIELD}\n")

#     with BedReader.from_path(tmp_path / ("test.bed" + ext)) as reader:
#         assert list(reader) == [bed]

# @pytest.mark.parametrize("ext", _GZIP_EXTENSIONS)
# def test_bed_writer_can_write_gzip_compressed(tmp_path: Path, ext: str) -> None:
#     """Test that the BED writer can write gzip compressed paths."""
#     bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

#     with BedWriter[Bed4](gzip.open(tmp_path / ("test.bed" + ext), "wt")) as writer:
#         writer.write(bed)

#     with BedReader[Bed4](gzip.open(tmp_path / ("test.bed" + ext), "rt")) as reader:
#         assert list(reader) == [bed]

# @pytest.mark.parametrize("ext", _GZIP_EXTENSIONS)
# def test_bed_writer_can_write_gzip_compressed_generic(tmp_path: Path, ext: str) -> None:
#     """Test that the BED writer can write gzip compressed paths."""
#     bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

#     with BedWriter.from_path(tmp_path / ("test.bed" + ext)) as writer:
#         writer.write(bed)

#     with BedReader[Bed4](gzip.open(tmp_path / ("test.bed" + ext), "rt")) as reader:
#         assert list(reader) == [bed]

# @pytest.mark.parametrize("ext", _BGZIP_EXTENSIONS)
# def test_bed_writer_can_write_block_gzip_compressed_generic(tmp_path: Path, ext: str) -> None:
#     """Test that the BED writer can write gzip compressed paths."""
#     bed: Bed4 = Bed4(refname="chr1", start=1, end=2, name=None)

#     with BedWriter.from_path(tmp_path / ("test.bed" + ext)) as writer:
#         writer.write(bed)

#     with BedReader[Bed4](gzip.open(tmp_path / ("test.bed" + ext), "rt")) as reader:
#         assert list(reader) == [bed]
