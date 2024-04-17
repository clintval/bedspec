from pathlib import Path

import pytest

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedType
from bedspec.io import BedWriter


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
        ],
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
        writer.write_comment("hello mom!")
        writer.write(Bed2(contig="chr1", start=1))
        writer.write_comment("hello dad!")
        writer.write(Bed2(contig="chr2", start=2))

    expected = "# hello mom!\nchr1\t1\n# hello dad!\nchr2\t2\n"
    assert Path(tmp_path / "test.bed").read_text() == expected
