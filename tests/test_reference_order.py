"""Tests of ordering records by reference sequence, then start, then end."""

import gzip
import struct
from dataclasses import dataclass
from pathlib import Path

import pytest

from bedspec import Bed3
from bedspec import Bed12
from bedspec import BedPE
from bedspec import ReferenceOrder
from bedspec import register_locatable

ORDER = ReferenceOrder(["chr2", "chr1", "chrX"])


@dataclass(frozen=True)
class OneBased:
    """A record with a refname, start and end in 1-based closed coordinates."""

    refname: str
    start: int
    end: int


@dataclass(frozen=True)
class Unplaced:
    """A record placed on no reference."""


def test_records_sort_by_reference_order_then_start_then_end() -> None:
    """Test that records sort by their reference's place in the order, then start, then end."""
    records = [
        Bed3("chr1", start=5, end=9),
        Bed3("chrX", start=0, end=1),
        Bed3("chr2", start=7, end=8),
        Bed3("chr1", start=5, end=6),
        Bed3("chr1", start=2, end=3),
    ]
    assert sorted(records, key=ORDER) == [
        Bed3("chr2", start=7, end=8),
        Bed3("chr1", start=2, end=3),
        Bed3("chr1", start=5, end=6),
        Bed3("chr1", start=5, end=9),
        Bed3("chrX", start=0, end=1),
    ]


def test_a_record_is_keyed_by_its_reference_place_start_and_end() -> None:
    """Test that a record's key is its reference's place in the order, then its start and end."""
    assert ORDER(Bed3("chr1", start=5, end=9)) == (1, 5, 9)


def test_a_bed12_is_keyed_by_its_whole_extent_and_a_bedpe_by_its_first_end() -> None:
    """Test that a record with a refname, start and end is keyed by them, and others by a span."""
    transcript = Bed12(
        "chrX",
        start=100,
        end=210,
        name=None,
        score=None,
        strand=None,
        thick_start=None,
        thick_end=None,
        item_rgb=None,
        block_count=2,
        block_sizes=(10, 10),
        block_starts=(0, 100),
    )
    pair = BedPE(
        refname1="chr1",
        start1=30,
        end1=40,
        refname2="chr2",
        start2=0,
        end2=10,
        name=None,
        score=None,
        strand1=None,
        strand2=None,
    )
    assert ORDER(transcript) == (2, 100, 210)
    assert ORDER(pair) == (1, 30, 40)


def test_a_registered_record_is_keyed_by_its_registered_spans() -> None:
    """Test that a registration beats a record's own refname, start and end."""
    with register_locatable(
        OneBased, lambda record: [Bed3(record.refname, start=record.start - 1, end=record.end)]
    ):
        assert ORDER(OneBased("chr1", 10, 20)) == (1, 9, 20)


def test_a_record_with_no_spans_sorts_after_every_other() -> None:
    """Test that a record found by nothing sorts after records on every reference."""
    unplaced = Unplaced()
    with register_locatable(Unplaced, lambda _: []):
        assert ORDER(unplaced) == (3, 0, 0)
        assert sorted([unplaced, Bed3("chrX", start=5, end=6)], key=ORDER)[-1] is unplaced


def test_a_reference_not_in_the_order_raises() -> None:
    """Test that a record on a reference the order does not name raises, naming the reference."""
    with pytest.raises(ValueError, match="'chr3' is not in this order"):
        _ = ORDER(Bed3("chr3", start=0, end=1))


def test_a_reference_named_twice_raises() -> None:
    """Test that an order naming a reference twice raises, naming the reference."""
    with pytest.raises(ValueError, match="'chr1' is named twice"):
        _ = ReferenceOrder(["chr1", "chr2", "chr1"])


def write_bam_header(path: Path, references: list[tuple[str, int]]) -> None:
    """Write a BAM file holding only a header with these references."""
    text = "@HD\tVN:1.6\n"
    body = b"BAM\x01" + struct.pack("<i", len(text)) + text.encode()
    body += struct.pack("<i", len(references))
    for name, length in references:
        encoded = name.encode() + b"\x00"
        body += struct.pack("<i", len(encoded)) + encoded + struct.pack("<i", length)
    _ = path.write_bytes(gzip.compress(body))


HEADER = "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:200\n"


@pytest.mark.parametrize(
    ("name", "content"),
    [
        ("reference.fa.fai", "chr2\t100\t6\t60\t61\nchr1\t200\t114\t60\t61\n"),
        ("reference.dict", HEADER),
        ("alignments.sam", HEADER + "read\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"),
        ("alignments.bam", None),
    ],
)
def test_an_order_is_read_from_an_index_a_dictionary_or_an_alignment_header(
    tmp_path: Path, name: str, content: str | None
) -> None:
    """Test that an order is read from a FASTA index, a sequence dictionary, or a SAM or BAM."""
    path = tmp_path / name
    if content is None:
        write_bam_header(path, [("chr2", 100), ("chr1", 200)])
    else:
        _ = path.write_text(content)
    order = ReferenceOrder.from_path(path)
    assert order(Bed3("chr2", start=0, end=1)) == (0, 0, 1)
    assert order(Bed3("chr1", start=0, end=1)) == (1, 0, 1)


@pytest.mark.parametrize("name", ["reference.fa", "reference.fasta.gz"])
def test_a_fasta_is_ordered_by_the_index_or_dictionary_beside_it(tmp_path: Path, name: str) -> None:
    """Test that a FASTA is ordered by its index beside it, or else its dictionary."""
    fasta = tmp_path / name
    _ = fasta.write_text(">chr2\nAC\n>chr1\nGT\n")
    with pytest.raises(ValueError, match="samtools faidx"):
        _ = ReferenceOrder.from_path(fasta)
    dictionary = tmp_path / "reference.dict"
    _ = dictionary.write_text(HEADER)
    assert ReferenceOrder.from_path(fasta)(Bed3("chr1", start=0, end=1)) == (1, 0, 1)
    dictionary.unlink()
    _ = (tmp_path / f"{name}.fai").write_text("chr1\t2\t6\t2\t3\nchr2\t2\t15\t2\t3\n")
    assert ReferenceOrder.from_path(fasta)(Bed3("chr1", start=0, end=1)) == (0, 0, 1)


def test_another_kind_of_file_or_a_file_that_is_not_a_bam_raises(tmp_path: Path) -> None:
    """Test that a file of another kind, or a BAM that is not one, raises naming it."""
    other = tmp_path / "names.txt"
    _ = other.write_text("chr1\n")
    with pytest.raises(ValueError, match=r"names\.txt"):
        _ = ReferenceOrder.from_path(other)
    fake = tmp_path / "fake.bam"
    _ = fake.write_bytes(gzip.compress(b"not a bam"))
    with pytest.raises(ValueError, match=r"fake\.bam is not a BAM file"):
        _ = ReferenceOrder.from_path(fake)
    cut = tmp_path / "cut.bam"
    _ = cut.write_bytes(gzip.compress(b"BAM\x01" + struct.pack("<i", 100) + b"@HD"))
    with pytest.raises(ValueError, match=r"cut\.bam ends within its header"):
        _ = ReferenceOrder.from_path(cut)


def test_a_bam_that_is_not_gzip_or_is_cut_short_raises_naming_it(tmp_path: Path) -> None:
    """Test that a BAM that is not gzip data, or whose compressed data ends early, raises."""
    plain = tmp_path / "plain.bam"
    _ = plain.write_bytes(b"BAM\x01" + struct.pack("<i", 0) + struct.pack("<i", 0))
    with pytest.raises(ValueError, match=r"plain\.bam is not a BAM file"):
        _ = ReferenceOrder.from_path(plain)
    cut = tmp_path / "cut.bam"
    write_bam_header(cut, [("chr1", 100)])
    _ = cut.write_bytes(cut.read_bytes()[:-12])
    with pytest.raises(ValueError, match=r"cut\.bam is not a BAM file"):
        _ = ReferenceOrder.from_path(cut)


@pytest.mark.parametrize("field", ["text", "references"])
def test_a_bam_header_with_a_negative_length_or_count_raises(tmp_path: Path, field: str) -> None:
    """Test that a negative header text length or reference count raises instead of being read."""
    text, references = (-1, 1) if field == "text" else (0, -1)
    corrupt = tmp_path / "corrupt.bam"
    body = b"BAM\x01" + struct.pack("<i", text) + struct.pack("<i", references)
    _ = corrupt.write_bytes(gzip.compress(body))
    with pytest.raises(ValueError, match=r"corrupt\.bam has a corrupt header"):
        _ = ReferenceOrder.from_path(corrupt)
