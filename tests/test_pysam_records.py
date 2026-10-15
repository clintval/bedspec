"""Tests that features naming their reference, start and end like pysam's records work unwrapped."""

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import pysam
import pytest

from bedspec import Bed3
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import Territory
from bedspec import locate
from bedspec.overlap import TabixDetector
from bedspec.overlap import TreeDetector


@dataclass
class ChromRecord:
    """A record naming its reference `chrom`, as pybedlite's and pybedtools' records do."""

    chrom: str
    start: int
    end: int


def variant(start: int, ref: str = "C", contig: str = "chr1") -> pysam.VariantRecord:
    """Return a variant record at a 0-based start, spanning its reference allele."""
    header = pysam.VariantHeader()
    header.add_line("##contig=<ID=chr1,length=1000>")
    header.add_line("##contig=<ID=chr2,length=1000>")
    return header.new_record(contig=contig, start=start, alleles=(ref, "T"))  # pyright: ignore[reportUnknownMemberType]


def alignment(start: int, length: int, mapped: bool = True) -> pysam.AlignedSegment:
    """Return a read aligned to chr1 at a 0-based start with a match of the given length."""
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})  # pyright: ignore[reportUnknownMemberType]
    read = pysam.AlignedSegment(header)
    read.query_name = "read"
    read.query_sequence = "A" * length
    if mapped:
        read.reference_name = "chr1"
        read.reference_start = start
        read.cigarstring = f"{length}M"
    else:
        read.is_unmapped = True
    return read


def gtf_record(start: int, end: int) -> Any:
    """Return a pysam tabix GTF record on chr1 for a 0-based half-open interval."""
    line = f'chr1\ttest\texon\t{start + 1}\t{end}\t.\t+\t.\tgene_id "g"; transcript_id "t";'
    handle = StringIO(f"{line}\n")
    records = pysam.tabix_generic_iterator(handle, parser=pysam.asGTF())  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]
    return next(iter(records))


@pytest.mark.parametrize(
    ("feature", "expected"),
    [
        pytest.param(Bed3("chr1", start=5, end=9), ("chr1", 5, 9), id="bed3"),
        pytest.param(variant(99, ref="CAG"), ("chr1", 99, 102), id="pysam-variant"),
        pytest.param(alignment(100, 50), ("chr1", 100, 150), id="pysam-alignment"),
        pytest.param(gtf_record(10, 20), ("chr1", 10, 20), id="pysam-tabix-gtf"),
        pytest.param(ChromRecord("chr2", 3, 7), ("chr2", 3, 7), id="chrom-dataclass"),
    ],
)
def test_locate_reads_the_reference_start_and_end_under_any_known_name(
    feature: Any, expected: tuple[str, int, int]
) -> None:
    """Test that locate() reads a feature's reference, start and end whatever names it uses."""
    assert locate(feature)(feature) == expected


def test_locate_refuses_a_feature_with_none_of_the_names() -> None:
    """Test that a feature without a reference name raises an error naming the ones looked for."""
    with pytest.raises(AttributeError, match="has none of the attributes refname, contig, chrom"):
        _ = locate(object())


def test_an_unmapped_read_is_refused() -> None:
    """Test that a read placed on no reference is refused rather than indexed under None."""
    read = alignment(0, 10, mapped=False)
    with pytest.raises(ValueError, match="must be placed on a reference"):
        _ = TreeDetector([read])
    with pytest.raises(ValueError, match="must be placed on a reference"):
        _ = TreeDetector([Bed3("chr1", start=0, end=10)]).overlaps(read)


def test_a_tree_detector_holds_and_returns_pysam_records() -> None:
    """Test that pysam records are added and returned as they are, without wrapping."""
    snv, deletion, elsewhere = variant(10), variant(20, ref="CAGT"), variant(10, contig="chr2")
    detector = TreeDetector([snv, deletion, elsewhere])
    assert list(detector.overlapping(Bed3("chr1", start=0, end=15))) == [snv]
    assert detector.overlaps(Bed3("chr1", start=22, end=23))
    assert not detector.overlaps(Bed3("chr1", start=24, end=30))
    assert list(detector.enclosing(Bed3("chr1", start=21, end=23))) == [deletion]
    assert list(detector.enclosed_by(Bed3("chr1", start=5, end=30))) in (
        [snv, deletion],
        [deletion, snv],
    )
    assert set(map(id, detector)) == {id(snv), id(deletion), id(elsewhere)}


def test_pysam_records_query_a_detector_of_bed_records() -> None:
    """Test that pysam records, as queries, find the BED records they overlap."""
    target = Bed3("chr1", start=100, end=200)
    detector = TreeDetector([target, Bed3("chr1", start=300, end=400)])
    assert list(detector.overlapping(variant(150))) == [target]
    assert list(detector.overlapping(alignment(180, 50))) == [target]
    assert list(detector.enclosing(alignment(120, 30))) == [target]
    assert not detector.overlaps(variant(250))


def test_a_detector_holding_bed12_spans_and_pysam_records_finds_both() -> None:
    """Test that pysam records sit beside multi-span BED records in one detector."""
    transcript = Bed12(
        "chr1",
        start=100,
        end=200,
        name="tx",
        score=None,
        strand=BedStrand.Positive,
        thick_start=None,
        thick_end=None,
        item_rgb=None,
        block_count=2,
        block_sizes=(10, 20),
        block_starts=(0, 80),
    )
    read = alignment(105, 10)
    detector: TreeDetector[Bed12 | pysam.AlignedSegment] = TreeDetector([transcript, read])
    assert {id(hit) for hit in detector.overlapping(variant(106))} == {id(transcript), id(read)}
    assert list(detector.overlapping(variant(150))) == []


def test_a_stranded_query_by_a_pysam_record_finds_nothing() -> None:
    """Test that a pysam record has no strand, so a stranded query by it matches nothing."""
    detector = TreeDetector([
        Bed6("chr1", start=0, end=100, name=None, score=None, strand=BedStrand.Positive)
    ])
    assert detector.overlaps(variant(10))
    assert not detector.overlaps(variant(10), stranded=True)


def test_a_territory_of_pysam_records() -> None:
    """Test that a territory is built from pysam records and other libraries' records."""
    territory = Territory([variant(10, ref="CAG"), alignment(11, 10), ChromRecord("chr2", 0, 5)])
    assert list(territory) == [Bed3("chr1", start=10, end=21), Bed3("chr2", start=0, end=5)]


def test_a_tabix_detector_answers_a_pysam_record_query(tmp_path: Path) -> None:
    """Test that a query by a pysam record reads an indexed BED file like any other query."""
    from pybgzf import IndexFormat

    path = tmp_path / "targets.bed.gz"
    target = Bed3("chr1", start=100, end=200)
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(target)
        writer.write(Bed3("chr1", start=300, end=400))
    with TabixDetector[Bed3](path) as detector:
        assert list(detector.overlapping(variant(150))) == [target]
        assert list(detector.overlapping(alignment(190, 20))) == [target]
        assert not detector.overlaps(variant(250))
