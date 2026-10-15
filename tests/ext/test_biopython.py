"""Tests of the function that reads Biopython's records as spans."""

from pathlib import Path

from Bio import Align
from Bio.Align import Alignment

import bedspec.ext.biopython
from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import register_locatable
from bedspec.ext.biopython import as_blocks
from bedspec.ext.biopython import as_interval
from bedspec.overlap import TreeDetector

BED12 = """chr1\t100\t300\tminus\t0\t-\t100\t300\t0\t2\t50,60,\t0,140,
chr2\t10\t20\tplus\t0\t+\t10\t20\t0\t1\t10,\t0,
"""


def test_as_blocks_of_bed12_lines_on_both_strands(tmp_path: Path) -> None:
    """Test that BED12 lines parsed by Bio.Align give their blocks on their strand."""
    path = tmp_path / "genes.bed"
    _ = path.write_text(BED12)
    minus, plus = list(Align.parse(str(path), "bed"))
    assert as_blocks(minus) == [
        Bed6("chr1", start=100, end=150, name=None, score=None, strand=BedStrand.Negative),
        Bed6("chr1", start=240, end=300, name=None, score=None, strand=BedStrand.Negative),
    ]
    assert as_blocks(plus) == [
        Bed6("chr2", start=10, end=20, name=None, score=None, strand=BedStrand.Positive)
    ]


def test_as_interval_of_bed12_lines_spans_their_blocks_on_both_strands(tmp_path: Path) -> None:
    """Test that a BED12 line's interval runs from its first block to its last, on its strand."""
    path = tmp_path / "genes.bed"
    _ = path.write_text(BED12)
    minus, plus = list(Align.parse(str(path), "bed"))
    assert as_interval(minus) == [
        Bed6("chr1", start=100, end=300, name=None, score=None, strand=BedStrand.Negative)
    ]
    assert as_interval(plus) == [
        Bed6("chr2", start=10, end=20, name=None, score=None, strand=BedStrand.Positive)
    ]


def test_registered_alignments_are_not_found_in_their_introns(tmp_path: Path) -> None:
    """Test that alignments registered by their blocks are found by them and not their introns."""
    path = tmp_path / "genes.bed"
    _ = path.write_text(BED12)
    alignments = list(Align.parse(str(path), "bed"))
    with register_locatable(Alignment, as_blocks):
        detector = TreeDetector(alignments)
        assert list(detector.overlapping(Bed3("chr1", start=250, end=251))) == [alignments[0]]
        assert not detector.overlaps(Bed3("chr1", start=200, end=210))


def test_bedspec_ext_biopython_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.biopython holds its functions and none of the names it imports."""
    module = bedspec.ext.biopython
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert public == sorted(module.__all__) == ["as_blocks", "as_interval"]


def test_as_blocks_joins_blocks_across_an_insertion(tmp_path: Path) -> None:
    """Test that an insertion in a SAM record parsed by Bio.Align does not split its blocks."""
    path = tmp_path / "read.sam"
    _ = path.write_text("@SQ\tSN:chr1\tLN:1000\nread\t0\tchr1\t101\t60\t50M2I50M\t*\t0\t0\t*\t*\n")
    (alignment,) = list(Align.parse(str(path), "sam"))
    assert as_blocks(alignment) == [
        Bed6("chr1", start=100, end=200, name=None, score=None, strand=BedStrand.Positive)
    ]
