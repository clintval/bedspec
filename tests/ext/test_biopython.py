"""Tests of the function that reads Biopython's records as spans."""

from pathlib import Path

from Bio import Align
from Bio.Align import Alignment

import bedspec.ext.biopython
from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import register_for_bedspec
from bedspec.ext.biopython import alignment_blocks
from bedspec.overlap import TreeDetector

BED12 = """chr1\t100\t300\tminus\t0\t-\t100\t300\t0\t2\t50,60,\t0,140,
chr2\t10\t20\tplus\t0\t+\t10\t20\t0\t1\t10,\t0,
"""


def test_alignment_blocks_of_bed12_lines_on_both_strands(tmp_path: Path) -> None:
    """Test that BED12 lines parsed by Bio.Align give their blocks on their strand."""
    path = tmp_path / "genes.bed"
    _ = path.write_text(BED12)
    minus, plus = list(Align.parse(str(path), "bed"))
    assert alignment_blocks(minus) == [
        Bed6("chr1", start=100, end=150, name=None, score=None, strand=BedStrand.Negative),
        Bed6("chr1", start=240, end=300, name=None, score=None, strand=BedStrand.Negative),
    ]
    assert alignment_blocks(plus) == [
        Bed6("chr2", start=10, end=20, name=None, score=None, strand=BedStrand.Positive)
    ]


def test_registered_alignments_are_not_found_in_their_introns(tmp_path: Path) -> None:
    """Test that alignments registered by their blocks are found by them and not their introns."""
    path = tmp_path / "genes.bed"
    _ = path.write_text(BED12)
    alignments = list(Align.parse(str(path), "bed"))
    with register_for_bedspec(Alignment, alignment_blocks):
        detector = TreeDetector(alignments)
        assert list(detector.overlapping(Bed3("chr1", start=250, end=251))) == [alignments[0]]
        assert not detector.overlaps(Bed3("chr1", start=200, end=210))


def test_bedspec_ext_biopython_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.biopython holds its functions and none of the names it imports."""
    module = bedspec.ext.biopython
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert public == sorted(module.__all__) == ["alignment_blocks"]
