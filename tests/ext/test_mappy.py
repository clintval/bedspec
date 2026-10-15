"""Tests of the function that reads mappy's records as spans."""

import random
from pathlib import Path

from mappy import Aligner
from mappy import Alignment

from bedspec import Bed6
from bedspec import BedStrand
from bedspec import register_for_bedspec
from bedspec.ext.mappy import alignment_blocks
from bedspec.overlap import TreeDetector


def aligned(tmp_path: Path) -> Alignment:
    """Return a real alignment of 200 bases, a 50-base deletion, then 200 bases, to chr1."""
    rng = random.Random(1)
    reference = "".join(rng.choice("ACGT") for _ in range(2000))
    fasta = tmp_path / "reference.fa"
    _ = fasta.write_text(f">chr1\n{reference}\n")
    query = reference[500:700] + reference[750:950]
    hit: Alignment = next(Aligner(str(fasta), preset="sr").map(query))
    return hit


def test_alignment_blocks_skip_a_deletion_on_the_alignment_strand(tmp_path: Path) -> None:
    """Test that a real alignment's blocks skip its deletion, on its strand."""
    alignment = aligned(tmp_path)
    assert alignment.cigar_str == "200M50D200M"
    assert alignment_blocks(alignment) == [
        Bed6("chr1", start=500, end=700, name=None, score=None, strand=BedStrand.Positive),
        Bed6("chr1", start=750, end=950, name=None, score=None, strand=BedStrand.Positive),
    ]


def test_registered_alignments_are_not_found_in_their_deletions(tmp_path: Path) -> None:
    """Test that alignments registered by their blocks are found by them and not their gaps."""
    alignment = aligned(tmp_path)
    with register_for_bedspec(Alignment, alignment_blocks):
        detector = TreeDetector([alignment])
        assert list(
            detector.overlapping(
                Bed6("chr1", start=800, end=801, name=None, score=None, strand=None)
            )
        ) == [alignment]
        assert not detector.overlaps(
            Bed6("chr1", start=720, end=730, name=None, score=None, strand=None)
        )
