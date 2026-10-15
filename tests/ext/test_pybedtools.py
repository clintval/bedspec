"""Tests of the function that reads pybedtools' records as spans."""

from pathlib import Path

from pybedtools import BedTool
from pybedtools import Interval
from pybedtools import create_interval_from_list

import bedspec.ext.pybedtools
from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import Territory
from bedspec import register_locatable
from bedspec.ext.pybedtools import interval_span
from bedspec.overlap import TreeDetector


def span(start: int, end: int, strand: BedStrand | None) -> Bed6:
    """Return a span on chr1 on a strand, or on none."""
    return Bed6("chr1", start=start, end=end, name=None, score=None, strand=strand)


def test_interval_span_reads_bed_and_gff_intervals_as_0_based() -> None:
    """Test that BED and GFF intervals are read 0-based, on their strand or on none."""
    bed = create_interval_from_list(["chr1", "10", "20", "name", "0", "-"])
    gff = create_interval_from_list(["chr1", "src", "gene", "11", "20", ".", "+", ".", "ID=g"])
    assert interval_span(bed) == [span(10, 20, BedStrand.Negative)]
    assert interval_span(gff) == [span(10, 20, BedStrand.Positive)]
    assert interval_span(Interval("chr1", 10, 20)) == [span(10, 20, None)]


def test_registered_intervals_are_found_on_their_strand() -> None:
    """Test that intervals registered by their spans are found as themselves, on their strand."""
    minus = create_interval_from_list(["chr1", "10", "20", "name", "0", "-"])
    with register_locatable(Interval, interval_span):
        detector = TreeDetector([minus])
        assert list(detector.overlapping(span(15, 16, BedStrand.Negative), stranded=True)) == [
            minus
        ]
        assert not detector.overlaps(span(15, 16, BedStrand.Positive), stranded=True)
        assert list(Territory([minus])) == [Bed3("chr1", start=10, end=20)]


def test_the_readme_example_finds_the_genes_of_a_gff_file(tmp_path: Path) -> None:
    """Test that the README's lambda registers intervals read from a GFF file."""
    path = tmp_path / "genes.gff"
    _ = path.write_text("chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n")
    with register_locatable(Interval, lambda i: [Bed3(i.chrom, start=i.start, end=i.end)]):
        genes = TreeDetector(BedTool(str(path)))
        assert genes.overlaps(Interval("chr1", 100, 200))
        assert genes.overlaps(Interval("chr1", 199, 250))
        assert not genes.overlaps(Interval("chr1", 200, 300))


def test_bedspec_ext_pybedtools_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.pybedtools holds its functions and none of the names it imports."""
    module = bedspec.ext.pybedtools
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert public == sorted(module.__all__) == ["interval_span"]
