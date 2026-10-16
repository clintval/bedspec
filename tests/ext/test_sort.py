"""Tests of sorting records on disk with spillsort, in a genome's reference order."""

from pathlib import Path

from pybgzf import IndexFormat
from spillsort import Sorter

from bedspec import Bed3
from bedspec import BedReader
from bedspec import BedWriter
from bedspec import ReferenceOrder
from bedspec.overlap import TabixDetector


def test_records_sorted_on_disk_in_reference_order_write_an_indexed_bed(tmp_path: Path) -> None:
    """Test that records spilled to disk come back in reference order and write an indexed BED."""
    order = ReferenceOrder(["chr2", "chr1"])
    records = [
        Bed3(refname, start=start, end=start + 10)
        for start in range(10_000, 0, -7)
        for refname in ("chr1", "chr2")
    ]
    path = tmp_path / "sorted.bed.gz"
    with (
        Sorter[Bed3](order, max_bytes=4_096, tmp_dir=tmp_path) as sorter,
        BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer,
    ):
        writer.write_all(sorter.sort(records))
    with BedReader.from_path[Bed3](path) as reader:
        assert list(reader) == sorted(records, key=order)
    with TabixDetector[Bed3](path) as detector:
        assert list(detector.overlapping(Bed3("chr1", start=0, end=10))) == [
            Bed3("chr1", start=4, end=14)
        ]
    assert sorted(child.name for child in tmp_path.iterdir()) == [
        "sorted.bed.gz",
        "sorted.bed.gz.tbi",
    ]
