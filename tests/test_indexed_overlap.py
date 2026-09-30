from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pytest
from pybgzf import IndexFormat
from typing_extensions import override

from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import ReferenceSpan
from bedspec.overlap import TabixDetector


@dataclass(frozen=True)
class Blocked(Bed3):
    """A BED3 whose territory is only its first and last ten bases."""

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield the first and last ten bases of this feature."""
        yield Bed3(self.refname, start=self.start, end=self.start + 10)
        yield Bed3(self.refname, start=self.end - 10, end=self.end)


@dataclass(frozen=True)
class Blocked6(Bed6):
    """A BED6 whose territory is only its first and last ten bases, without a strand."""

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield the first and last ten bases of this feature."""
        yield Bed3(self.refname, start=self.start, end=self.start + 10)
        yield Bed3(self.refname, start=self.end - 10, end=self.end)


def write_sorted(path: Path, records: list[Bed3], index: IndexFormat) -> None:
    """Write sorted BED3 records to an indexed BGZF file."""
    with BedWriter.from_path[Bed3](path, index=index) as writer:
        for record in records:
            writer.write(record)


def test_zero_length_features_overlap_the_bases_beside_them(tmp_path: Path) -> None:
    """Test that a zero-length feature is found by queries of either base beside it."""
    insertion = Bed3("chr1", start=5, end=5)
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [insertion], IndexFormat.TBI)

    with TabixDetector[Bed3](path) as detector:
        assert list(detector.overlapping(Bed3("chr1", start=4, end=5))) == [insertion]
        assert list(detector.overlapping(Bed3("chr1", start=5, end=6))) == [insertion]
        assert list(detector.overlapping(Bed3("chr1", start=5, end=5))) == [insertion]
        assert list(detector.overlapping(Bed3("chr1", start=3, end=4))) == []
        assert list(detector.overlapping(Bed3("chr1", start=6, end=7))) == []


def test_a_zero_length_feature_at_the_start_of_a_reference_is_never_found(tmp_path: Path) -> None:
    """Test that a zero-length feature at position 0 is not found, as tabix never returns it."""
    path = tmp_path / "features.bed.gz"
    with pytest.warns(UserWarning, match="never returned by an index query"):
        write_sorted(path, [Bed3("chr1", start=0, end=0)], IndexFormat.TBI)

    with TabixDetector[Bed3](path) as detector:
        assert not detector.overlaps(Bed3("chr1", start=0, end=1))


def test_an_index_path_and_threads_may_be_given(tmp_path: Path) -> None:
    """Test that the index may be found at any path and read on several threads."""
    bed = Bed3("chr1", start=1, end=5)
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [bed], IndexFormat.CSI)
    index_path = (tmp_path / "features.bed.gz.csi").rename(tmp_path / "elsewhere.csi")

    with TabixDetector[Bed3](path, index_path=index_path, threads=2) as detector:
        assert list(detector.overlapping(bed)) == [bed]


def test_the_file_is_closed_when_the_context_ends(tmp_path: Path) -> None:
    """Test that the file is closed when the context ends."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    with TabixDetector[Bed3](path) as detector:
        assert not detector.closed
    assert detector.closed


def test_a_record_type_is_required(tmp_path: Path) -> None:
    """Test that the detector must be subscripted with a BED type, and closes its file if not."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    detector: TabixDetector[Bed3] = TabixDetector(path)
    with pytest.raises(TypeError, match=r"TabixDetector\[Bed3\]"), detector:
        pass
    assert detector.closed


def test_a_paired_bed_cannot_be_queried(tmp_path: Path) -> None:
    """Test that BEDPE, which has two intervals, is refused, and the file is closed."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    detector = TabixDetector[BedPE](path)  # type: ignore[type-var]  # pyright: ignore[reportInvalidTypeArguments]
    with pytest.raises(TypeError, match="BedPE"), detector:
        pass
    assert detector.closed


def test_every_span_of_a_territory_is_found_and_the_gap_between_them_is_not(
    tmp_path: Path,
) -> None:
    """Test that a query of any span of a feature's territory finds it, and a gap does not."""
    blocked = Blocked("chr1", start=0, end=100)
    path = tmp_path / "features.bed.gz"
    with BedWriter.from_path[Blocked](path, index=IndexFormat.TBI) as writer:
        writer.write(blocked)

    with TabixDetector[Blocked](path) as detector:
        assert list(detector.overlapping(Bed3("chr1", start=5, end=6))) == [blocked]
        assert list(detector.overlapping(Bed3("chr1", start=95, end=96))) == [blocked]
        assert list(detector.overlapping(Bed3("chr1", start=0, end=100))) == [blocked]
        assert list(detector.overlapping(Bed3("chr1", start=40, end=60))) == []
        assert list(detector.enclosing(Bed3("chr1", start=92, end=98))) == [blocked]
        assert list(detector.enclosing(Bed3("chr1", start=5, end=95))) == []


def test_a_feature_is_enclosed_only_when_every_span_of_its_territory_is(tmp_path: Path) -> None:
    """Test that enclosed_by needs every span of a feature's territory inside the query."""
    blocked = Blocked("chr1", start=0, end=100)
    path = tmp_path / "features.bed.gz"
    with BedWriter.from_path[Blocked](path, index=IndexFormat.TBI) as writer:
        writer.write(blocked)

    with TabixDetector[Blocked](path) as detector:
        assert list(detector.enclosed_by(Bed3("chr1", start=0, end=20))) == []
        assert list(detector.enclosed_by(Bed3("chr1", start=80, end=100))) == []
        assert list(detector.enclosed_by(Bed3("chr1", start=0, end=100))) == [blocked]


def test_a_span_without_a_strand_takes_the_strand_of_its_feature(tmp_path: Path) -> None:
    """Test that a stranded query compares a strandless span with its feature's strand."""
    blocked = Blocked6("chr1", start=0, end=100, name="a", score=0, strand=BedStrand.Positive)
    path = tmp_path / "features.bed.gz"
    with BedWriter.from_path[Blocked6](path, index=IndexFormat.TBI) as writer:
        writer.write(blocked)

    same = Bed6("chr1", start=95, end=96, name="q", score=0, strand=BedStrand.Positive)
    other = Bed6("chr1", start=95, end=96, name="q", score=0, strand=BedStrand.Negative)
    with TabixDetector[Blocked6](path) as detector:
        assert list(detector.overlapping(same, stranded=True)) == [blocked]
        assert list(detector.overlapping(other, stranded=True)) == []
        whole = Bed6("chr1", start=0, end=100, name="q", score=0, strand=BedStrand.Positive)
        assert list(detector.enclosed_by(whole, stranded=True)) == [blocked]
