from pathlib import Path

import pytest
from pybgzf import IndexFormat

from bedspec import Bed3
from bedspec import BedPE
from bedspec import BedWriter
from bedspec.overlap import TabixDetector


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
    """Test that the detector must be subscripted with a BED type."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    with pytest.raises(TypeError, match=r"TabixDetector\[Bed3\]"):
        with TabixDetector(path):
            pass


def test_a_paired_bed_cannot_be_queried(tmp_path: Path) -> None:
    """Test that BEDPE, which has two intervals, is refused."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    with pytest.raises(TypeError, match="BedPE"):
        with TabixDetector[BedPE](path):  # type: ignore[type-var]  # pyright: ignore[reportInvalidTypeArguments]
            pass
