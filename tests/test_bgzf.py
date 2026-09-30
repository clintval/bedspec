import bz2
import gzip
from pathlib import Path

import pytest
from pybgzf import IndexedReader
from pybgzf import IndexFormat
from typeline import Comment

from bedspec import Bed2
from bedspec import Bed3
from bedspec import BedPE
from bedspec import BedReader
from bedspec import BedWriter

BEDS = [
    Bed3(refname="chr1", start=1, end=5),
    Bed3(refname="chr1", start=3, end=9),
    Bed3(refname="chr2", start=0, end=2),
]


def is_bgzf(path: Path) -> bool:
    """Return True if a file starts with a BGZF block header."""
    header = path.read_bytes()[:16]
    return header[:4] == b"\x1f\x8b\x08\x04" and header[12:14] == b"BC"


@pytest.mark.parametrize("suffix", [".gz", ".bgz"])
def test_a_compressed_path_is_written_as_bgzf(suffix: str, tmp_path: Path) -> None:
    """Test that a path ending in `.gz` or `.bgz` is written as BGZF and read back."""
    path = tmp_path / f"test.bed{suffix}"
    with BedWriter.from_path[Bed3](path) as writer:
        for bed in BEDS:
            writer.write(bed)

    assert is_bgzf(path)
    assert gzip.decompress(path.read_bytes()) == b"chr1\t1\t5\nchr1\t3\t9\nchr2\t0\t2\n"
    with BedReader.from_path[Bed3](path) as reader:
        assert list(reader) == BEDS


def test_bgzf_is_written_on_several_threads(tmp_path: Path) -> None:
    """Test that BGZF may be compressed on several threads."""
    path = tmp_path / "test.bed.gz"
    with BedWriter.from_path[Bed3](path, threads=4) as writer:
        for bed in BEDS:
            writer.write(bed)

    with BedReader.from_path[Bed3](path) as reader:
        assert list(reader) == BEDS


def test_other_paths_are_written_as_before(tmp_path: Path) -> None:
    """Test that a path without a BGZF suffix is written as plain text, or as its compression."""
    with BedWriter.from_path[Bed3](tmp_path / "test.bed") as writer:
        writer.write(BEDS[0])
    with BedWriter.from_path[Bed3](tmp_path / "test.bed.bz2") as writer:
        writer.write(BEDS[0])

    assert (tmp_path / "test.bed").read_text() == "chr1\t1\t5\n"
    assert bz2.decompress((tmp_path / "test.bed.bz2").read_bytes()) == b"chr1\t1\t5\n"


@pytest.mark.parametrize("index", list(IndexFormat))
def test_an_index_is_written_beside_the_bgzf_file(index: IndexFormat, tmp_path: Path) -> None:
    """Test that a tabix or CSI index is written beside the file when asked for."""
    path = tmp_path / "test.bed.gz"
    with BedWriter.from_path[Bed3](path, index=index) as writer:
        for bed in BEDS:
            writer.write(bed)

    assert (tmp_path / f"test.bed.gz.{index.value}").is_file()
    with IndexedReader(path) as reader:
        assert list(reader.query("chr1", 5, 6)) == ["chr1\t3\t9"]


def test_a_point_bed_is_indexed_as_one_base(tmp_path: Path) -> None:
    """Test that a BED2 record is indexed as the one base it describes."""
    path = tmp_path / "test.bed.gz"
    with BedWriter.from_path[Bed2](path, index=IndexFormat.TBI) as writer:
        writer.write(Bed2(refname="chr1", start=4))

    with IndexedReader(path) as reader:
        assert list(reader.query("chr1", 4, 5)) == ["chr1\t4"]
        assert list(reader.query("chr1", 5, 6)) == []


@pytest.mark.parametrize("path", ["test.bed", "test.bed.bz2"])
def test_an_index_needs_a_bgzf_path(path: str, tmp_path: Path) -> None:
    """Test that an index is refused, before any file is made, unless the path is BGZF."""
    with pytest.raises(ValueError, match=r"\.gz or \.bgz"):
        _ = BedWriter.from_path[Bed3](tmp_path / path, index=IndexFormat.TBI)
    assert list(tmp_path.iterdir()) == []


def test_threads_need_a_bgzf_path(tmp_path: Path) -> None:
    """Test that compression threads are refused, before any file is made, unless BGZF."""
    with pytest.raises(ValueError, match=r"\.gz or \.bgz"):
        _ = BedWriter.from_path[Bed3](tmp_path / "test.bed", threads=2)
    assert list(tmp_path.iterdir()) == []


def test_a_paired_bed_cannot_be_indexed(tmp_path: Path) -> None:
    """Test that an index of BEDPE is refused, before any file is made."""
    with pytest.raises(ValueError, match="BedPE"):
        _ = BedWriter.from_path[BedPE](tmp_path / "test.bedpe.gz", index=IndexFormat.TBI)
    assert list(tmp_path.iterdir()) == []


def test_comments_are_not_indexed(tmp_path: Path) -> None:
    """Test that comments are written to an indexed file and skipped by the index."""
    path = tmp_path / "test.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write_comment("hello")
        writer.write_comment(Comment(1, "#world"))
        writer.write(BEDS[0])

    assert gzip.decompress(path.read_bytes()) == b"# hello\n#world\nchr1\t1\t5\n"
    with IndexedReader(path) as reader:
        assert list(reader.query("chr1", 0, 10)) == ["chr1\t1\t5"]


@pytest.mark.parametrize(
    "comment", ["track name=genes", "browser hide all", "#fine\ntrack", Comment(1, "track")]
)
def test_track_and_browser_lines_are_refused_when_indexing(
    comment: str | Comment, tmp_path: Path
) -> None:
    """Test that a track or browser line is refused, since an index would read it as a feature."""
    with BedWriter.from_path[Bed3](tmp_path / "test.bed.gz", index=IndexFormat.TBI) as writer:
        with pytest.raises(ValueError, match="track or browser"):
            writer.write_comment(comment)
        writer.write(BEDS[0])


def test_track_and_browser_lines_are_written_without_an_index(tmp_path: Path) -> None:
    """Test that track and browser lines are written to BGZF when there is no index."""
    path = tmp_path / "test.bed.gz"
    with BedWriter.from_path[Bed3](path) as writer:
        writer.write_comment("track name=genes")
        writer.write(BEDS[0])

    assert gzip.decompress(path.read_bytes()) == b"track name=genes\nchr1\t1\t5\n"
