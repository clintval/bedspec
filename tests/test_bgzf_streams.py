import gzip
import os
import sys
import threading
from collections.abc import Callable
from io import StringIO
from pathlib import Path
from typing import Any
from typing import BinaryIO

import pybgzf
import pytest
from pybgzf import IndexFormat
from typeline import TsvWriter

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed6
from bedspec import Bed12N
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import NarrowPeak
from bedspec.overlap import TabixDetector

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="pipes and FIFOs need POSIX")

TIMEOUT = 30.0
"""The seconds to wait for the other end of a pipe before failing instead of hanging."""

RECORDS: dict[type[Any], list[Any]] = {
    Bed2: [Bed2("chr1", start=start) for start in range(0, 20_000, 3)],
    Bed6: [
        Bed6("chr1", start=start, end=start + 7, name="f", score=5, strand=BedStrand.Negative)
        for start in range(0, 20_000, 3)
    ],
    Bed12N: [
        Bed12N(
            "chr2",
            start=start,
            end=start + 4,
            name=None,
            score=None,
            strand=None,
            thick_start=None,
            thick_end=None,
            item_rgb=None,
            block_count=2,
            block_sizes=(1, 1),
            block_starts=(0, 3),
            extra=("x",),
        )
        for start in range(0, 20_000, 3)
    ],
    NarrowPeak: [
        NarrowPeak(
            "chr1",
            start=start,
            end=start,
            name=None,
            score=None,
            strand=BedStrand.Positive,
            signal_value=1.5,
            p_value=-1.0,
            q_value=2.0,
            peak=-1,
        )
        for start in range(1, 20_000, 3)
    ],
}
"""Enough sorted records of a few BED types to fill several BGZF blocks."""

RECORD_TYPES = pytest.mark.parametrize("record_type", list(RECORDS), ids=lambda t: t.__name__)
THREADS = pytest.mark.parametrize("threads", [1, 4])


class Background:
    """Run a function in a thread, and return its result or raise its error when joined."""

    def __init__(self, function: Callable[[], Any]) -> None:
        """Start running a function in a thread."""
        self._result: list[Any] = []
        self._error: list[BaseException] = []

        def run() -> None:
            try:
                self._result.append(function())
            except BaseException as error:
                self._error.append(error)

        self._thread: threading.Thread = threading.Thread(target=run, daemon=True)
        self._thread.start()

    def join(self) -> Any:
        """Wait for the function, failing if it takes longer than `TIMEOUT`."""
        self._thread.join(TIMEOUT)
        assert not self._thread.is_alive(), "the other end of the pipe never finished"
        if self._error:
            raise self._error[0]
        return self._result[0]


def read_all(source: BinaryIO | Path) -> Callable[[], bytes]:
    """Return a function that reads every byte from a binary file or a path, then closes it."""

    def read() -> bytes:
        with source if not isinstance(source, Path) else source.open("rb") as handle:
            return handle.read()

    return read


def parse(record_type: type[Any], data: bytes) -> list[Any]:
    """Decompress BGZF bytes and read the BED records in them."""
    handle = StringIO(gzip.decompress(data).decode())
    return list(BedReader[record_type](handle))  # type: ignore[valid-type]


def write_all(writer: TsvWriter[Any], records: list[Any]) -> None:
    """Write a comment and then every record."""
    writer.write_comment("streamed")
    for record in records:
        writer.write(record)


@RECORD_TYPES
@THREADS
def test_bgzf_is_written_into_a_pipe(record_type: type[Any], threads: int) -> None:
    """Test that BGZF BED is written into a pipe while it is read from the other end."""
    read_end, write_end = os.pipe()
    reader = Background(read_all(os.fdopen(read_end, "rb")))
    with os.fdopen(write_end, "wb") as sink:
        handle = pybgzf.open_writer(sink, newline="", threads=threads)
        with BedWriter[record_type](handle) as writer:  # type: ignore[valid-type]
            write_all(writer, RECORDS[record_type])

    assert parse(record_type, reader.join()) == RECORDS[record_type]


@RECORD_TYPES
@THREADS
def test_bgzf_is_written_into_a_fifo_by_path(
    record_type: type[Any], threads: int, tmp_path: Path
) -> None:
    """Test that BGZF BED is written to a FIFO named like a BGZF file while it is read."""
    fifo = tmp_path / "stream.bed.gz"
    os.mkfifo(fifo)
    reader = Background(read_all(fifo))
    with BedWriter.from_path[record_type](fifo, threads=threads) as writer:
        write_all(writer, RECORDS[record_type])

    assert parse(record_type, reader.join()) == RECORDS[record_type]


@THREADS
@pytest.mark.parametrize("index", [IndexFormat.CSI, IndexFormat.TBI], ids=["CSI", "TBI"])
def test_a_fifo_is_indexed_to_an_index_path(
    index: IndexFormat, threads: int, tmp_path: Path
) -> None:
    """Test that BGZF streamed into a FIFO is indexed to a given path, and the index queries it."""
    fifo = tmp_path / "stream.bed.gz"
    os.mkfifo(fifo)
    index_path = tmp_path / "saved.bed.gz.index"
    reader = Background(read_all(fifo))
    with BedWriter.from_path[Bed6](
        fifo, index=index, index_path=index_path, threads=threads
    ) as writer:
        write_all(writer, RECORDS[Bed6])

    saved = tmp_path / "saved.bed.gz"
    _ = saved.write_bytes(reader.join())
    query = Bed3("chr1", start=100, end=103)
    with TabixDetector[Bed6](saved, index_path=index_path) as detector:
        assert list(detector.overlapping(query)) == [
            record for record in RECORDS[Bed6] if record.start < 103 and record.end > 100
        ]


def test_a_fifo_is_not_indexed_without_an_index_path(tmp_path: Path) -> None:
    """Test that indexing into a FIFO without an index path is refused before the FIFO is opened."""
    fifo = tmp_path / "stream.bed.gz"
    os.mkfifo(fifo)
    read_end = os.open(fifo, os.O_RDONLY | os.O_NONBLOCK)
    try:
        with pytest.raises(ValueError, match="index_path"):
            _ = BedWriter.from_path[Bed6](fifo, index=IndexFormat.TBI)
        assert os.read(read_end, 1) == b""
    finally:
        os.close(read_end)
    assert list(tmp_path.iterdir()) == [fifo]


@RECORD_TYPES
@THREADS
def test_bgzf_is_read_from_a_pipe(record_type: type[Any], threads: int, tmp_path: Path) -> None:
    """Test that BGZF BED is read from a pipe while it is written to the other end."""
    path = tmp_path / "features.bed.gz"
    with BedWriter.from_path[record_type](path) as writer:
        write_all(writer, RECORDS[record_type])
    read_end, write_end = os.pipe()

    def write() -> None:
        with os.fdopen(write_end, "wb") as sink:
            _ = sink.write(path.read_bytes())

    writer_thread = Background(write)
    with (
        os.fdopen(read_end, "rb") as source,
        BedReader[record_type](pybgzf.open_reader(source, threads=threads)) as reader,  # type: ignore[valid-type]
    ):
        assert list(reader) == RECORDS[record_type]
    writer_thread.join()
