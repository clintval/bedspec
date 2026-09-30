from inspect import unwrap
from io import StringIO
from pathlib import Path
from typing import TextIO

import pybgzf
from pybgzf import Columns
from pybgzf import IndexFormat
from typeline import Comment
from typeline import SubscriptableClassmethod
from typeline import TsvWriter
from typeline import WriterOptions
from typing_extensions import Self
from typing_extensions import Unpack
from typing_extensions import override

from bedspec._bedspec import COMMENT_PREFIXES
from bedspec._bedspec import MISSING_FIELD
from bedspec._bedspec import BedType
from bedspec._bedspec import PointBed
from bedspec._bedspec import SimpleBed
from bedspec._codecs import BED_CODECS

BED_COMMENT_PREFIXES: tuple[str, ...] = ("#", *sorted(COMMENT_PREFIXES - {"#"}))
"""The BED comment prefixes, with the one added to comment lines that have none first."""

BGZF_SUFFIXES: tuple[str, ...] = (".bgz", ".gz")
"""The file extensions written as BGZF, which is also valid gzip."""


class BedWriter(TsvWriter[BedType]):
    """A writer for writing dataclasses into BED text data."""

    @override
    def __init__(self, handle: TextIO, /, **options: Unpack[WriterOptions]) -> None:
        """Instantiate a new BED writer.

        Args:
            handle: a file-like object to write delimited data to.
            options: the options of the writer, with BED defaults for any not given.
        """
        _ = options.setdefault("codecs", BED_CODECS)
        _ = options.setdefault("comment_prefixes", BED_COMMENT_PREFIXES)
        _ = options.setdefault("none_field", MISSING_FIELD)
        _ = options.setdefault("quoting", False)
        super().__init__(handle, **options)
        self._indexed: bool = False

    @override
    def write_comment(self, comment: str | Comment) -> None:
        """Write a comment, refusing a track or browser line when the file is being indexed.

        An index skips only lines starting with `#`, so it would read any other line as a feature.
        """
        if self._indexed:
            prefixes = self._comment_prefixes
            lines = (
                [comment.text]
                if isinstance(comment, Comment)
                else [
                    line if line.startswith(prefixes) else f"{prefixes[0]} {line}"
                    for line in comment.splitlines()
                ]
            )
            for line in lines:
                if not line.startswith("#"):
                    raise ValueError(
                        "Cannot write a track or browser line to an indexed file, which would"
                        + f" index it as a feature! Start it with '#' instead: {line!r}"
                    )
        super().write_comment(comment)

    @SubscriptableClassmethod
    @classmethod
    def from_path(  # pyright: ignore[reportIncompatibleVariableOverride]
        cls,
        path: Path | str,
        /,
        *,
        index: IndexFormat | None = None,
        threads: int = 1,
        **options: Unpack[WriterOptions],
    ) -> Self:
        """Construct a BED writer from a file path.

        A path ending in `.gz` or `.bgz` is written as BGZF, which any gzip reader can read, and
        can be indexed with tabix or CSI as it is written, on as many threads as given.
        Features must then be sorted by start within each reference, and each reference must be
        contiguous.
        Other paths are written as UTF-8, and compressed when they end in `.bz2` or `.xz`.
        The writer is checked before the file is opened, so a refused writer leaves a file alone.

        Args:
            path: the path to the file to write BED to.
            index: the kind of index to write beside a BGZF file, or None to write none.
            threads: the number of threads compressing a BGZF file.
            options: the options of the writer, with BED defaults for any not given.
        """
        path = Path(path).expanduser()
        if path.suffix not in BGZF_SUFFIXES:
            if index is not None or threads != 1:
                raise ValueError(
                    f"An index and threads need a BGZF path ending in .gz or .bgz, not: {path}"
                )
            plain: Self = unwrap(super().from_path)(cls, path, **options)
            return plain
        _ = cls(StringIO(), **options)
        columns = cls._index_columns() if index is not None else None
        handle = pybgzf.open_writer(path, columns=columns, index=index, newline="", threads=threads)
        try:
            writer = cls(handle, **options)
        except BaseException:
            handle.close()
            raise
        writer._indexed = index is not None
        return writer

    @classmethod
    def _index_columns(cls) -> Columns:
        """Return the columns an index reads each feature's location from, refusing BEDPE."""
        record_type = cls._parameterized_record_type
        if record_type is not None and issubclass(record_type, PointBed):
            return Columns.BED2
        if record_type is not None and issubclass(record_type, SimpleBed):
            return Columns.BED
        name = "records" if record_type is None else record_type.__name__
        raise ValueError(f"Cannot index {name} records, which do not each describe one interval!")
