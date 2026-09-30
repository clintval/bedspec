from collections.abc import Iterator
from contextlib import AbstractContextManager
from functools import cached_property
from io import StringIO
from pathlib import Path
from types import TracebackType
from typing import Any
from typing import Generic
from typing import TypeVar
from typing import cast
from typing import get_args

from pybgzf import IndexedReader
from typing_extensions import Self
from typing_extensions import override

from bedspec._bedspec import PointBed
from bedspec._bedspec import ReferenceSpan
from bedspec._bedspec import SimpleBed
from bedspec._reader import BedReader
from bedspec.overlap._overlap import _closed  # pyright: ignore[reportPrivateUsage]
from bedspec.overlap._overlap import _strand  # pyright: ignore[reportPrivateUsage]

IntervalBedType = TypeVar("IntervalBedType", bound=PointBed | SimpleBed)
"""A type variable for a BED record type that describes one interval."""


class IndexedOverlapDetector(
    AbstractContextManager["IndexedOverlapDetector[IntervalBedType]"], Generic[IntervalBedType]
):
    """Detects overlaps with the features of a BGZF BED file, read through its tabix or CSI index.

    Queries are answered exactly as `OverlapDetector` answers them, reading only the parts of the
    file the index points to, and features are read as the record type the detector is
    subscripted with, e.g. `IndexedOverlapDetector[Bed6](path)`.
    A zero-length feature at the start of a reference is never found, since tabix never returns
    it.

    Args:
        path: the BGZF BED file.
        index_path: its index; defaults to `path` plus `.csi` or, if there is none, `.tbi`.
        threads: the number of threads decompressing the file.
    """

    def __init__(
        self, path: Path | str, *, index_path: Path | str | None = None, threads: int = 1
    ) -> None:
        self._reader: IndexedReader = IndexedReader(path, index_path=index_path, threads=threads)

    @cached_property
    def _record_type(self) -> type[IntervalBedType]:
        """The record type this detector is subscripted with, refusing any other kind of BED."""
        name = type(self).__name__
        args: tuple[Any, ...] = get_args(getattr(self, "__orig_class__", None))
        if not args:
            raise TypeError(f"{name} must be subscripted with a BED type, e.g. {name}[Bed3]!")
        record_type: Any = args[0]
        if not isinstance(record_type, type) or not issubclass(record_type, PointBed | SimpleBed):
            raise TypeError(f"{name} needs a BED type of one interval each, not {record_type}!")
        return cast("type[IntervalBedType]", record_type)

    @property
    def closed(self) -> bool:
        """True once the file is closed."""
        return self._reader.closed

    @override
    def __enter__(self) -> Self:
        """Enter this context, checking the record type."""
        _ = super().__enter__()
        _ = self._record_type
        return self

    @override
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Exit this context, closing the file."""
        self.close()

    def close(self) -> None:
        """Close the file."""
        self._reader.close()

    def overlapping(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[IntervalBedType]:
        """Yields all the overlapping features for a given query feature."""
        strand = _strand(feature)
        if stranded and strand is None:
            return
        start, end = _closed(feature)
        lines = self._reader.query(feature.refname, max(start - 1, 0), end + 2)
        text = StringIO("".join(f"{line}\n" for line in lines))
        for record in BedReader[self._record_type](text):  # type: ignore[name-defined]
            span = next(record.territory())
            found_start, found_end = _closed(span)
            if found_start > end or found_end < start:
                continue
            if stranded and _strand(record) is not strand:
                continue
            yield record

    def overlaps(self, feature: ReferenceSpan, *, stranded: bool = False) -> bool:
        """Determine if a query feature overlaps any other features."""
        return next(self.overlapping(feature, stranded=stranded), None) is not None

    def enclosing(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[IntervalBedType]:
        """Yields all the overlapping features that completely enclose the given query feature."""
        for overlap in self.overlapping(feature, stranded=stranded):
            span = next(overlap.territory())
            if feature.start >= span.start and feature.end <= span.end:
                yield overlap

    def enclosed_by(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[IntervalBedType]:
        """Yields all the overlapping features that are enclosed by the given query feature."""
        for overlap in self.overlapping(feature, stranded=stranded):
            span = next(overlap.territory())
            if feature.start <= span.start and feature.end >= span.end:
                yield overlap
