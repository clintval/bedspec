import inspect
import io
from dataclasses import asdict as as_dict
from pathlib import Path
from types import FrameType
from types import TracebackType
from typing import ContextManager
from typing import Generic
from typing import Iterable
from typing import Iterator
from typing import TypeVar
from typing import _GenericAlias  # type: ignore[attr-defined]
from typing import cast

import bedspec._typing

####################################################################################################

BED_EXTENSION: str = ".bed"
"""The default file extension for BED files."""

BEDGRAPH_EXTENSION: str = ".bedgraph"
"""The default file extension for bedGraph files."""

BEDPE_EXTENSION: str = ".bedpe"
"""The default file extension for BedPE files."""

TRACK_EXTENSION: str = ".track"
"""The default file extension for UCSC track files."""

####################################################################################################

COMMENT_PREFIXES: set[str] = {"#", "browser", "track"}
"""The set of BED comment prefixes that this library supports."""

MISSING_FIELD: str = "."
"""The string used to indicate a missing field in a BED record."""

####################################################################################################

BedKind = TypeVar("BedKind", bound=bedspec._types.BedType)
"""A type variable for any kind of BED record type."""


class BedWriter(Generic[BedKind], ContextManager):
    """A writer of BED records.

    Args:
        handle: An open file-like object to write to.

    Attributes:
        bed_kind: The kind of BED type that this writer will write.

    """

    bed_kind: type[BedKind] | None

    def __class_getitem__(cls, key: object) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return _GenericAlias(cls, key)  # type: ignore[no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedWriter[BedKind]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        typelevel = signature.f_locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __enter__(self) -> "BedWriter[BedKind]":
        """Enter this context."""
        return self

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED writer without knowing yet what BED types we will write."""
        self._handle = handle

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return super().__exit__(__exc_type, __exc_value, __traceback)

    @bedspec._typing.classmethod_generic
    def from_path(cls, path: Path | str) -> "BedWriter[BedKind]":
        """Open a BED writer from a file path."""
        reader = cls(handle=Path(path).open("w"))  # type: ignore[operator]
        reader.bed_kind = None if len(cls.__args__) == 0 else cls.__args__[0]  # type: ignore[attr-defined]
        return cast("BedWriter[BedKind]", reader)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()

    def write_comment(self, comment: str) -> None:
        """Write a comment to the BED output."""
        for line in comment.splitlines():
            if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
                self._handle.write(f"{comment}\n")
            else:
                self._handle.write(f"# {comment}\n")

    def write(self, bed: BedKind) -> None:
        """Write a BED record to the BED output."""
        if self.bed_kind is not None:
            if type(bed) is not self.bed_kind:
                raise TypeError(
                    f"BedWriter can only continue to write features of the same type."
                    f" Will not write a {type(bed).__name__} after a {self.bed_kind.__name__}."
                )
        else:
            self.bed_kind = type(bed)

        self._handle.write(f"{"\t".join(map(str, as_dict(bed).values()))}\n")

    def write_all(self, beds: Iterable[BedKind]) -> None:
        """Write all the BED records to the BED output."""
        for bed in beds:
            self.write(bed)


class BedReader(Generic[BedKind], ContextManager, Iterable[BedKind]):
    """A reader of BED records.

    This reader is capable of reading BED records but must be typed at runtime:

    ```python
    from bedspec import BedReader, Bed3

    with BedReader[Bed3](path) as reader:
        print(list(reader)
    ```

    Args:
        handle: An open file-like object to read from.

    Attributes:
        bed_kind: The kind of BED type that this reader will read.

    """

    bed_kind: type[BedKind] | None

    def __class_getitem__(cls, key: object) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return _GenericAlias(cls, key)  # type: ignore[no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedReader[BedKind]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        typelevel = signature.f_locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED reader without knowing yet what BED types we will write."""
        self._handle = handle

    def __enter__(self) -> "BedReader[BedKind]":
        """Enter this context."""
        return self

    def __iter__(self) -> Iterator[BedKind]:
        """Iterate through the BED records of this IO handle."""
        # TODO: Implement __next__ and type this class as an iterator.
        if self.bed_kind is None:
            raise NotImplementedError("Untyped reading is not yet supported!")
        for line in map(lambda line: line.strip(), self._handle):
            if line == "" or any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
                continue
            yield cast(BedKind, self.bed_kind.decode(line))

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return super().__exit__(__exc_type, __exc_value, __traceback)

    @bedspec._typing.classmethod_generic
    def from_path(cls, path: Path | str) -> "BedReader[BedKind]":
        """Open a BED reader from a file path."""
        reader = cls(handle=Path(path).open("r"))  # type: ignore[operator]
        reader.bed_kind = None if len(cls.__args__) == 0 else cls.__args__[0]  # type: ignore[attr-defined]
        return cast("BedReader[BedKind]", reader)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()
