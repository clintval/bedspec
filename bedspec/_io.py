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

BGZ_EXTENSION: str = ".bgz"
"""The default file extension for block-compressed gzip files (`.bgz`)."""

BGZIP_EXTENSION: str = ".bgzip"
"""The default file extension for block-compressed gzip files (`.bgzip`)."""

GZ_EXTENSION: str = ".gz"
"""The default file extension for compressed gzip files (`.gz`)."""

GZIP_EXTENSION: str = ".gzip"
"""The default file extension for compressed gzip files (`.gzip`)."""

TRACK_EXTENSION: str = ".track"
"""The default file extension for UCSC track files."""

_BGZIP_EXTENSIONS: set[str] = {BGZ_EXTENSION, BGZIP_EXTENSION}
"""All supported block-compressed gzip file extensions."""

_GZIP_EXTENSIONS: set[str] = {GZ_EXTENSION, GZIP_EXTENSION}
"""All supported compressed gzip file extensions."""

_ALL_GZIP_COMPATIBLE_EXTENSIONS: set[str] = _BGZIP_EXTENSIONS.union(_GZIP_EXTENSIONS)
"""All supported compressed and block-compressed gzip file extensions."""

####################################################################################################

COMMENT_PREFIXES: set[str] = {"#", "browser", "track"}
"""The set of BED comment prefixes that this library supports."""

MISSING_FIELD: str = "."
"""The string used to indicate a missing field in a BED record."""

####################################################################################################

_BedKind = TypeVar("_BedKind", bound=bedspec._types.BedType)
"""A type variable for any kind of BED record type."""


class BedWriter(ContextManager, Generic[_BedKind]):
    """A writer of BED records.

    Args:
        handle: An open file-like object to write to.

    Attributes:
        bed_kind: The kind of BED type that this writer will write.

    """

    bed_kind: type[_BedKind] | None

    def __class_getitem__(cls, key: object) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return _GenericAlias(cls, key)  # type: ignore[no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedWriter[_BedKind]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        typelevel = signature.f_locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __enter__(self) -> "BedWriter[_BedKind]":
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
    def from_path(cls, path: Path | str) -> "BedWriter[_BedKind]":
        """Open a BED writer from a file path."""
        reader = cls(handle=Path(path).open("w"))  # type: ignore[operator]
        reader.bed_kind = None if len(cls.__args__) == 0 else cls.__args__[0]  # type: ignore[attr-defined]
        return cast("BedWriter[_BedKind]", reader)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()

    def write_comment(self, comment: str) -> None:
        """Write a comment to the BED output."""
        for line in comment.splitlines():
            prefix = "" if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES) else "# "
            self._handle.write(f"{prefix}{comment}\n")

    def write(self, bed: _BedKind) -> None:
        """Write a BED record to the BED output."""
        if self.bed_kind is not None:
            if type(bed) is not self.bed_kind:
                raise TypeError(
                    f"BedWriter can only continue to write features of the same type."
                    f" Will not write a {type(bed).__name__} after a {self.bed_kind.__name__}."
                )
        else:
            self.bed_kind = type(bed)

        mapping: dict[str, object] = {
            key: (value if value is not None else MISSING_FIELD)
            for key, value in as_dict(bed).items()
        }

        self._handle.write(f"{"\t".join(map(str, mapping.values()))}\n")

    def write_all(self, beds: Iterable[_BedKind]) -> None:
        """Write all the BED records to the BED output."""
        for bed in beds:
            self.write(bed)


class BedReader(ContextManager, Iterable[_BedKind], Generic[_BedKind]):
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

    bed_kind: type[_BedKind] | None

    def __class_getitem__(cls, key: object) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return _GenericAlias(cls, key)  # type: ignore[no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedReader[_BedKind]":
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

    def __enter__(self) -> "BedReader[_BedKind]":
        """Enter this context."""
        return self

    def __iter__(self) -> Iterator[_BedKind]:
        """Iterate through the BED records of this IO handle."""
        if self.bed_kind is None:
            raise NotImplementedError("Untyped reading is not yet supported!")
        for line in map(lambda line: line.strip(), self._handle):
            if line == "" or any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
                continue
            yield cast(_BedKind, self.bed_kind.decode(line))

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
    def from_path(cls, path: Path | str) -> "BedReader[_BedKind]":
        """Open a BED reader from a plaintext or gzip compressed file path."""
        reader = cls(handle=Path(path).open("r"))  # type: ignore[operator]
        reader.bed_kind = None if len(cls.__args__) == 0 else cls.__args__[0]  # type: ignore[attr-defined]
        return cast("BedReader[_BedKind]", reader)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()
