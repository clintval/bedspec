import io
import json
from csv import DictReader
from csv import DictWriter
from dataclasses import asdict as as_dict
from pathlib import Path
from types import NoneType
from types import TracebackType
from typing import Any
from typing import ContextManager
from typing import Generic
from typing import Iterable
from typing import Iterator
from typing import TypeAlias
from typing import TypeVar
from typing import cast
from typing import get_args
from typing import get_origin

from msgspec import convert
from msgspec import to_builtins
from typing_extensions import override

from bedspec._bedspec import BedColor
from bedspec._bedspec import BedLike
from bedspec._bedspec import BedStrand
from bedspec._bedspec import header
from bedspec._bedspec import types

BedType = TypeVar("BedType", bound=BedLike)
"""A type variable for any kind of BED record type."""

JsonType: TypeAlias = dict[str, "JsonType"] | list["JsonType"] | str | int | float | bool | None
"""A JSON-like data type."""

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


class BedWriter(ContextManager, Generic[BedType]):
    """
    A writer of BED records.

    Args:
        handle: An open file-like object to write to.

    Attributes:
        bed_type: The BED type that this writer will write.

    """

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED writer without knowing yet what BED types we will write."""
        self._bed_type: type[BedType] | None = None
        self._handle: io.TextIOWrapper = handle
        self._writer: DictWriter | None = None

    @property
    def bed_type(self) -> type[BedType] | None:
        return self._bed_type

    @bed_type.setter
    def bed_type(self, value: type[BedType]) -> None:
        self._bed_type: type[BedType] = value  # type: ignore[no-redef]
        self._header: list[str] = header(cast(BedLike, value))
        self._types: list[type | str | Any] = types(cast(BedLike, value))

    @override
    def __enter__(self) -> "BedWriter[BedType]":
        """Enter this context."""
        return self

    @override
    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return None

    def _maybe_setup_with(self, bed: BedType) -> None:
        """Perform post-initialization and record validation."""
        if self.bed_type is None:
            self.bed_type = type(bed)

        if self._writer is None:
            self._writer = DictWriter(self._handle, delimiter="\t", fieldnames=self._header)

        if not isinstance(bed, self.bed_type):
            raise TypeError(
                f"BedWriter can only continue to write features of the same type."
                f" Will not write a {type(bed).__name__} after a"
                f" {self.bed_type.__name__}."
            )

    def _bed_to_dict(self, bed: BedType) -> dict[str, object]:
        """Convert a BED record into a shallow dictionary."""
        shallow = {name: self._encode(getattr(bed, name)) for name in self._header}
        return cast(dict[str, object], to_builtins(shallow, order="deterministic"))

    @staticmethod
    def _encode(obj: Any) -> Any:
        """A callback for special encoding of custom types."""
        if obj is None:
            return "."
        if isinstance(obj, (list, set, tuple)):
            return ",".join(map(str, obj))
        if isinstance(obj, BedColor):
            return str(obj)
        return obj

    def write(self, bed: BedType) -> None:
        """Write a BED record to the BED output."""
        self._maybe_setup_with(bed)
        encoded = self._bed_to_dict(bed)
        self._writer.writerow(encoded)

    def write_comment(self, comment: str) -> None:
        """Write a comment to the BED output."""
        for line in comment.splitlines():
            prefix = "" if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES) else "# "
            self._handle.write(f"{prefix}{line}\n")

    @classmethod
    def from_path(cls, path: Path | str) -> "BedWriter":
        """Open a BED writer from a plain text file path."""
        writer: BedWriter = cls(handle=Path(path).open("w"))
        return writer

    def close(self) -> bool | None:
        """Close the underlying IO handle."""
        self._handle.close()
        return None


class BedReader(ContextManager, Iterable[BedType], Generic[BedType]):
    """
    A reader of BED records.

    This reader is capable of reading BED records but must be typed at runtime:

    ```python
    from bedspec import BedReader, Bed3

    with BedReader.from_path(path, Bed3) as reader:
        print(list(reader))
    ```

    Args:
        handle: An open file-like object to read from.

    Attributes:
        bed_type: The type of BED record that this reader will read.

    """

    def __init__(self, handle: io.TextIOWrapper, bed_type: type[BedType]) -> None:
        """Initialize a BED reader without knowing yet what BED types we will write."""
        self.bed_type: type[BedType] = bed_type
        self._handle: io.TextIOWrapper = handle
        self._header: list[str] = header(cast(BedLike, bed_type))
        self._types: list[type | str | Any] = types(cast(BedLike, bed_type))

        def without_comments() -> Iterable[str]:
            for line in self._handle:
                line = line.strip()
                if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
                    continue
                else:
                    yield line

        self._reader: DictReader = DictReader(
            without_comments(),
            delimiter="\t",
            fieldnames=self._header,
        )

    @override
    def __enter__(self) -> "BedReader[BedType]":
        """Enter this context."""
        return self

    @override
    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return None

    @override
    def __iter__(self) -> Iterator[BedType]:
        """Iterate through the BED records of this IO handle."""
        for bed in self._reader:
            yield convert(
                self._csv_dict_to_json(bed),
                self.bed_type,
                strict=False,
            )

        self.close()

    @staticmethod
    def _pre_decode(kind: type, obj: Any) -> Any:
        if obj == MISSING_FIELD and NoneType in get_args(kind):
            return None
        if kind is BedColor or BedColor in get_args(kind):
            if obj == "0":
                return None
            return json.dumps(as_dict(BedColor.from_string(cast(str, obj))))
        if kind is BedStrand or BedStrand in get_args(kind):
            return f'"{obj}"'
        return obj

    def _csv_dict_to_json(self, record: dict[str, str]) -> JsonType:
        """Convert a CSV dictionary record to JSON using the known field types."""
        key_values: list[str] = []
        for (name, value), field_type in zip(record.items(), self._types, strict=True):
            pre_decoded: str = self._pre_decode(cast(type, field_type), value)
            origin_type = get_origin(field_type)
            if pre_decoded is None:
                key_values.append(f'"{name}":null')
            elif origin_type is list:
                key_values.append(f'"{name}":[{pre_decoded.rstrip(",")}]')
            elif field_type is str or str in get_args(field_type):
                key_values.append(f'"{name}":"{pre_decoded}"')
            else:
                key_values.append(f'"{name}":{pre_decoded}')
        json_string: JsonType = json.loads(f"{{{','.join(key_values)}}}")
        return json_string

    @classmethod
    def from_path(cls, path: Path | str, bed_type: type[BedType]) -> "BedReader":
        """Open a BED reader from a plain text file path."""
        reader = cls(handle=Path(path).open("r"), bed_type=bed_type)
        return reader

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()
        return None
