import dataclasses
import inspect
import io
from abc import ABC
from abc import abstractmethod
from dataclasses import asdict as as_dict
from dataclasses import dataclass
from dataclasses import fields
from enum import StrEnum
from enum import unique
from functools import update_wrapper
from pathlib import Path
from types import FrameType
from types import TracebackType
from types import UnionType
from typing import Any
from typing import Callable
from typing import ClassVar
from typing import ContextManager
from typing import Generic
from typing import Iterable
from typing import Iterator
from typing import Protocol
from typing import Type
from typing import TypeVar
from typing import Union
from typing import _BaseGenericAlias  # type: ignore[attr-defined]
from typing import _GenericAlias  # type: ignore[attr-defined]
from typing import cast
from typing import get_args
from typing import get_origin
from typing import get_type_hints
from typing import runtime_checkable

COMMENT_PREFIXES: set[str] = {"#", "browser", "track"}
"""The set of BED comment prefixes supported by this implementation."""

MISSING_FIELD: str = "."
"""The string used to indicate a missing field in a BED record."""

BED_EXTENSION: str = ".bed"
"""The specification defined file extension for BED files."""

BEDPE_EXTENSION: str = ".bedpe"
"""The specification defined file extension for BedPE files."""


def is_union(annotation: Type) -> bool:
    """Test if we have a union type annotation or not."""
    return get_origin(annotation) in {Union, UnionType}


def is_optional(annotation: Type) -> bool:
    """Return if this type annotation is optional (a union type with None) or not."""
    return is_union(annotation) and type(None) in get_args(annotation)


def singular_non_optional_type(annotation: Type) -> Type:
    """Return the non-optional version of a singular type annotation."""
    if not is_optional(annotation):
        return annotation

    not_none: list[Type] = [arg for arg in get_args(annotation) if arg is not type(None)]
    if len(not_none) == 1:
        return not_none[0]
    else:
        raise TypeError(f"Complex non-optional types are not supported! Found: {not_none}")


class MethodType:
    def __init__(self, func: Callable, obj: object) -> None:
        self.__func__ = func
        self.__self__ = obj

    def __call__(self, *args: object, **kwargs: object) -> object:
        func = self.__func__
        obj = self.__self__
        return func(obj, *args, **kwargs)


class classmethod_generic:
    def __init__(self, f: Callable) -> None:
        self.f = f
        update_wrapper(self, f)

    def __get__(self, obj: object, cls: object | None = None) -> Callable:
        if cls is None:
            cls = type(obj)
        method = MethodType(self.f, cls)
        method._generic_classmethod = True  # type: ignore[attr-defined]
        return method


def __getattr__(self: object, name: str | None = None) -> object:
    if hasattr(obj := orig_getattr(self, name), "_generic_classmethod"):
        obj.__self__ = self
    return obj


orig_getattr = _BaseGenericAlias.__getattr__
_BaseGenericAlias.__getattr__ = __getattr__


@unique
class BedStrand(StrEnum):
    """Valid BED strands for forward, reverse, and unknown directions."""

    POSITIVE = "+"
    NEGATIVE = "-"

    def opposite(self) -> "BedStrand":
        """Return the opposite strand."""
        match self:
            case BedStrand.POSITIVE:
                return BedStrand.NEGATIVE
            case BedStrand.NEGATIVE:
                return BedStrand.POSITIVE


@dataclass
class BedColor:
    """The color of a BED record in red, green, and blue values."""

    r: int
    g: int
    b: int

    def __str__(self) -> str:
        """Return a string representation of this BED color."""
        return f"{self.r},{self.g},{self.b}"


@runtime_checkable
class DataclassProtocol(Protocol):
    """A protocol for objects that are dataclass instances."""

    __dataclass_fields__: ClassVar[dict[str, Any]]


@runtime_checkable
class Locatable(Protocol):
    """A protocol for 0-based half-open objects located on a reference sequence."""

    contig: str
    start: int
    end: int


@runtime_checkable
class Stranded(Protocol):
    """A protocol for stranded BED types."""

    strand: BedStrand | None


class BedType(ABC, DataclassProtocol):
    """An abstract base class for all types of BED records."""

    def __new__(cls, *args: object, **kwargs: object) -> "BedType":
        if not dataclasses.is_dataclass(cls):
            raise TypeError("You must mark custom BED records with @dataclass!")
        return cast("BedType", object.__new__(cls))

    @classmethod
    def decode(cls, line: str) -> "BedType":
        """Decode a line of text into a BED record."""
        row: list[str] = line.strip().split()
        coerced: dict[str, object] = {}

        try:
            zipped = list(zip(fields(cls), row, strict=True))
        except ValueError:
            raise ValueError(
                f"Expected {len(fields(cls))} fields but found {len(row)} in record:"
                f" '{' '.join(row)}'"
            ) from None

        hints: dict[str, Type] = get_type_hints(cls)

        for field, value in zipped:
            try:
                if is_optional(hints[field.name]) and value == MISSING_FIELD:
                    coerced[field.name] = None
                else:
                    coerced[field.name] = singular_non_optional_type(field.type)(value)
            except ValueError:
                raise TypeError(
                    f"Tried to build the BED field '{field.name}' (of type '{field.type.__name__}')"
                    f" from the value '{value}' but couldn't for record '{' '.join(row)}'"
                ) from None

        return cls(**coerced)

    @abstractmethod
    def territory(self) -> Iterator[Locatable]:
        """Return intervals that describe the territory of this BED record."""
        pass


class PointBed(BedType, ABC):
    """An abstract class for a BED record that describes a 0-based 1-length point."""

    contig: str
    start: int

    @property
    def length(self) -> int:
        """The length of this record."""
        return 1

    def territory(self) -> Iterator[Locatable]:
        """Return the territory of a single point BED record which is 1-length."""
        yield Bed3(contig=self.contig, start=self.start, end=self.start + 1)


class SimpleBed(BedType, ABC, Locatable):
    """An abstract class for a BED record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        """Validate this dataclass."""
        if self.start >= self.end or self.start < 0:
            raise ValueError("start must be greater than 0 and less than end!")

    @property
    def length(self) -> int:
        """The length of this record."""
        return self.end - self.start

    def territory(self) -> Iterator[Locatable]:
        """Return the territory of a simple BED record which is just itself."""
        yield self


class PairBed(BedType, ABC):
    """An abstract base class for a BED record that describes a pair of intervals."""

    contig1: str
    start1: int
    end1: int
    contig2: str
    start2: int
    end2: int

    def __post_init__(self) -> None:
        """Validate this dataclass."""
        if self.start1 >= self.end1 or self.start1 < 0:
            raise ValueError("start1 must be greater than 0 and less than end1!")
        if self.start2 >= self.end2 or self.start2 < 0:
            raise ValueError("start2 must be greater than 0 and less than end2!")

    @property
    def bed1(self) -> SimpleBed:
        """The first of the two intervals."""
        return Bed3(contig=self.contig1, start=self.start1, end=self.end1)

    @property
    def bed2(self) -> SimpleBed:
        """The second of the two intervals."""
        return Bed3(contig=self.contig2, start=self.start2, end=self.end2)

    def territory(self) -> Iterator[Locatable]:
        """Return the territory of this BED record which are two intervals."""
        yield self.bed1
        yield self.bed2


@dataclass
class Bed2(PointBed):
    """A BED2 record that describes a single 0-based 1-length point."""

    contig: str
    start: int


@dataclass
class Bed3(SimpleBed):
    """A BED3 record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int


@dataclass
class Bed4(SimpleBed):
    """A BED4 record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int
    name: str | None


@dataclass
class Bed5(SimpleBed):
    """A BED5 record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int
    name: str | None
    score: int | None


@dataclass
class Bed6(SimpleBed, Stranded):
    """A BED6 record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int
    name: str | None
    score: int | None
    strand: BedStrand | None


# @dataclass
# class Bed12(SimpleBed, Stranded):
#     """A BED12 record that describes a simple contiguous interval."""
#     contig: str
#     start: int
#     end: int
#     name: str
#     score: int
#     strand: BedStrand
#     thickStart: int
#     thickEnd: int
#     itemRgb: BedColor | None
#     blockCount: int
#     blockSizes: list[int]
#     blockStarts: list[int]

# TODO: Implement BED detail format? https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
# TODO: Implement bedGraph format? https://genome.ucsc.edu/goldenPath/help/bedgraph.html


@dataclass
class BedPE(PairBed):
    """A BED record that describes a pair of BED records as per the bedtools spec."""

    contig1: str
    start1: int
    end1: int
    contig2: str
    start2: int
    end2: int
    name: str | None
    score: int | None
    strand1: BedStrand | None
    strand2: BedStrand | None

    @property
    def bed1(self) -> Bed6:
        """The first of the two intervals as a BED6 record."""
        return Bed6(
            contig=self.contig1,
            start=self.start1,
            end=self.end1,
            name=self.name,
            score=self.score,
            strand=self.strand1,
        )

    @property
    def bed2(self) -> Bed6:
        """The second of the two intervals as a BED6 record."""
        return Bed6(
            contig=self.contig2,
            start=self.start2,
            end=self.end2,
            name=self.name,
            score=self.score,
            strand=self.strand2,
        )


BedKind = TypeVar("BedKind", bound=BedType)


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
        """Initialize a BED writer wihout knowing yet what BED types we will write."""
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

    @classmethod_generic
    def from_path(cls, path: Path | str) -> "BedWriter[BedKind]":
        """Open a BED reader from a file path."""
        reader = cls(handle=Path(path).open("w")) # type: ignore[operator]
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
        """Initialize a BED reader wihout knowing yet what BED types we will write."""
        self._handle = handle

    def __enter__(self) -> "BedReader[BedKind]":
        """Enter this context."""
        return self

    def __iter__(self) -> Iterator[BedKind]:
        """Iterate through the BED records of this IO handle."""
        # TODO: Implement __next__ and type this class as an iterator.
        if self.bed_kind is None:
            raise NotImplementedError("Untyped reading is not yet supported!")
        for line in self._handle:
            if line.strip() == "":
                continue
            if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
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

    @classmethod_generic
    def from_path(cls, path: Path | str) -> "BedReader[BedKind]":
        """Open a BED reader from a file path."""
        reader = cls(handle=Path(path).open()) # type: ignore[operator]
        reader.bed_kind = None if len(cls.__args__) == 0 else cls.__args__[0]  # type: ignore[attr-defined]
        return cast("BedReader[BedKind]", reader)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()
