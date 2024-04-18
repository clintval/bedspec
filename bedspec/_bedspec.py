import dataclasses
import inspect
import io
import typing
from abc import ABC
from abc import abstractmethod
from csv import DictWriter
from dataclasses import Field
from dataclasses import asdict as as_dict
from dataclasses import dataclass
from dataclasses import fields
from enum import StrEnum
from types import FrameType
from types import TracebackType
from typing import Any
from typing import ClassVar
from typing import ContextManager
from typing import Generic
from typing import Iterator
from typing import Protocol
from typing import TypeVar
from typing import cast

COMMENT_PREFIXES: set[str] = {"#", "browser", "track"}
"""The set of BED comment prefixes supported by this implementation."""

MISSING_FIELD: str = "."
"""The string used to indicate a missing field in a BED record."""


class DataclassProtocol(Protocol):
    """A protocol for objects that are dataclass instances."""

    __dataclass_fields__: ClassVar[dict[str, Any]]


class Locatable(Protocol):
    """A protocol for 0-based half-open objects located on a reference sequence."""

    contig: str
    start: int
    end: int


class BedStrand(StrEnum):
    """Valid BED strands for forward, reverse, and unknown directions."""

    POSITIVE = "+"
    NEGATIVE = "-"
    UNKNOWN = MISSING_FIELD


class BedColor:
    """The color of a BED record in red, green, and blue values."""

    def __init__(self, r: int, g: int, b: int):
        """Build a new BED color from red, green, and blue values."""
        self.r = r
        self.g = g
        self.b = b

    def __str__(self) -> str:
        """Return a string representation of this BED color."""
        return f"{self.r},{self.g},{self.b}"


class BedType(ABC, DataclassProtocol):
    """An abstract base class for all types of BED records."""

    def __new__(cls, *args: Any, **kwargs: Any) -> "BedType":
        if not dataclasses.is_dataclass(cls):
            raise TypeError("You must mark custom BED records with @dataclass!")
        return cast("BedType", object.__new__(cls))

    @classmethod
    def fieldnames(cls) -> list[str]:
        """Return the names of this BED record's fields."""
        return [field.name for field in fields(cls)]

    @classmethod
    def decode(cls, line: str) -> "BedType":
        """Decode the target BED record from a line of text."""
        row: list[str] = line.strip().split()

        try:
            zipped: list[tuple[Field, str]] = list(zip(fields(cls), row, strict=True))
        except ValueError:
            raise ValueError(
                f"Expected {len(fields(cls))} fields but found {len(row)} in record:"
                f" '{' '.join(row)}'"
            ) from None

        coerced: dict[str, Any] = {}

        for field, value in zipped:
            try:
                coerced[field.name] = field.type(value)
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


class Stranded(Protocol):
    """A structural type for stranded BED types."""

    strand: BedStrand


class PointBed(BedType, ABC):
    """A BED record that describes a single 0-based point."""

    contig: str
    start: int

    def length(self) -> int:
        """The length of this record."""
        return 1

    def territory(self) -> Iterator[Locatable]:
        """Return the territory of a single point BED record which is 1-length."""
        yield Bed3(self.contig, self.start, self.start + 1)


class SimpleBed(BedType, ABC, Locatable):
    """A BED record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int

    def length(self) -> int:
        """The length of this record."""
        return self.end - self.start

    def territory(self) -> Iterator[Locatable]:
        """Return the territory of a simple BED record."""
        yield self


class PairBed(BedType, ABC):
    """A BED record that describes a pair of possibly-overlapping intervals."""

    contig1: str
    start1: int
    end1: int
    contig2: str
    start2: int
    end2: int

    @property
    def bed1(self) -> SimpleBed:
        """The first of the two intervals."""
        return Bed3(contig=self.contig1, start=self.start1, end=self.end1)

    @property
    def bed2(self) -> SimpleBed:
        """The second of the two intervals."""
        return Bed3(contig=self.contig2, start=self.start2, end=self.end2)

    def territory(self) -> Iterator[Locatable]:
        """Return the territory of this BED record which are both intervals."""
        yield self.bed1
        yield self.bed2


@dataclass
class Bed2(PointBed):
    """A BED2 record that describes a single 0-based point."""

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
    name: str


@dataclass
class Bed5(SimpleBed):
    """A BED5 record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int
    name: str
    score: int


@dataclass
class Bed6(SimpleBed, Stranded):
    """A BED6 record that describes a simple contiguous interval."""

    contig: str
    start: int
    end: int
    name: str
    score: int
    strand: BedStrand


# @dataclass
# class Bed12(SimpleBed):
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
    """A BED record that describes a pair of BED records per the bedtools spec."""

    contig1: str
    start1: int
    end1: int
    contig2: str
    start2: int
    end2: int
    name: str
    score: int
    strand1: BedStrand
    strand2: BedStrand

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


T = TypeVar("T", bound=BedType)


class BedWriter(Generic[T], ContextManager):
    bed_kind: type[T] | None

    def __class_getitem__(cls, key: Any) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return typing._GenericAlias(cls, key)  # type: ignore[attr-defined,no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedWriter[T]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        argvalues = inspect.getargvalues(signature)
        typelevel = argvalues.locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __enter__(self) -> "BedWriter[T]":
        """Enter this context."""
        return self

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED writer wihout knowing yet what BED types we will write."""
        self._handle = handle
        self._writer: DictWriter | None = None

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return super().__exit__(__exc_type, __exc_value, __traceback)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()

    def write(self, bed: T) -> None:
        """Write a BED record to the BED output."""
        if self.bed_kind is not None:
            if type(bed) is not self.bed_kind:
                raise TypeError(
                    f"BedWriter can only continue to write features of the same type."
                    f" Will not write a {type(bed).__name__} after a {self.bed_kind.__name__}."
                )
        else:
            self.bed_kind = type(bed)

        if self._writer is None:
            self._writer = DictWriter(
                self._handle, delimiter="\t", fieldnames=self.bed_kind.fieldnames()
            )

        self._writer.writerow(as_dict(bed))

    def write_comment(self, comment: str) -> None:
        """Write a comment to the BED output."""
        for line in comment.splitlines():
            if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
                self._handle.write(f"{comment}\n")
            else:
                self._handle.write(f"# {comment}\n")


class BedReader(Generic[T], ContextManager):
    bed_kind: type[T] | None

    def __class_getitem__(cls, key: Any) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return typing._GenericAlias(cls, key)  # type: ignore[attr-defined,no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedReader[T]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        argvalues = inspect.getargvalues(signature)
        typelevel = argvalues.locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED reader wihout knowing yet what BED types we will write."""
        self._handle = handle

    def __enter__(self) -> "BedReader[T]":
        """Enter this context."""
        return self

    def __iter__(self) -> Iterator[T]:
        """Iterate through the BED records of this IO handle."""
        for line in self._handle:
            if line.strip() == "":
                continue
            if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES):
                continue
            yield self._decode(line)

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return super().__exit__(__exc_type, __exc_value, __traceback)

    def _decode(self, line: str) -> T:
        if self.bed_kind is None:
            raise NotImplementedError("Untyped reading is not yet supported!")
        return cast(T, self.bed_kind.decode(line))

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()
