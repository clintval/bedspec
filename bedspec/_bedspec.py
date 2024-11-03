import dataclasses
from abc import ABC
from abc import abstractmethod
from collections.abc import Iterator
from dataclasses import Field
from dataclasses import dataclass
from dataclasses import field
from dataclasses import fields
from enum import StrEnum
from enum import unique
from typing import Any
from typing import ClassVar
from typing import Protocol
from typing import final
from typing import runtime_checkable

from typing_extensions import override


@runtime_checkable
class DataclassInstance(Protocol):
    """A protocol for objects that are dataclass instances."""

    __dataclass_fields__: ClassVar[dict[str, Field[Any]]]


@unique
class BedStrand(StrEnum):
    """BED strands for forward and reverse orientations."""

    Positive = "+"
    """The positive BED strand."""

    Negative = "-"
    """The negative BED strand."""

    def opposite(self) -> "BedStrand":
        """Return the opposite BED strand."""
        return BedStrand.Negative if self is BedStrand.Positive else BedStrand.Positive


@runtime_checkable
class GenomicSpan(Protocol):
    """A structural protocol for 0-based half-open objects located on a reference sequence."""

    refname: str
    start: int
    end: int


@runtime_checkable
class Named(Protocol):
    """A structural protocol for a named BED type."""

    name: str | None


@runtime_checkable
class Stranded(Protocol):
    """A structural protocol for stranded BED types."""

    strand: BedStrand | None


class BedLike(ABC, DataclassInstance):
    """An abstract base class for all types of BED records."""

    def __new__(cls, *_: object, **__: object) -> "BedLike":
        if not dataclasses.is_dataclass(cls):
            raise TypeError("You must annotate custom BED class definitions with @dataclass!")
        instance: BedLike = object.__new__(cls)
        return instance

    @abstractmethod
    def territory(self) -> Iterator[GenomicSpan]:
        """Return intervals that describe the territory of this BED record."""


def header(bed: BedLike | type[BedLike]) -> list[str]:
    """Return the list of field names for this BED record."""
    return [field.name for field in fields(bed)]


def types(bed: BedLike | type[BedLike]) -> list[type | str | Any]:
    """Return the list of field types for this BED record."""
    return [field.type for field in fields(bed)]


class PointBed(BedLike, ABC):
    """An abstract class for a BED record that describes a 0-based 1-length point."""

    refname: str
    start: int

    @final
    @property
    def length(self) -> int:
        """The length of this record."""
        return 1

    @override
    def territory(self) -> Iterator[GenomicSpan]:
        """Return the territory of a single point BED record which is 1-length."""
        yield Bed3(refname=self.refname, start=self.start, end=self.start + 1)


class SimpleBed(BedLike, GenomicSpan, ABC):
    """An abstract class for a BED record that describes a contiguous linear interval."""

    refname: str
    start: int
    end: int

    def __post_init__(self) -> None:
        """Validate this linear BED record."""
        if self.start >= self.end or self.start < 0:
            raise ValueError("start must be greater than 0 and less than end!")

    @final
    @property
    def length(self) -> int:
        """The length of this record."""
        return self.end - self.start

    @override
    def territory(self) -> Iterator[GenomicSpan]:
        """Return the territory of a linear BED record which is just itself."""
        yield self


class PairBed(BedLike, ABC):
    """An abstract base class for a BED record that describes a pair of linear linear intervals."""

    refname1: str
    start1: int
    end1: int
    refname2: str
    start2: int
    end2: int

    def __post_init__(self) -> None:
        """Validate this pair of BED records."""
        if self.start1 >= self.end1 or self.start1 < 0:
            raise ValueError("start1 must be greater than 0 and less than end1!")
        if self.start2 >= self.end2 or self.start2 < 0:
            raise ValueError("start2 must be greater than 0 and less than end2!")

    @property
    def bed1(self) -> SimpleBed:
        """The first of the two intervals."""
        return Bed3(refname=self.refname1, start=self.start1, end=self.end1)

    @property
    def bed2(self) -> SimpleBed:
        """The second of the two intervals."""
        return Bed3(refname=self.refname2, start=self.start2, end=self.end2)

    def territory(self) -> Iterator[GenomicSpan]:
        """Return the territory of this BED record which are two intervals."""
        yield self.bed1
        yield self.bed2


@dataclass(eq=True, frozen=True)
class BedColor:
    """The color of a BED record in red, green, and blue color values."""

    r: int
    g: int
    b: int

    def __post_init__(self) -> None:
        """Validate that all color values are well-formatted."""
        if any(value > 255 or value < 0 for value in (self.r, self.g, self.b)):
            raise ValueError(f"RGB color values must be in the range [0, 255] but found: {self}")

    @classmethod
    def from_string(cls, string: str) -> "BedColor":
        """Build a BED color instance from a string."""
        try:
            r, g, b = map(int, string.split(","))
        except ValueError as error:
            raise ValueError(f"Invalid string '{string}'. Expected 'int,int,int'!") from error
        return cls(r, g, b)

    def __str__(self) -> str:
        """Return a comma-delimited string representation of this BED color."""
        return f"{self.r},{self.g},{self.b}"


@dataclass(eq=True, frozen=True)
class Bed2(PointBed):
    """A BED2 record that describes a single 0-based 1-length point."""

    refname: str
    start: int


@dataclass(eq=True, frozen=True)
class Bed3(SimpleBed):
    """A BED3 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)


@dataclass(eq=True, frozen=True)
class Bed4(SimpleBed):
    """A BED4 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)


@dataclass(eq=True, frozen=True)
class Bed5(SimpleBed, Named):
    """A BED5 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)


@dataclass(eq=True, frozen=True)
class Bed6(SimpleBed, Named, Stranded):
    """A BED6 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)
    strand: BedStrand | None = field(kw_only=True)


@dataclass(eq=True, frozen=True)
class Bed12(SimpleBed, Named, Stranded):
    """A BED12 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)
    strand: BedStrand | None = field(kw_only=True)
    thick_start: int | None = field(kw_only=True)
    thick_end: int | None = field(kw_only=True)
    item_rgb: BedColor | None = field(kw_only=True)
    block_count: int | None = field(kw_only=True)
    block_sizes: list[int] = field(kw_only=True)
    block_starts: list[int] = field(kw_only=True)

    def __post_init__(self) -> None:
        """Validate this BED12 record."""
        super().__post_init__()
        if (self.thick_start is None) != (self.thick_end is None):
            raise ValueError("thick_start and thick_end must both be None or both be set!")
        if self.block_count is None:
            if self.block_sizes is not None or self.block_starts is not None:
                raise ValueError("block_count, block_sizes, block_starts must all be set or unset!")
        else:
            if self.block_sizes is None or self.block_starts is None:
                raise ValueError("block_count, block_sizes, block_starts must all be set or unset!")
            if self.block_count <= 0:
                raise ValueError("When set, block_count must be greater than or equal to 1!")
            if self.block_count != len(self.block_sizes) or self.block_count != len(
                self.block_starts
            ):
                raise ValueError("Length of block_sizes and block_starts must equal block_count!")
            if self.block_starts[0] != 0:
                raise ValueError("block_starts must start with 0!")
            if any(size <= 0 for size in self.block_sizes):
                raise ValueError("All sizes in block_size must be greater than or equal to one!")
            if (self.start + self.block_starts[-1] + self.block_sizes[-1]) != self.end:
                raise ValueError("The last defined block's end must be equal to the BED end!")


@dataclass(eq=True, frozen=True)
class BedGraph(SimpleBed):
    """A bedGraph feature for continuous-valued data."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    value: float = field(kw_only=True)


@dataclass(eq=True, frozen=True)
class BedPE(PairBed, Named):
    """A BED record that describes a pair of BED records as per the bedtools spec."""

    refname1: str = field(kw_only=True)
    start1: int = field(kw_only=True)
    end1: int = field(kw_only=True)
    refname2: str = field(kw_only=True)
    start2: int = field(kw_only=True)
    end2: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)
    strand1: BedStrand | None = field(kw_only=True)
    strand2: BedStrand | None = field(kw_only=True)

    @property
    def bed1(self) -> Bed6:
        """The first of the two intervals as a BED6 record."""
        return Bed6(
            refname=self.refname1,
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
            refname=self.refname2,
            start=self.start2,
            end=self.end2,
            name=self.name,
            score=self.score,
            strand=self.strand2,
        )
