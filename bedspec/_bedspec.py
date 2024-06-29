from abc import ABC
from dataclasses import dataclass
from typing import Iterator

from bedspec._types import BedStrand
from bedspec._types import BedType
from bedspec._types import GenomicSpan
from bedspec._types import Named
from bedspec._types import Stranded


class PointBed(BedType, ABC):
    """An abstract class for a BED record that describes a 0-based 1-length point."""

    contig: str
    start: int

    @property
    def length(self) -> int:
        """The length of this record."""
        return 1

    def territory(self) -> Iterator[GenomicSpan]:
        """Return the territory of a single point BED record which is 1-length."""
        yield Bed3(contig=self.contig, start=self.start, end=self.start + 1)


class SimpleBed(BedType, ABC, GenomicSpan):
    """An abstract class for a BED record that describes a contiguous linear interval."""

    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        """Validate this linear BED record."""
        if self.start >= self.end or self.start < 0:
            raise ValueError("Start must be greater than 0 and less than end!")

    @property
    def length(self) -> int:
        """The length of this record."""
        return self.end - self.start

    def territory(self) -> Iterator[GenomicSpan]:
        """Return the territory of a linear BED record which is just itself."""
        yield self


class PairBed(BedType, ABC):
    """An abstract base class for a BED record that describes a pair of linear linear intervals."""

    contig1: str
    start1: int
    end1: int
    contig2: str
    start2: int
    end2: int

    def __post_init__(self) -> None:
        """Validate this pair of BED records."""
        if self.start1 >= self.end1 or self.start1 < 0:
            raise ValueError("Start1 must be greater than 0 and less than end1!")
        if self.start2 >= self.end2 or self.start2 < 0:
            raise ValueError("Start2 must be greater than 0 and less than end2!")

    @property
    def bed1(self) -> SimpleBed:
        """The first of the two intervals."""
        return Bed3(contig=self.contig1, start=self.start1, end=self.end1)

    @property
    def bed2(self) -> SimpleBed:
        """The second of the two intervals."""
        return Bed3(contig=self.contig2, start=self.start2, end=self.end2)

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

    def __str__(self) -> str:
        """Return a comma-delimited string representation of this BED color."""
        return f"{self.r},{self.g},{self.b}"

    def __post_init__(self) -> None:
        """Validate that all color values are well-formatted."""
        if any(value > 255 or value < 0 for value in (self.r, self.g, self.b)):
            raise ValueError(f"RGB color values must be in the range [0, 255] but found: {self}")


@dataclass(eq=True, frozen=True)
class Bed2(PointBed):
    """A BED2 record that describes a single 0-based 1-length point."""

    contig: str
    start: int


@dataclass(eq=True, frozen=True)
class Bed3(SimpleBed):
    """A BED3 record that describes a contiguous linear interval."""

    contig: str
    start: int
    end: int


@dataclass(eq=True, frozen=True)
class Bed4(SimpleBed):
    """A BED4 record that describes a contiguous linear interval."""

    contig: str
    start: int
    end: int
    name: str | None


@dataclass(eq=True, frozen=True)
class Bed5(SimpleBed, Named):
    """A BED5 record that describes a contiguous linear interval."""

    contig: str
    start: int
    end: int
    name: str | None
    score: int | None


@dataclass(eq=True, frozen=True)
class Bed6(SimpleBed, Stranded, Named):
    """A BED6 record that describes a contiguous linear interval."""

    contig: str
    start: int
    end: int
    name: str | None
    score: int | None
    strand: BedStrand | None


@dataclass(eq=True, frozen=True)
class Bed12(SimpleBed, Stranded, Named):
    """A BED12 record that describes a contiguous linear interval."""

    contig: str
    start: int
    end: int
    name: str | None
    score: int | None
    strand: BedStrand | None
    thickStart: int | None
    thickEnd: int | None
    itemRgb: BedColor | None
    blockCount: int | None
    blockSizes: list[int]
    blockStarts: list[int]


@dataclass(eq=True, frozen=True)
class BedGraph(SimpleBed):
    """A bedGraph feature for continuous-valued data."""

    contig: str
    start: int
    end: int
    value: float


@dataclass(eq=True, frozen=True)
class BedPE(PairBed, Named):
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
