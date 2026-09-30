import dataclasses
from abc import ABC
from abc import abstractmethod
from collections.abc import Iterator
from dataclasses import Field
from dataclasses import dataclass
from dataclasses import field
from enum import Enum
from enum import unique
from typing import Any
from typing import ClassVar
from typing import Protocol
from typing import TypeVar
from typing import final
from typing import runtime_checkable

from typeline import ExtraColumns
from typing_extensions import Self
from typing_extensions import override

COMMENT_PREFIXES: set[str] = {"#", "browser", "track"}
"""The set of BED comment prefixes that this library supports."""

MISSING_FIELD: str = "."
"""The string used to indicate a missing field in a BED record."""


def _check_span(refname: str, start: int, end: int, suffix: str = "") -> None:
    """Check that a span is named, starts at or after 0, and ends at or after its start."""
    if not refname:
        raise ValueError(f"refname{suffix} must not be empty!")
    if start < 0:
        raise ValueError(f"start{suffix} must be greater than or equal to 0!")
    if end < start:
        raise ValueError(f"end{suffix} must be greater than or equal to start{suffix}!")


def _check_name(name: str | None) -> None:
    """Check that a name, if given, is 1 to 255 characters long."""
    if name is not None and not 1 <= len(name) <= 255:
        raise ValueError("name must be 1 to 255 characters long!")


def _check_thick(start: int, end: int, thick_start: int | None, thick_end: int | None) -> None:
    """Check that a thick part is given whole or not at all, and sits within its feature."""
    if (thick_start is None) != (thick_end is None):
        raise ValueError("thick_start and thick_end must both be None or both be set!")
    if thick_start is not None and thick_end is not None:
        if not start <= thick_start <= thick_end <= end:
            raise ValueError(
                "thick_start and thick_end must satisfy start <= thick_start <= thick_end <= end!"
            )


def _check_score(score: int | None) -> None:
    """Check that a score, if given, is between 0 and 1000."""
    if score is not None and not 0 <= score <= 1000:
        raise ValueError("score must be between 0 and 1000!")


@runtime_checkable
class DataclassInstance(Protocol):
    """A protocol for objects that are dataclass instances."""

    __dataclass_fields__: ClassVar[dict[str, Field[Any]]]


@unique
class BedStrand(str, Enum):
    """BED strands for forward and reverse orientations."""

    Positive = "+"
    """The positive BED strand."""

    Negative = "-"
    """The negative BED strand."""

    def opposite(self) -> "BedStrand":
        """Return the opposite BED strand."""
        if self is BedStrand.Positive:
            return BedStrand.Negative
        else:
            return BedStrand.Positive

    @override
    def __str__(self) -> str:
        """Return this strand as a string."""
        return self.value


@runtime_checkable
class ReferenceSpan(Protocol):
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

    @abstractmethod
    def territory(self) -> Iterator[ReferenceSpan]:
        """Return intervals that describe the territory of this BED record."""


BedType = TypeVar("BedType", bound=BedLike)
"""A type variable for any kind of BED record type."""


@dataclass(frozen=True)
class PointBed(BedLike, ABC):
    """An abstract class for a BED record that describes a 0-based 1-length point."""

    refname: str
    start: int

    def __init_subclass__(cls) -> None:
        if not dataclasses.is_dataclass(cls):
            raise TypeError(
                "You must annotate custom BED class definitions with @dataclass(frozen=True)!"
            )
        return super().__init_subclass__()

    def __post_init__(self) -> None:
        """Validate this point BED record."""
        if not self.refname:
            raise ValueError("refname must not be empty!")
        if self.start < 0:
            raise ValueError("start must be greater than or equal to 0!")

    @final
    def __len__(self) -> int:
        """The length of this record."""
        return 1

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Return the territory of a single point BED record which is 1-length."""
        yield Bed3(refname=self.refname, start=self.start, end=self.start + 1)


@dataclass(frozen=True)
class SimpleBed(BedLike, ReferenceSpan, ABC):
    """An abstract class for a BED record that describes a contiguous linear interval."""

    refname: str
    start: int
    end: int

    def __init_subclass__(cls) -> None:
        if not dataclasses.is_dataclass(cls):
            raise TypeError(
                "You must annotate custom BED class definitions with @dataclass(frozen=True)!"
            )
        return super().__init_subclass__()

    def __post_init__(self) -> None:
        """Validate this linear BED record."""
        _check_span(self.refname, self.start, self.end)

    @final
    def __len__(self) -> int:
        """The length of this record."""
        return self.end - self.start

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Return the territory of a linear BED record which is just itself."""
        yield self


@dataclass(frozen=True)
class PairBed(BedLike, ABC):
    """An abstract base class for a BED record that describes a pair of linear linear intervals."""

    refname1: str
    start1: int
    end1: int
    refname2: str
    start2: int
    end2: int

    def __init_subclass__(cls) -> None:
        if not dataclasses.is_dataclass(cls):
            raise TypeError(
                "You must annotate custom BED class definitions with @dataclass(frozen=True)!"
            )
        return super().__init_subclass__()

    def __post_init__(self) -> None:
        """Validate this pair of BED records."""
        _check_span(self.refname1, self.start1, self.end1, suffix="1")
        _check_span(self.refname2, self.start2, self.end2, suffix="2")

    @property
    def bed1(self) -> SimpleBed:
        """The first of the two intervals."""
        return Bed3(refname=self.refname1, start=self.start1, end=self.end1)

    @property
    def bed2(self) -> SimpleBed:
        """The second of the two intervals."""
        return Bed3(refname=self.refname2, start=self.start2, end=self.end2)

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Return the territory of this BED record which are two intervals."""
        yield self.bed1
        yield self.bed2


@dataclass(slots=True, frozen=True)
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
    def from_string(cls, string: str) -> Self:
        """Build a BED color instance from a string."""
        try:
            r, g, b = map(int, string.split(","))
        except ValueError as error:
            raise ValueError(f"Invalid string '{string}'. Expected 'int,int,int'!") from error
        return cls(r, g, b)

    @override
    def __str__(self) -> str:
        """Return a comma-delimited string representation of this BED color."""
        return f"{self.r},{self.g},{self.b}"


@dataclass(slots=True, frozen=True)
class Bed2(PointBed):
    """A BED2 record that describes a single 0-based 1-length point."""

    refname: str
    start: int


@dataclass(slots=True, frozen=True)
class Bed3(SimpleBed):
    """A BED3 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)


@dataclass(slots=True, frozen=True)
class Bed4(SimpleBed):
    """A BED4 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)

    @override
    def __post_init__(self) -> None:
        """Validate this BED4 record."""
        super(Bed4, self).__post_init__()
        _check_name(self.name)


@dataclass(slots=True, frozen=True)
class Bed5(SimpleBed, Named):
    """A BED5 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)

    @override
    def __post_init__(self) -> None:
        """Validate this BED5 record."""
        super(Bed5, self).__post_init__()
        _check_name(self.name)
        _check_score(self.score)


@dataclass(slots=True, frozen=True)
class Bed6(SimpleBed, Named, Stranded):
    """A BED6 record that describes a contiguous linear interval."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)
    strand: BedStrand | None = field(kw_only=True)

    @override
    def __post_init__(self) -> None:
        """Validate this BED6 record."""
        super(Bed6, self).__post_init__()
        _check_name(self.name)
        _check_score(self.score)


@dataclass(slots=True, frozen=True)
class Bed9(SimpleBed, Named, Stranded):
    """A BED9 record that describes a contiguous linear interval with a thick part and a color."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    name: str | None = field(kw_only=True)
    score: int | None = field(kw_only=True)
    strand: BedStrand | None = field(kw_only=True)
    thick_start: int | None = field(kw_only=True)
    thick_end: int | None = field(kw_only=True)
    item_rgb: BedColor | None = field(kw_only=True)

    @override
    def __post_init__(self) -> None:
        """Validate this BED9 record."""
        super(Bed9, self).__post_init__()
        _check_name(self.name)
        _check_score(self.score)
        _check_thick(self.start, self.end, self.thick_start, self.thick_end)


@dataclass(slots=True, frozen=True)
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
    block_sizes: tuple[int, ...] | None = field(kw_only=True)
    block_starts: tuple[int, ...] | None = field(kw_only=True)

    def __post_init__(self) -> None:
        """Validate this BED12 record."""
        super(Bed12, self).__post_init__()
        _check_name(self.name)
        _check_score(self.score)
        _check_thick(self.start, self.end, self.thick_start, self.thick_end)
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
            if any(
                start < previous + size
                for previous, size, start in zip(
                    self.block_starts, self.block_sizes, self.block_starts[1:], strict=False
                )
            ):
                raise ValueError("Blocks must be in ascending order and must not overlap!")
            if (self.start + self.block_starts[-1] + self.block_sizes[-1]) != self.end:
                raise ValueError("The last defined block's end must be equal to the BED end!")


@dataclass(slots=True, frozen=True)
class Bed3N(Bed3):
    """A BED3+N record: a BED3 record followed by any number of extra columns, kept as text."""

    extra: ExtraColumns = field(default=(), kw_only=True)


@dataclass(slots=True, frozen=True)
class Bed4N(Bed4):
    """A BED4+N record: a BED4 record followed by any number of extra columns, kept as text."""

    extra: ExtraColumns = field(default=(), kw_only=True)


@dataclass(slots=True, frozen=True)
class Bed5N(Bed5):
    """A BED5+N record: a BED5 record followed by any number of extra columns, kept as text."""

    extra: ExtraColumns = field(default=(), kw_only=True)


@dataclass(slots=True, frozen=True)
class Bed6N(Bed6):
    """A BED6+N record: a BED6 record followed by any number of extra columns, kept as text."""

    extra: ExtraColumns = field(default=(), kw_only=True)


@dataclass(slots=True, frozen=True)
class Bed9N(Bed9):
    """A BED9+N record: a BED9 record followed by any number of extra columns, kept as text."""

    extra: ExtraColumns = field(default=(), kw_only=True)


@dataclass(slots=True, frozen=True)
class Bed12N(Bed12):
    """A BED12+N record: a BED12 record followed by any number of extra columns, kept as text."""

    extra: ExtraColumns = field(default=(), kw_only=True)


@dataclass(slots=True, frozen=True)
class NarrowPeak(Bed6):
    """An ENCODE narrowPeak (BED6+4) record of a peak with a summit.

    ENCODE writes -1 for a p-value, q-value, or summit that is not given.
    """

    signal_value: float = field(kw_only=True)
    """The overall enrichment of the peak."""

    p_value: float = field(kw_only=True)
    """The -log10 p-value of the peak, or -1."""

    q_value: float = field(kw_only=True)
    """The -log10 q-value of the peak, or -1."""

    peak: int = field(kw_only=True)
    """The summit, as a 0-based offset from start, or -1."""

    @override
    def __post_init__(self) -> None:
        """Validate this narrowPeak record."""
        super(NarrowPeak, self).__post_init__()
        if self.peak != -1 and not 0 <= self.peak < self.end - self.start:
            raise ValueError("peak must be -1 or an offset within the feature!")


@dataclass(slots=True, frozen=True)
class BroadPeak(Bed6):
    """An ENCODE broadPeak (BED6+3) record of a broad region of enrichment.

    ENCODE writes -1 for a p-value or q-value that is not given.
    """

    signal_value: float = field(kw_only=True)
    """The overall enrichment of the region."""

    p_value: float = field(kw_only=True)
    """The -log10 p-value of the region, or -1."""

    q_value: float = field(kw_only=True)
    """The -log10 q-value of the region, or -1."""


@dataclass(slots=True, frozen=True)
class GappedPeak(Bed12):
    """An ENCODE gappedPeak (BED12+3) record of a peak made of blocks.

    ENCODE writes -1 for a p-value or q-value that is not given.
    """

    signal_value: float = field(kw_only=True)
    """The overall enrichment of the peak."""

    p_value: float = field(kw_only=True)
    """The -log10 p-value of the peak, or -1."""

    q_value: float = field(kw_only=True)
    """The -log10 q-value of the peak, or -1."""


@dataclass(slots=True, frozen=True)
class BedGraph(SimpleBed):
    """A bedGraph feature for continuous-valued data."""

    refname: str
    start: int = field(kw_only=True)
    end: int = field(kw_only=True)
    value: float = field(kw_only=True)


@dataclass(slots=True, frozen=True)
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
    def _bed_score(self) -> int | None:
        """This pair's score if BED can hold it, since a BEDPE score is not limited to 0-1000."""
        return self.score if self.score is None or 0 <= self.score <= 1000 else None

    @property
    @override
    def bed1(self) -> Bed6:
        """The first of the two intervals as a BED6 record, without a score BED can't hold."""
        return Bed6(
            refname=self.refname1,
            start=self.start1,
            end=self.end1,
            name=self.name,
            score=self._bed_score,
            strand=self.strand1,
        )

    @property
    @override
    def bed2(self) -> Bed6:
        """The second of the two intervals as a BED6 record, without a score BED can't hold."""
        return Bed6(
            refname=self.refname2,
            start=self.start2,
            end=self.end2,
            name=self.name,
            score=self._bed_score,
            strand=self.strand2,
        )

    @classmethod
    def from_bed6(
        cls, bed1: Bed6, bed2: Bed6, name: str | None = None, score: int | None = None
    ) -> Self:
        return cls(
            refname1=bed1.refname,
            start1=bed1.start,
            end1=bed1.end,
            refname2=bed2.refname,
            start2=bed2.start,
            end2=bed2.end,
            name=name,
            score=score,
            strand1=bed1.strand,
            strand2=bed2.strand,
        )
