import dataclasses
from abc import ABC
from abc import abstractmethod
from array import array
from bisect import bisect_right
from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import Field
from dataclasses import dataclass
from dataclasses import field
from enum import Enum
from enum import unique
from itertools import chain
from typing import Any
from typing import ClassVar
from typing import Protocol
from typing import TypeAlias
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

    def __init_subclass__(cls) -> None:
        if "territory" in vars(cls):
            raise TypeError(
                "Override spans() in custom BED class definitions, not territory(), which is built"
                + " from spans()!"
            )
        return super().__init_subclass__()

    @abstractmethod
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield the pieces of this BED record on reference sequences, each with its strand."""

    @final
    def territory(self) -> "Territory":
        """Return the bases the spans of this BED record hold, without their strands."""
        return Territory(self.spans())


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
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield the one 1-length span of this point BED record."""
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
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield the one span of a linear BED record, which is itself."""
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
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield the two intervals of this BED record, each with its own strand."""
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

    @override
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield each block of this record on its strand, or the record if it has no blocks."""
        if self.block_starts is None or self.block_sizes is None:
            yield self
            return
        for block_start, block_size in zip(self.block_starts, self.block_sizes, strict=True):
            start = self.start + block_start
            yield Bed6(
                refname=self.refname,
                start=start,
                end=start + block_size,
                name=self.name,
                score=self.score,
                strand=self.strand,
            )


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


_Keep: TypeAlias = Callable[[bool, bool], bool]
"""Whether a base is kept, given whether it is in the first and in the second territory."""

_NONE: "array[int]" = array("q")
"""The bounds of a reference that holds no bases."""


def _spans(feature: BedLike | ReferenceSpan) -> Iterable[ReferenceSpan]:
    """Return the spans of a BED record, or the feature if it is not a BED record."""
    return feature.spans() if isinstance(feature, BedLike) else (feature,)


def _merged(spans: list[tuple[int, int]]) -> "array[int]":
    """Return the start and end of each of the fewest spans holding the bases of the given spans."""
    bounds = array("q")
    for start, end in sorted(spans):
        if bounds and start <= bounds[-1]:
            bounds[-1] = max(bounds[-1], end)
        else:
            bounds.extend((start, end))
    return bounds


def _swept(first: "array[int]", second: "array[int]", keep: _Keep) -> "array[int]":
    """Return the bounds of the bases kept, sweeping the bounds of two references in step."""
    bounds = array("q")
    i = j = 0
    inside = False
    while i < len(first) or j < len(second):
        if j == len(second) or (i < len(first) and first[i] <= second[j]):
            position = first[i]
        else:
            position = second[j]
        if i < len(first) and first[i] == position:
            i += 1
        if j < len(second) and second[j] == position:
            j += 1
        if keep(i % 2 == 1, j % 2 == 1) is not inside:
            inside = not inside
            bounds.append(position)
    return bounds


@final
class Territory(Iterable[Bed3]):
    """An immutable set of bases on reference sequences, built from spans and BED records.

    A BED record adds every one of its `spans()`, so a `Bed2` adds its base, a `BedPE` both of its
    ends, and a `Bed12` its blocks and not its introns, and any other object with a `refname`,
    `start` and `end` adds itself. A BED record's `territory()` is the territory of its spans.
    Spans that overlap or abut are joined, since they hold the same bases either way, a
    zero-length span holds no bases and adds none, and strands are not kept.

    Territories combine like sets of bases: `a | b` holds the bases of either, `a & b` the bases
    of both, and `a - b` the bases of `a` that are not in `b`. Territories are equal when they
    hold the same bases.

    Iterating yields the fewest `Bed3` spans holding exactly these bases, by start within each
    reference, with references in the order they were first added, those of the left operand
    first.

    Building a territory sorts the spans of each reference, and combining two territories sweeps
    their spans once.
    """

    __slots__ = ("_bounds",)

    def __init__(self, features: Iterable[BedLike | ReferenceSpan] = ()) -> None:
        """Build the territory of the given spans and BED records.

        Raises:
            ValueError: If a span has no reference name, starts before 0, or ends before it starts.
        """
        spans: dict[str, list[tuple[int, int]]] = {}
        for feature in features:
            for span in _spans(feature):
                _check_span(span.refname, span.start, span.end)
                if span.start < span.end:
                    spans.setdefault(span.refname, []).append((span.start, span.end))
        self._bounds: dict[str, array[int]] = {
            refname: _merged(starts_and_ends) for refname, starts_and_ends in spans.items()
        }

    @property
    def length(self) -> int:
        """The number of bases in this territory."""
        return sum(sum(bounds[1::2]) - sum(bounds[::2]) for bounds in self._bounds.values())

    def contains(self, refname: str, position: int) -> bool:
        """Return whether the base at a 0-based position on a reference is in this territory."""
        bounds = self._bounds.get(refname)
        return bounds is not None and bisect_right(bounds, position) % 2 == 1

    def __contains__(self, feature: object) -> bool:
        """Return whether this territory holds every base of every span of a feature.

        A zero-length span is held when either base beside it is.
        """
        if not isinstance(feature, BedLike | ReferenceSpan):
            return False
        for span in _spans(feature):
            bounds = self._bounds.get(span.refname, _NONE)
            index = bisect_right(bounds, span.start)
            if index % 2 == 1:
                if span.end > bounds[index]:
                    return False
            elif not (span.start == span.end and index > 0 and bounds[index - 1] == span.start):
                return False
        return True

    @override
    def __iter__(self) -> Iterator[Bed3]:
        """Yield the fewest spans holding exactly these bases, by reference and then by start."""
        for refname, bounds in self._bounds.items():
            for start, end in zip(bounds[::2], bounds[1::2], strict=True):
                yield Bed3(refname, start=start, end=end)

    def __bool__(self) -> bool:
        """Return whether this territory holds any bases."""
        return bool(self._bounds)

    def __or__(self, other: "Territory") -> "Territory":
        """Return the bases in this territory, the other, or both."""
        return self._combined(other, lambda first, second: first or second)

    def __and__(self, other: "Territory") -> "Territory":
        """Return the bases in both this territory and the other."""
        return self._combined(other, lambda first, second: first and second)

    def __sub__(self, other: "Territory") -> "Territory":
        """Return the bases in this territory that are not in the other."""
        return self._combined(other, lambda first, second: first and not second)

    def _combined(self, other: object, keep: _Keep) -> "Territory":
        """Return the bases kept, reference by reference, or NotImplemented for a non-territory."""
        if not isinstance(other, Territory):
            return NotImplemented  # type: ignore[no-any-return]
        combined = Territory()
        for refname in dict.fromkeys(chain(self._bounds, other._bounds)):
            first, second = self._bounds.get(refname, _NONE), other._bounds.get(refname, _NONE)
            bounds = _swept(first, second, keep)
            if bounds:
                combined._bounds[refname] = bounds
        return combined

    @override
    def __eq__(self, other: object) -> bool:
        """Return whether two territories hold the same bases."""
        if not isinstance(other, Territory):
            return NotImplemented
        return self._bounds == other._bounds

    @override
    def __hash__(self) -> int:
        """Return a hash of the bases in this territory."""
        return hash(
            frozenset((refname, bounds.tobytes()) for refname, bounds in self._bounds.items())
        )

    @override
    def __repr__(self) -> str:
        """Return a representation of this territory that builds it again."""
        return f"{type(self).__name__}({list(self)!r})"
