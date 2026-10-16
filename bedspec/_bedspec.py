import dataclasses
from abc import ABC
from abc import abstractmethod
from array import array
from bisect import bisect_right
from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from contextlib import AbstractContextManager
from dataclasses import Field
from dataclasses import dataclass
from dataclasses import field
from enum import Enum
from enum import unique
from itertools import chain
from types import TracebackType
from typing import Any
from typing import ClassVar
from typing import Protocol
from typing import TypeAlias
from typing import TypeVar
from typing import cast
from typing import final
from typing import runtime_checkable

from typeline import ExtraColumns
from typing_extensions import Self
from typing_extensions import TypeAliasType
from typing_extensions import is_protocol
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


def _thick_territory(
    record: "BedLike", refname: str, thick_start: int | None, thick_end: int | None
) -> "Territory":
    """Return the bases of a record's spans from its thick start to its thick end."""
    if thick_start is None or thick_end is None:
        return Territory()
    return record.territory() & Territory([Bed3(refname, start=thick_start, end=thick_end)])


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

    def thick_territory(self) -> "Territory":
        """Return the bases of this record in its thick part, none if the thick part is unset."""
        return _thick_territory(self, self.refname, self.thick_start, self.thick_end)


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

    def thick_territory(self) -> "Territory":
        """Return the bases of this record's blocks in its thick part, none if it is unset."""
        return _thick_territory(self, self.refname, self.thick_start, self.thick_end)

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


RecordType = TypeVar("RecordType")
"""A type variable for a class of records registered with `register_locatable`."""

SpansOf: TypeAlias = Callable[[Any], Iterable[ReferenceSpan]]
"""A function returning the 0-based half-open spans a record is found by."""

Locatable = TypeAliasType("Locatable", object)
"""A record bedspec can locate: a BED record, an object with a `refname`, `start` and `end`, a
record with a `spans()` method, or a record of a class registered with `register_locatable`.

It is a type alias of `object`, since a registration happens at runtime, and not a class.
"""

_REGISTRY: dict[type, SpansOf] = {}
"""The function registered for each class with `register_locatable`."""

_RESOLVED: dict[type, SpansOf | None] = {}
"""The registered function each class resolves to through its MRO, or None, until it changes."""

_UNREGISTERED: set[type] = set()
"""The classes known to resolve to no registered function, until the registry changes."""


def _registered_for(kind: type[Any]) -> SpansOf | None:
    """Return the function registered for a class or the nearest class in its MRO, if any."""
    try:
        return _RESOLVED[kind]
    except KeyError:
        spans_of: SpansOf | None = None
        if not issubclass(kind, BedLike):
            spans_of = next((_REGISTRY[base] for base in kind.__mro__ if base in _REGISTRY), None)
        _RESOLVED[kind] = spans_of
        if spans_of is None:
            _UNREGISTERED.add(kind)
        return spans_of


def _forget_resolutions() -> None:
    """Forget which function each class resolves to, after the registry changes."""
    _RESOLVED.clear()
    _UNREGISTERED.clear()


_INSTALLERS: dict[type, "_Registration"] = {}
"""The registration that installed the function registered for each class."""

_Installed: TypeAlias = tuple[SpansOf | None, "_Registration | None"]
"""A function registered for a class, or None, and the registration that installed it, or None."""


def _install(kind: type, installed: _Installed) -> None:
    """Register a function for a class with the registration that installed it, or neither."""
    spans_of, installer = installed
    if spans_of is None:
        _ = _REGISTRY.pop(kind, None)
    else:
        _REGISTRY[kind] = spans_of
    if installer is None:
        _ = _INSTALLERS.pop(kind, None)
    else:
        _INSTALLERS[kind] = installer
    _forget_resolutions()


class _Registration(AbstractContextManager[None]):
    """A registration that, used in a `with` statement, restores the one before it on exit."""

    def __init__(self, kind: type, spans_of: SpansOf, displaced: _Installed) -> None:
        self._kind: type = kind
        self._spans_of: SpansOf = spans_of
        self._displaced_by_call: _Installed | None = displaced
        self._displaced_by_blocks: list[_Installed] = []

    @override
    def __enter__(self) -> None:
        """Enter the block this registration holds for, registering it again if it was replaced."""
        _ = super().__enter__()
        kind = self._kind
        if _INSTALLERS.get(kind) is self and self._displaced_by_call is not None:
            self._displaced_by_blocks.append(self._displaced_by_call)
            self._displaced_by_call = None
        else:
            self._displaced_by_blocks.append((_REGISTRY.get(kind), _INSTALLERS.get(kind)))
            _install(kind, (self._spans_of, self))

    @override
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Restore what was registered for the class before the block, including nothing."""
        if self._displaced_by_blocks:
            _install(self._kind, self._displaced_by_blocks.pop())
        elif self._displaced_by_call is not None:
            _install(self._kind, self._displaced_by_call)
            self._displaced_by_call = None


def register_locatable(
    kind: type[RecordType], spans: Callable[[RecordType], Iterable[ReferenceSpan]]
) -> AbstractContextManager[None]:
    """Register how detectors and territories find the spans of a class of records.

    `spans` returns the 0-based half-open spans a record is found by, such as `Bed3` records, or
    `Bed6` records when strand matters; returning none means the record is found by nothing. It
    is used for the class and its subclasses, before anything else, replacing any earlier
    registration of the class. Used in a `with` statement, the registration holds only until the
    block exits, when what was registered before the block, or nothing, is restored; one replaced
    before its block is entered holds again for the block. A BED record is always found by its
    own `spans()`, even when a parent class of it is registered.

    Raises:
        TypeError: If `kind` is not a class, is `object`, a `Protocol`, or a BED class, or if
            `spans` is not callable.
    """
    if not isinstance(kind, type) or kind is object or is_protocol(kind):  # pyright: ignore[reportUnnecessaryIsInstance]
        raise TypeError(f"Only a class of records can be registered, not {kind!r}!")
    if issubclass(kind, BedLike):
        raise TypeError(
            "A BED record is found by its own spans: override spans() in a subclass of"
            + f" {kind.__name__} instead of registering it!"
        )
    if not callable(spans):
        raise TypeError(f"spans must be a function returning spans, not {spans!r}!")  # pyright: ignore[reportUnreachable]
    registration = _Registration(kind, spans, (_REGISTRY.get(kind), _INSTALLERS.get(kind)))
    _install(kind, (spans, registration))
    return registration


class _UnlocatableError(TypeError):
    """Raised for a record whose spans cannot be found."""


def _spans(feature: Locatable) -> Iterable[ReferenceSpan]:
    """Return the spans a feature is found by.

    In order: the function registered for its class, the spans of a BED record, the feature itself
    if it has a `refname`, `start` and `end`, or what its `spans()` method returns.

    Raises:
        TypeError: If none of these finds the spans of the feature.
    """
    if _REGISTRY:
        kind = type(feature)
        spans_of = None if kind in _UNREGISTERED else _registered_for(kind)
        if spans_of is not None:
            return spans_of(feature)
    if isinstance(feature, BedLike):
        return feature.spans()
    record: Any = feature
    try:
        _ = record.refname, record.start, record.end
        return (record,)
    except (AttributeError, KeyError):
        pass
    try:
        spans: Any = record.spans
    except (AttributeError, KeyError):
        spans = None
    if callable(spans):
        return cast("Iterable[ReferenceSpan]", spans())
    raise _UnlocatableError(
        f"Cannot find the spans of a {type(feature).__name__}: register a function returning them"
        + " with register_locatable(), or use one from bedspec.ext!"
    )


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
    ends, and a `Bed12` its blocks and not its introns, any other object with a `refname`,
    `start` and `end` adds itself, and a record registered with `register_locatable` adds the
    spans its function returns. A BED record's `territory()` is the territory of its spans.
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

    def __init__(self, features: Iterable[Locatable] = ()) -> None:
        """Build the territory of the given spans, BED records, and registered records.

        Raises:
            TypeError: If the spans of a feature cannot be found.
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

        A zero-length span is held when either base beside it is. A feature with no spans, like an
        unmapped read, or whose spans cannot be found, is not held.
        """
        try:
            spans = tuple(_spans(feature))
        except _UnlocatableError:
            return False
        if not spans:
            return False
        for span in spans:
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
