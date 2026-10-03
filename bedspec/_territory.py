from array import array
from bisect import bisect_right
from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from itertools import chain
from typing import TypeAlias
from typing import final

from typing_extensions import override

from bedspec._bedspec import Bed3
from bedspec._bedspec import BedLike
from bedspec._bedspec import ReferenceSpan
from bedspec._bedspec import _check_span  # pyright: ignore[reportPrivateUsage]

Keep: TypeAlias = Callable[[bool, bool], bool]
"""Whether a base is kept, given whether it is in the first and in the second territory."""

_NONE: "array[int]" = array("q")
"""The bounds of a reference that holds no bases."""


def _spans(feature: BedLike | ReferenceSpan) -> Iterable[ReferenceSpan]:
    """Return every span of a BED record's territory, or the feature if it is not a BED record."""
    return feature.territory() if isinstance(feature, BedLike) else (feature,)


def _merged(spans: list[tuple[int, int]]) -> "array[int]":
    """Return the start and end of each of the fewest spans holding the bases of the given spans."""
    bounds = array("q")
    for start, end in sorted(spans):
        if bounds and start <= bounds[-1]:
            bounds[-1] = max(bounds[-1], end)
        else:
            bounds.extend((start, end))
    return bounds


def _swept(first: "array[int]", second: "array[int]", keep: Keep) -> "array[int]":
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

    A BED record adds every span of its `territory()`, so a `Bed2` adds its base and a `BedPE`
    both of its ends, and any other object with a `refname`, `start` and `end` adds itself.
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

    def _combined(self, other: object, keep: Keep) -> "Territory":
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
