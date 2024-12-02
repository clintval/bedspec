from collections import defaultdict
from collections.abc import Iterable
from collections.abc import Iterator
from itertools import chain
from typing import Generic
from typing import TypeAlias
from typing import TypeVar

from superintervals import (  # type: ignore[import-untyped]  # pyright: ignore[reportMissingTypeStubs]
    IntervalSet,  # pyright: ignore[reportUnknownVariableType]
)
from typing_extensions import override

from bedspec._bedspec import ReferenceSpan

ReferenceSpanType = TypeVar("ReferenceSpanType", bound=ReferenceSpan)
"""Type variable for features stored within the overlap detector."""

Refname: TypeAlias = str
"""A type alias for a reference sequence name string."""

IntervalTree: TypeAlias = IntervalSet  # pyright: ignore[reportUnknownVariableType]
"""A type alias for the untyped interval set."""


class OverlapDetector(Iterable[ReferenceSpanType], Generic[ReferenceSpanType]):
    """Detects and returns overlaps between a collection of reference features and query feature.

    The overlap detector may be built with any feature-like Python object that has the following
    properties:

      * `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based half-open end position

    This detector is most efficiently used when all features to be queried are added ahead of time.
    """

    def __init__(self, features: Iterable[ReferenceSpanType] | None = None) -> None:
        self._refname_to_features: dict[Refname, list[ReferenceSpanType]] = defaultdict(list)
        self._refname_to_tree: dict[Refname, IntervalTree] = defaultdict(IntervalTree)  # pyright: ignore[reportUnknownArgumentType]
        self._refname_to_is_indexed: dict[Refname, bool] = defaultdict(lambda: False)
        if features is not None:
            self.add(*features)

    @override
    def __iter__(self) -> Iterator[ReferenceSpanType]:
        """Iterate over the features in the overlap detector."""
        return chain(*self._refname_to_features.values())

    def add(self, *features: ReferenceSpanType) -> None:
        """Add a feature to this overlap detector."""
        for feature in features:
            refname: Refname = feature.refname
            feature_index: int = len(self._refname_to_features[refname])

            self._refname_to_features[refname].append(feature)
            self._refname_to_tree[refname].add(feature.start, feature.end - 1, feature_index)  # pyright: ignore[reportUnknownMemberType]
            self._refname_to_is_indexed[refname] = False  # mark that this tree needs re-indexing

    def overlapping(self, feature: ReferenceSpan) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features for a given query feature."""
        refname: Refname = feature.refname

        if refname in self._refname_to_tree.keys() and not self._refname_to_is_indexed[refname]:  # pyright: ignore[reportUnknownMemberType]
            self._refname_to_tree[refname].index()  # pyright: ignore[reportUnknownMemberType]

        index: int
        for index in self._refname_to_tree[refname].find_overlaps(feature.start, feature.end - 1):  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
            yield self._refname_to_features[refname][index]

    def overlaps(self, feature: ReferenceSpan) -> bool:
        """Determine if a query feature overlaps any other features."""
        return next(self.overlapping(feature), None) is not None

    def enclosing(self, feature: ReferenceSpan) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features that completely enclose the given query feature."""
        for overlap in self.overlapping(feature):
            if feature.start >= overlap.start and feature.end <= overlap.end:
                yield overlap

    def enclosed_by(self, feature: ReferenceSpan) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features that are enclosed by the given query feature."""
        for overlap in self.overlapping(feature):
            if feature.start <= overlap.start and feature.end >= overlap.end:
                yield overlap
