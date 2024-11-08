from collections import defaultdict
from collections.abc import Iterable
from collections.abc import Iterator
from itertools import chain
from typing import Generic
from typing import TypeAlias
from typing import TypeVar

from typing_extensions import override

import cgranges as cr
from bedspec._bedspec import ReferenceSpan

ReferenceSpanType = TypeVar("ReferenceSpanType", bound=ReferenceSpan)
"""Type variable for features stored within the overlap detector."""

Refname: TypeAlias = str
"""A type alias for a reference sequence name string."""


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
        self._refname_to_tree: dict[Refname, cr.cgranges] = defaultdict(cr.cgranges)  # type: ignore[attr-defined,name-defined]  # pyright: ignore[reportUnknownArgumentType, reportUnknownMemberType]
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
            feature_idx: int = len(self._refname_to_features[refname])

            self._refname_to_features[refname].append(feature)
            self._refname_to_tree[refname].add(refname, feature.start, feature.end, feature_idx)  # pyright: ignore[reportUnknownMemberType]
            self._refname_to_is_indexed[refname] = False  # mark that this tree needs re-indexing

    def overlapping(self, feature: ReferenceSpan) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features for a given query feature."""
        refname: Refname = feature.refname

        if refname in self._refname_to_tree.keys() and not self._refname_to_is_indexed[refname]:  # pyright: ignore[reportUnknownMemberType]
            self._refname_to_tree[refname].index()  # pyright: ignore[reportUnknownMemberType]

        idx: int
        for *_, idx in self._refname_to_tree[refname].overlap(refname, feature.start, feature.end):  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
            yield self._refname_to_features[refname][idx]

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
