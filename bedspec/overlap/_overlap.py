from collections import defaultdict
from collections.abc import Iterable
from collections.abc import Iterator
from itertools import chain
from typing import Any
from typing import Generic
from typing import TypeAlias
from typing import TypeVar

from superintervals import IntervalMap
from typing_extensions import override

from bedspec._bedspec import BedStrand
from bedspec._bedspec import ReferenceSpan

ReferenceSpanType = TypeVar("ReferenceSpanType", bound=ReferenceSpan)
"""Type variable for features stored within the overlap detector."""

Refname: TypeAlias = str
"""A type alias for a reference sequence name string."""

IntervalTree: TypeAlias = IntervalMap
"""A type alias for the untyped interval map."""


def _strand(feature: Any) -> BedStrand | None:
    """Return the strand of a feature, or None if it has none."""
    strand: BedStrand | None = getattr(feature, "strand", None)
    return strand


def _closed(feature: ReferenceSpan) -> tuple[int, int]:
    """Return the closed interval of bases a feature covers, or flanks if it is zero-length."""
    if feature.start == feature.end:
        return max(feature.start - 1, 0), feature.start
    return feature.start, feature.end - 1


class OverlapDetector(Iterable[ReferenceSpanType], Generic[ReferenceSpanType]):
    """Detects and returns overlaps between a collection of reference features and query feature.

    The overlap detector may be built with any feature-like Python object that has the following
    properties:

      * `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based half-open end position

    A zero-length feature, such as an insertion, overlaps features holding either base beside it.

    Every query may be limited to features on the same strand as the query with `stranded=True`.
    A feature without a strand never matches a stranded query.

    This detector is most efficiently used when all features to be queried are added ahead of time.
    """

    def __init__(self, features: Iterable[ReferenceSpanType] | None = None) -> None:
        self._refname_to_features: dict[Refname, list[ReferenceSpanType]] = defaultdict(list)
        self._refname_to_tree: dict[Refname, IntervalTree] = defaultdict(IntervalTree)
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
            self._refname_to_tree[refname].add(*_closed(feature), feature_index)
            self._refname_to_is_indexed[refname] = False  # mark that this tree needs re-indexing

    def overlapping(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features for a given query feature."""
        refname: Refname = feature.refname

        if refname not in self._refname_to_tree:
            return

        if not self._refname_to_is_indexed[refname]:
            self._refname_to_tree[refname].build()
            self._refname_to_is_indexed[refname] = True

        features = self._refname_to_features[refname]
        indices: list[int] = self._refname_to_tree[refname].search_values(*_closed(feature))
        if not stranded:
            for index in indices:
                yield features[index]
            return

        strand = _strand(feature)
        if strand is None:
            return
        for index in indices:
            if _strand(features[index]) is strand:
                yield features[index]

    def overlaps(self, feature: ReferenceSpan, *, stranded: bool = False) -> bool:
        """Determine if a query feature overlaps any other features."""
        return next(self.overlapping(feature, stranded=stranded), None) is not None

    def enclosing(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features that completely enclose the given query feature."""
        for overlap in self.overlapping(feature, stranded=stranded):
            if feature.start >= overlap.start and feature.end <= overlap.end:
                yield overlap

    def enclosed_by(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features that are enclosed by the given query feature."""
        for overlap in self.overlapping(feature, stranded=stranded):
            if feature.start <= overlap.start and feature.end >= overlap.end:
                yield overlap
