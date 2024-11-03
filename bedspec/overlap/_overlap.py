from collections import defaultdict
from itertools import chain
from typing import Generic
from typing import Hashable
from typing import Iterable
from typing import Iterator
from typing import Protocol
from typing import TypeAlias
from typing import TypeVar
from typing import runtime_checkable

import cgranges as cr


@runtime_checkable
class Span(Hashable, Protocol):
    """A span with a start and an end. 0-based open-ended."""

    @property
    def start(self) -> int:
        """A 0-based start position."""
        raise NotImplementedError

    @property
    def end(self) -> int:
        """A 0-based open-ended position."""
        raise NotImplementedError


@runtime_checkable
class ReferenceSpan(Span, Protocol):
    """A feature on a reference sequence."""

    @property
    def refname(self) -> str:
        """A reference sequence name."""
        raise NotImplementedError


QueryReferenceSpanType = TypeVar("QueryReferenceSpanType", bound=ReferenceSpan)
"""Type variable for features being queried against the overlap detector."""

GenericReferenceSpanType = TypeVar("GenericReferenceSpanType", bound=ReferenceSpan)
"""Type variable for features stored within the overlap detector."""

Refname: TypeAlias = str
"""A type alias for a reference sequence name string."""


class OverlapDetector(Iterable[GenericReferenceSpanType], Generic[GenericReferenceSpanType]):
    """
    Detects and returns overlaps between a collection of reference features and query feature.

    The overlap detector may be built with any feature-like Python object that has the following
    properties:

      * `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based half-open end position

    This detector is most efficiently used when all features to be queried are added ahead of time.
    """

    def __init__(self, features: Iterable[GenericReferenceSpanType] | None = None) -> None:
        self._refname_to_features: dict[Refname, list[GenericReferenceSpanType]] = defaultdict(list)
        self._refname_to_tree: dict[Refname, cr.cgranges] = defaultdict(cr.cgranges)  # type: ignore[attr-defined,name-defined]
        self._refname_to_is_indexed: dict[Refname, bool] = defaultdict(lambda: False)
        if features is not None:
            self.add_all(features)

    def __iter__(self) -> Iterator[GenericReferenceSpanType]:
        """Iterate over the features in the overlap detector."""
        return chain(*self._refname_to_features.values())

    def add(self, feature: GenericReferenceSpanType) -> None:
        """Add a feature to this overlap detector."""
        if not isinstance(feature, Hashable):
            raise ValueError(f"Genomic feature is not hashable but should be: {feature}")

        refname: Refname = feature.refname
        feature_idx: int = len(self._refname_to_features[refname])

        self._refname_to_features[refname].append(feature)
        self._refname_to_tree[refname].add(refname, feature.start, feature.end, feature_idx)
        self._refname_to_is_indexed[refname] = False  # mark that this tree needs re-indexing

    def add_all(self, features: Iterable[GenericReferenceSpanType]) -> None:
        """Adds one or more features to this overlap detector."""
        for feature in features:
            self.add(feature)

    def overlapping(self, feature: QueryReferenceSpanType) -> Iterator[GenericReferenceSpanType]:
        """Yields all the overlapping features for a given query feature."""
        refname: Refname = feature.refname

        if refname in self._refname_to_tree and not self._refname_to_is_indexed[refname]:
            self._refname_to_tree[refname].index()  # index the tree if we find it is not indexed

        for *_, idx in self._refname_to_tree[refname].overlap(refname, feature.start, feature.end):
            yield self._refname_to_features[refname][idx]

    def overlaps(self, feature: QueryReferenceSpanType) -> bool:
        """Determine if a query feature overlaps any other features."""
        return next(self.overlapping(feature), None) is not None

    def enclosing(self, feature: QueryReferenceSpanType) -> Iterator[GenericReferenceSpanType]:
        """Yields all the overlapping features that completely enclose the given query feature."""
        for overlap in self.overlapping(feature):
            if feature.start >= overlap.start and feature.end <= overlap.end:
                yield overlap

    def enclosed_by(self, feature: QueryReferenceSpanType) -> Iterator[GenericReferenceSpanType]:
        """Yields all the overlapping features that are enclosed by the given query feature."""
        for overlap in self.overlapping(feature):
            if feature.start <= overlap.start and feature.end >= overlap.end:
                yield overlap
