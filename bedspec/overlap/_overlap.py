from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import dataclass
from dataclasses import field
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


@dataclass(slots=True)
class _RefIndex:
    """Per-reference interval index state."""

    tree: IntervalTree = field(default_factory=IntervalMap)
    is_built: bool = False


class TreeDetector(Iterable[ReferenceSpanType], Generic[ReferenceSpanType]):
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
    The `overlaps()` method is the cheapest way to test for any overlap, short-circuiting on the
    first hit.
    """

    def __init__(self, features: Iterable[ReferenceSpanType] | None = None) -> None:
        self._refname_to_index: dict[Refname, _RefIndex] = {}
        self._changes: int = 0
        if features is not None:
            self.add(*features)

    @override
    def __iter__(self) -> Iterator[ReferenceSpanType]:
        """Iterate over the features by reference, in the order first added, then by start.

        Queries may be made while iterating, but adding features raises a `RuntimeError`.
        """
        changes = self._changes
        trees = [self._built(index) for index in self._refname_to_index.values()]
        for tree in trees:
            for i in range(len(tree)):
                if self._changes != changes:
                    raise RuntimeError("TreeDetector changed during iteration")
                yield tree.data_at(i)

    def add(self, *features: ReferenceSpanType) -> None:
        """Add a feature to this overlap detector."""
        for feature in features:
            index = self._refname_to_index.get(feature.refname)
            if index is None:
                index = self._refname_to_index[feature.refname] = _RefIndex()
            index.tree.add(*_closed(feature), feature)
            index.is_built = False  # mark that this tree needs re-indexing
            self._changes += 1

    def _tree_for(self, refname: Refname) -> IntervalTree | None:
        """Return the built interval tree for a reference, or None if it has no features."""
        index = self._refname_to_index.get(refname)
        return None if index is None else self._built(index)

    @staticmethod
    def _built(index: _RefIndex) -> IntervalTree:
        """Return a reference's interval tree, building it first if features were added."""
        if not index.is_built:
            index.tree.build()
            index.is_built = True
        return index.tree

    def overlapping(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[ReferenceSpanType]:
        """Yields all the overlapping features for a given query feature."""
        tree = self._tree_for(feature.refname)
        if tree is None:
            return

        start, end = _closed(feature)
        if not stranded:
            yield from tree.search_values(start, end)
            return

        strand = _strand(feature)
        if strand is None:
            return
        for hit in tree.search_values(start, end):
            if _strand(hit) is strand:
                yield hit

    def overlaps(self, feature: ReferenceSpan, *, stranded: bool = False) -> bool:
        """Determine if a query feature overlaps any other features."""
        tree = self._tree_for(feature.refname)
        if tree is None:
            return False

        start, end = _closed(feature)
        if not stranded:
            # NB: superintervals' has_overlaps misses some nested overlaps; being fixed upstream.
            return next(tree.iter_idxs(start, end), None) is not None

        strand = _strand(feature)
        if strand is None:
            return False
        return any(_strand(hit) is strand for hit in tree.iter_values(start, end))

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
