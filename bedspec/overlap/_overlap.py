from bisect import bisect_left
from bisect import bisect_right
from collections.abc import Callable
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

from bedspec._bedspec import BedLike
from bedspec._bedspec import BedStrand
from bedspec._bedspec import ReferenceSpan

FeatureType = TypeVar("FeatureType", bound=BedLike | ReferenceSpan)
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


def _territory(feature: Any) -> Iterable[ReferenceSpan]:
    """Return the spans a feature covers: its territory if it has one, or else the feature."""
    territory: Callable[[], Iterable[ReferenceSpan]] | None = getattr(feature, "territory", None)
    return (feature,) if territory is None else territory()


@dataclass(slots=True)
class _RefIndex:
    """Per-reference interval index state."""

    tree: IntervalTree = field(default_factory=IntervalMap)
    is_built: bool = False


class TreeDetector(Iterable[FeatureType], Generic[FeatureType]):
    """Detects and returns overlaps between a collection of reference features and query feature.

    The overlap detector may be built with any BED record, indexed by every span of its
    `territory()`, or with any feature-like Python object that has the following properties:

      * `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based half-open end position

    A feature matches a query when any of its spans on the query's reference does, so a BEDPE
    record is found by either end, but it is enclosed by a query only when all of its spans are.
    A query yields each matching feature once, in the order the index finds it, which is
    repeatable but is not the order the features were added.

    A zero-length feature, such as an insertion, overlaps features holding either base beside it.

    Every query may be limited to features on the same strand as the query with `stranded=True`.
    The strand compared is that of each matching span, or of its feature if the span has none.
    A feature without a strand never matches a stranded query.

    This detector is most efficiently used when all features to be queried are added ahead of time.
    The `overlaps()` method is the cheapest way to test for any overlap, short-circuiting on the
    first hit.
    """

    def __init__(self, features: Iterable[FeatureType] | None = None) -> None:
        self._features: list[FeatureType] = []
        self._spans: list[ReferenceSpan] = []
        self._span_owners: list[int] = []
        self._one_span_each: bool = True
        self._refname_to_index: dict[Refname, _RefIndex] = {}
        if features is not None:
            self.add(*features)

    @override
    def __iter__(self) -> Iterator[FeatureType]:
        """Iterate over the features in the order they were added.

        Queries may be made while iterating, but adding features raises a `RuntimeError`.
        """
        count = len(self._features)
        for feature in self._features:
            if len(self._features) != count:
                raise RuntimeError("TreeDetector changed during iteration")
            yield feature

    def add(self, *features: FeatureType) -> None:
        """Add features to this overlap detector, indexing every span of each."""
        added = self._features
        spans = self._spans
        span_owners = self._span_owners
        refname_to_index = self._refname_to_index
        for feature in features:
            number = len(added)
            added.append(feature)
            for span in _territory(feature):
                index = refname_to_index.get(span.refname)
                if index is None:
                    index = refname_to_index[span.refname] = _RefIndex()
                index.tree.add(*_closed(span), len(spans))
                index.is_built = False  # mark that this tree needs re-indexing
                spans.append(span)
                span_owners.append(number)
            if len(spans) != len(added):
                self._one_span_each = False

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

    def _strand_of(self, span: int) -> BedStrand | None:
        """Return the strand of a span, or of its feature if the span has none."""
        return _strand(self._spans[span]) or _strand(self._features[self._span_owners[span]])

    def _spans_of(self, number: int) -> range:
        """Return the spans of a feature, which are numbered consecutively as they were added."""
        span_owners = self._span_owners
        return range(bisect_left(span_owners, number), bisect_right(span_owners, number))

    def _is_inside(self, span: int, feature: ReferenceSpan, strand: BedStrand | None) -> bool:
        """Return whether a span is inside a query feature, and on a strand if one is given."""
        inner = self._spans[span]
        return (
            inner.refname == feature.refname
            and feature.start <= inner.start
            and feature.end >= inner.end
            and (strand is None or self._strand_of(span) is strand)
        )

    def _hits(self, feature: ReferenceSpan, stranded: bool) -> list[int]:
        """Return the spans that overlap a query feature, on its strand if stranded."""
        tree = self._tree_for(feature.refname)
        if tree is None:
            return []

        hits: list[int] = tree.search_values(*_closed(feature))
        if not stranded:
            return hits

        strand = _strand(feature)
        if strand is None:
            return []
        return [hit for hit in hits if self._strand_of(hit) is strand]

    def _features_of(self, spans: list[int]) -> list[FeatureType]:
        """Return the features that own the given spans, once each, in the order of the spans."""
        features = self._features
        if self._one_span_each:  # each span's number is then its feature's number
            return [features[span] for span in spans]
        span_owners = self._span_owners
        return [features[number] for number in dict.fromkeys(span_owners[span] for span in spans)]

    def overlapping(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[FeatureType]:
        """Yields all the features with a span that overlaps the given query feature."""
        yield from self._features_of(self._hits(feature, stranded))

    def overlaps(self, feature: ReferenceSpan, *, stranded: bool = False) -> bool:
        """Determine if a query feature overlaps any other features."""
        tree = self._tree_for(feature.refname)
        if tree is None:
            return False

        start, end = _closed(feature)
        if not stranded:
            return tree.has_overlaps(start, end)

        strand = _strand(feature)
        if strand is None:
            return False
        return any(self._strand_of(span) is strand for span in tree.iter_values(start, end))

    def enclosing(self, feature: ReferenceSpan, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with a span that completely encloses the given query feature."""
        spans = self._spans
        yield from self._features_of([
            hit
            for hit in self._hits(feature, stranded)
            if feature.start >= spans[hit].start and feature.end <= spans[hit].end
        ])

    def enclosed_by(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[FeatureType]:
        """Yields all the features with every span enclosed by the given query feature."""
        spans = self._spans
        hits = self._hits(feature, stranded)
        if self._one_span_each:
            yield from self._features_of([
                hit
                for hit in hits
                if feature.start <= spans[hit].start and feature.end >= spans[hit].end
            ])
            return

        features = self._features
        strand = _strand(feature) if stranded else None
        for number in dict.fromkeys(self._span_owners[hit] for hit in hits):
            if all(self._is_inside(span, feature, strand) for span in self._spans_of(number)):
                yield features[number]
