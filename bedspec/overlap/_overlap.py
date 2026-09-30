from array import array
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
from bedspec._bedspec import SimpleBed

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


def _is_own_span(feature: Any) -> bool:
    """Return whether a feature is its own only span, as a record that is not BED-like is."""
    territory = getattr(feature.__class__, "territory", None)
    return territory is None or territory is SimpleBed.territory or not isinstance(feature, BedLike)


def _territory(feature: Any) -> tuple[ReferenceSpan, ...] | None:
    """Return the spans of a feature, or None if the feature is its own only span."""
    return None if _is_own_span(feature) else tuple(feature.territory())


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

    Queries must have the same properties, so a BED2 record can be added but cannot be a query.

    A feature overlaps or encloses a query when any of its spans on the query's reference does, so
    a BEDPE record is found by either end, but it is enclosed by a query only when all of its spans
    are. A query yields each matching feature once, in the order the index finds it, which is
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
        self._refname_to_index: dict[Refname, _RefIndex] = {}
        self._holds_spans: bool = False
        self._span_refnames: list[Refname] = []
        self._span_starts: array[int] = array("q")
        self._span_ends: array[int] = array("q")
        self._span_features: array[int] = array("q")
        self._span_strands: list[BedStrand | None] = []
        self._first_spans: array[int] = array("q", [0])
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
        """Add features to this overlap detector, indexing every span of each.

        Raises:
            ValueError: If a feature has no spans, in which case none of the features are added.
        """
        if not self._holds_spans and all(map(_is_own_span, features)):
            self._add_features(features)
            return
        territories = [_territory(feature) for feature in features]
        for feature, territory in zip(features, territories, strict=True):
            if territory == ():
                raise ValueError(f"A feature must have at least one span: {feature!r}")
        if not self._holds_spans:
            self._hold_spans()
        self._add_spans(features, territories)

    def _add_features(self, features: Iterable[Any]) -> None:
        """Add features that are each their own only span, holding the features in the trees."""
        added = self._features
        refname_to_index = self._refname_to_index
        for feature in features:
            start, end = _closed(feature)
            index = refname_to_index.get(feature.refname)
            if index is None:
                index = refname_to_index[feature.refname] = _RefIndex()
            index.tree.add(start, end, feature)
            index.is_built = False  # mark that this tree needs re-indexing
            added.append(feature)

    def _add_spans(
        self, features: Iterable[Any], territories: Iterable[tuple[ReferenceSpan, ...] | None]
    ) -> None:
        """Add features by every span, holding the number of each span in the trees."""
        added = self._features
        span_refnames = self._span_refnames
        span_starts = self._span_starts
        span_ends = self._span_ends
        span_features = self._span_features
        span_strands = self._span_strands
        first_spans = self._first_spans
        refname_to_index = self._refname_to_index
        for feature, territory in zip(features, territories, strict=True):
            number = len(added)
            strand = _strand(feature)
            for span in (feature,) if territory is None else territory:
                start, end = _closed(span)
                index = refname_to_index.get(span.refname)
                if index is None:
                    index = refname_to_index[span.refname] = _RefIndex()
                index.tree.add(start, end, len(span_refnames))
                index.is_built = False  # mark that this tree needs re-indexing
                span_refnames.append(span.refname)
                span_starts.append(span.start)
                span_ends.append(span.end)
                span_features.append(number)
                span_strands.append(_strand(span) or strand)
            first_spans.append(len(span_refnames))
            added.append(feature)

    def _hold_spans(self) -> None:
        """Rebuild the trees to hold span numbers instead of features, which have one span each."""
        features = self._features
        self._features = []
        self._refname_to_index = {}
        self._holds_spans = True
        self._add_spans(features, [None] * len(features))

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

    def _span_hits(self, feature: ReferenceSpan, stranded: bool) -> list[int]:
        """Return the spans overlapping a query feature, on its strand if stranded."""
        tree = self._tree_for(feature.refname)
        if tree is None:
            return []

        hits: list[int] = tree.search_values(*_closed(feature))
        if not stranded:
            return hits

        strand = _strand(feature)
        if strand is None:
            return []
        span_strands = self._span_strands
        return [hit for hit in hits if span_strands[hit] is strand]

    def _features_of(self, spans: list[int]) -> list[FeatureType]:
        """Return the features that own the given spans, once each, in the order of the spans."""
        numbers = dict.fromkeys(map(self._span_features.__getitem__, spans))
        return list(map(self._features.__getitem__, numbers))

    def _is_enclosed(self, number: int, feature: ReferenceSpan, strand: BedStrand | None) -> bool:
        """Return whether every span of a feature is inside a query, and on a strand if given."""
        span_refnames = self._span_refnames
        span_starts = self._span_starts
        span_ends = self._span_ends
        span_strands = self._span_strands
        for span in range(self._first_spans[number], self._first_spans[number + 1]):
            if (
                span_refnames[span] != feature.refname
                or span_starts[span] < feature.start
                or span_ends[span] > feature.end
                or (strand is not None and span_strands[span] is not strand)
            ):
                return False
        return True

    def overlapping(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[FeatureType]:
        """Yields all the features with a span that overlaps the given query feature."""
        if self._holds_spans:
            yield from self._features_of(self._span_hits(feature, stranded))
            return

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
            return tree.has_overlaps(start, end)

        strand = _strand(feature)
        if strand is None:
            return False
        if not self._holds_spans:
            return any(_strand(hit) is strand for hit in tree.iter_values(start, end))
        span_strands = self._span_strands
        return any(span_strands[hit] is strand for hit in tree.iter_values(start, end))

    def enclosing(self, feature: ReferenceSpan, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with a span that completely encloses the given query feature."""
        if not self._holds_spans:
            hits: Iterator[Any] = self.overlapping(feature, stranded=stranded)
            for hit in hits:
                if feature.start >= hit.start and feature.end <= hit.end:
                    yield hit
            return

        starts, ends = self._span_starts, self._span_ends
        yield from self._features_of([
            hit
            for hit in self._span_hits(feature, stranded)
            if feature.start >= starts[hit] and feature.end <= ends[hit]
        ])

    def enclosed_by(
        self, feature: ReferenceSpan, *, stranded: bool = False
    ) -> Iterator[FeatureType]:
        """Yields all the features with every span enclosed by the given query feature."""
        if not self._holds_spans:
            hits: Iterator[Any] = self.overlapping(feature, stranded=stranded)
            for hit in hits:
                if feature.start <= hit.start and feature.end >= hit.end:
                    yield hit
            return

        starts, ends = self._span_starts, self._span_ends
        span_features = self._span_features
        features = self._features
        strand = _strand(feature) if stranded else None
        checked: set[int] = set()
        for hit in self._span_hits(feature, stranded):
            if feature.start <= starts[hit] and feature.end >= ends[hit]:
                number = span_features[hit]
                if number not in checked:
                    checked.add(number)
                    if self._is_enclosed(number, feature, strand):
                        yield features[number]
