from array import array
from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import dataclass
from dataclasses import field
from itertools import chain
from typing import Any
from typing import Generic
from typing import TypeAlias
from typing import TypeVar
from typing import cast

from superintervals import IntervalMap
from typing_extensions import override

from bedspec._bedspec import _NATIVE_LOCATE  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import BedLike
from bedspec._bedspec import BedStrand
from bedspec._bedspec import Locatable
from bedspec._bedspec import Locate
from bedspec._bedspec import ReferenceSpan
from bedspec._bedspec import SimpleBed
from bedspec._bedspec import locate

FeatureType = TypeVar("FeatureType", bound=BedLike | Locatable)
"""Type variable for features stored within the overlap detector."""

Refname: TypeAlias = str
"""A type alias for a reference sequence name string."""

IntervalTree: TypeAlias = IntervalMap
"""A type alias for the untyped interval map."""

SpanTest: TypeAlias = Callable[[int, int, ReferenceSpan], bool]
"""A test of a span, by its start and end, against a query on the same reference."""


def _strand(feature: Any) -> BedStrand | None:
    """Return the strand of a feature, or None if it has none."""
    strand: BedStrand | None = getattr(feature, "strand", None)
    return strand


def span_strand(span: ReferenceSpan, strand: BedStrand | None) -> BedStrand | None:
    """Return the strand of a span, or the given strand of its feature if the span has none."""
    span_strand = _strand(span)
    return strand if span_strand is None else span_strand


@dataclass(slots=True)
class _Span:
    """A query read from a feature that names its reference, start or end another way."""

    refname: str
    start: int
    end: int
    strand: BedStrand | None


def as_span(feature: Locatable) -> ReferenceSpan:
    """Return a feature as a span with a refname, start and end: itself, if it has them."""
    getter = locate(feature)
    if getter is _NATIVE_LOCATE:
        return cast(ReferenceSpan, feature)
    refname, start, end = getter(feature)
    return _Span(refname, start, end, _strand(feature))


NOT_NATIVE = (AttributeError, KeyError)
"""What reading `refname`, `start` or `end` of a feature raises when it names them otherwise.

Each method reads them first and falls back to `locate()` only on these, so a BED record pays
nothing; pysam's tabix records raise KeyError for an attribute they do not have.
"""


def _closed(start: int, end: int) -> tuple[int, int]:
    """Return the closed interval of bases a span covers, or flanks if it is zero-length."""
    if start == end:
        return max(start - 1, 0), start
    return start, end - 1


def touches(start: int, end: int, query: ReferenceSpan) -> bool:
    """Return whether a span and a query share a base, a zero-length one by either base beside it.

    This is the test the interval trees answer, since they hold the closed interval of each span.
    """
    first, last = _closed(start, end)
    query_first, query_last = _closed(query.start, query.end)
    return first <= query_last and query_first <= last


def encloses(start: int, end: int, query: ReferenceSpan) -> bool:
    """Return whether a span encloses a query."""
    return start <= query.start and query.end <= end


def is_enclosed_by(start: int, end: int, query: ReferenceSpan) -> bool:
    """Return whether a span is enclosed by a query."""
    return query.start <= start and end <= query.end


def span_matches(
    refname: Refname,
    start: int,
    end: int,
    strand: BedStrand | None,
    query: ReferenceSpan,
    test: SpanTest,
    required_strand: BedStrand | None,
) -> bool:
    """Return whether a span is on the query's reference and required strand and passes a test.

    Both detectors match every span with this, so a stranded query requires the query's strand and
    an unstranded one requires none.
    """
    return (
        refname == query.refname
        and (required_strand is None or strand is required_strand)
        and test(start, end, query)
    )


def _is_own_span(feature: Any) -> bool:
    """Return whether a feature is its own only span, as a record that is not BED-like is."""
    spans = getattr(feature.__class__, "spans", None)
    return spans is None or spans is SimpleBed.spans or not isinstance(feature, BedLike)


def _spans(feature: Any) -> tuple[ReferenceSpan, ...] | None:
    """Return the spans of a feature, or None if the feature is its own only span."""
    return None if _is_own_span(feature) else tuple(feature.spans())


@dataclass(slots=True)
class _RefIndex:
    """Per-reference interval index state."""

    tree: IntervalTree = field(default_factory=IntervalMap)
    is_built: bool = False


class TreeDetector(Iterable[FeatureType], Generic[FeatureType]):
    """Detects and returns overlaps between a collection of reference features and query feature.

    The overlap detector may be built with any BED record, indexed by every one of its `spans()`,
    or with any feature-like Python object that has the following properties, under the first of
    the names given that it has:

      * `refname`, `contig`, `chrom` or `reference_name`: The reference sequence name
      * `start` or `reference_start`: A 0-based start position
      * `end`, `stop` or `reference_end`: A 0-based half-open end position

    So pysam's records are added and queried as they are. The names must hold 0-based half-open
    coordinates: gffutils' 1-based closed `Feature` is read wrongly without an error.

    Queries must have the same properties, so a BED2 record can be added but cannot be a query.

    A feature overlaps or encloses a query when any of its spans on the query's reference does, so
    a BEDPE record is found by either end, but it is enclosed by a query only when all of its spans
    are. A BED12 record's spans are its blocks, so a query that touches only an intron does not
    find it; to find a BED12 record by its whole span, add it as a BED6 record. A query yields
    each matching feature once, in the order the index finds it, which is repeatable but is not
    the order the features were added.

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
        spans = [_spans(feature) for feature in features]
        for feature, feature_spans in zip(features, spans, strict=True):
            if feature_spans == ():
                raise ValueError(f"A feature must have at least one span: {feature!r}")
        if not self._holds_spans:
            self._hold_spans()
        self._add_spans(features, spans)

    def _add_features(self, features: Iterable[Any]) -> None:
        """Add features that are each their own only span, holding the features in the trees."""
        added = self._features
        refname_to_index = self._refname_to_index
        remaining = iter(features)
        for feature in remaining:
            try:
                refname, feature_start, feature_end = feature.refname, feature.start, feature.end
            except NOT_NATIVE:
                self._add_located(chain((feature,), remaining))
                return
            start, end = _closed(feature_start, feature_end)
            index = refname_to_index.get(refname)
            if index is None:
                index = refname_to_index[refname] = _RefIndex()
            index.tree.add(start, end, feature)
            index.is_built = False  # mark that this tree needs re-indexing
            added.append(feature)

    def _add_located(self, features: Iterable[Any]) -> None:
        """Add features that are each their own only span, reading each class's names once."""
        added = self._features
        refname_to_index = self._refname_to_index
        kind: type[Any] | None = None
        getter: Locate = _NATIVE_LOCATE
        for feature in features:
            feature_kind: type[Any] = feature.__class__
            if feature_kind is not kind:
                kind, getter = feature_kind, locate(feature)
            refname, feature_start, feature_end = getter(feature)
            start, end = _closed(feature_start, feature_end)
            index = refname_to_index.get(refname)
            if index is None:
                index = refname_to_index[refname] = _RefIndex()
            index.tree.add(start, end, feature)
            index.is_built = False  # mark that this tree needs re-indexing
            added.append(feature)

    def _add_spans(
        self, features: Iterable[Any], spans: Iterable[tuple[ReferenceSpan, ...] | None]
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
        for feature, feature_spans in zip(features, spans, strict=True):
            number = len(added)
            strand = _strand(feature)
            for span in (feature,) if feature_spans is None else feature_spans:
                try:
                    refname, span_start, span_end = span.refname, span.start, span.end
                except NOT_NATIVE:
                    refname, span_start, span_end = locate(span)(span)
                start, end = _closed(span_start, span_end)
                index = refname_to_index.get(refname)
                if index is None:
                    index = refname_to_index[refname] = _RefIndex()
                index.tree.add(start, end, len(span_refnames))
                index.is_built = False  # mark that this tree needs re-indexing
                span_refnames.append(refname)
                span_starts.append(span_start)
                span_ends.append(span_end)
                span_features.append(number)
                span_strands.append(span_strand(span, strand))
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

    def _span_hits(self, feature: Locatable, stranded: bool) -> list[int]:
        """Return the spans overlapping a query feature, on its strand if stranded."""
        native: Any = feature
        try:
            refname, start, end = native.refname, native.start, native.end
        except NOT_NATIVE:
            return self._span_hits(as_span(feature), stranded)
        tree = self._tree_for(refname)
        if tree is None:
            return []

        hits: list[int] = tree.search_values(*_closed(start, end))
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
        """Return whether every span of a feature is enclosed by a query, on a strand if given."""
        refnames, starts, ends = self._span_refnames, self._span_starts, self._span_ends
        strands = self._span_strands
        for span in range(self._first_spans[number], self._first_spans[number + 1]):
            if not span_matches(
                refnames[span],
                starts[span],
                ends[span],
                strands[span],
                feature,
                is_enclosed_by,
                strand,
            ):
                return False
        return True

    def overlapping(self, feature: Locatable, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with a span that overlaps the given query feature."""
        if self._holds_spans:
            yield from self._features_of(self._span_hits(feature, stranded))
            return

        native: Any = feature
        try:
            refname, start, end = native.refname, native.start, native.end
        except NOT_NATIVE:
            yield from self.overlapping(as_span(feature), stranded=stranded)
            return
        tree = self._tree_for(refname)
        if tree is None:
            return

        start, end = _closed(start, end)
        if not stranded:
            yield from tree.search_values(start, end)
            return

        strand = _strand(feature)
        if strand is None:
            return
        for hit in tree.search_values(start, end):
            if _strand(hit) is strand:
                yield hit

    def overlaps(self, feature: Locatable, *, stranded: bool = False) -> bool:
        """Determine if a query feature overlaps any other features."""
        native: Any = feature
        try:
            refname, start, end = native.refname, native.start, native.end
        except NOT_NATIVE:
            return self.overlaps(as_span(feature), stranded=stranded)
        tree = self._tree_for(refname)
        if tree is None:
            return False

        start, end = _closed(start, end)
        if not stranded:
            return tree.has_overlaps(start, end)

        strand = _strand(feature)
        if strand is None:
            return False
        if not self._holds_spans:
            return any(_strand(hit) is strand for hit in tree.iter_values(start, end))
        span_strands = self._span_strands
        return any(span_strands[hit] is strand for hit in tree.iter_values(start, end))

    def enclosing(self, feature: Locatable, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with a span that completely encloses the given query feature."""
        query = as_span(feature)
        if not self._holds_spans:
            hits: Iterator[Any] = self.overlapping(query, stranded=stranded)
            for hit in hits:
                _, start, end = locate(hit)(hit)
                if encloses(start, end, query):
                    yield hit
            return

        starts, ends = self._span_starts, self._span_ends
        yield from self._features_of([
            hit
            for hit in self._span_hits(query, stranded)
            if encloses(starts[hit], ends[hit], query)
        ])

    def enclosed_by(self, feature: Locatable, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with every span enclosed by the given query feature."""
        query = as_span(feature)
        if not self._holds_spans:
            hits: Iterator[Any] = self.overlapping(query, stranded=stranded)
            for hit in hits:
                _, start, end = locate(hit)(hit)
                if is_enclosed_by(start, end, query):
                    yield hit
            return

        starts, ends = self._span_starts, self._span_ends
        span_features = self._span_features
        features = self._features
        strand = _strand(query) if stranded else None
        checked: set[int] = set()
        for hit in self._span_hits(query, stranded):
            if is_enclosed_by(starts[hit], ends[hit], query):
                number = span_features[hit]
                if number not in checked:
                    checked.add(number)
                    if self._is_enclosed(number, query, strand):
                        yield features[number]
