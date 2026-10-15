from array import array
from collections import Counter
from collections.abc import Callable
from collections.abc import Hashable
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

from bedspec._bedspec import _REGISTRY  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import _UNREGISTERED  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import BedLike
from bedspec._bedspec import BedStrand
from bedspec._bedspec import Locatable
from bedspec._bedspec import ReferenceSpan
from bedspec._bedspec import SimpleBed
from bedspec._bedspec import ToBed
from bedspec._bedspec import _check_span  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import _registered_for  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import _spans  # pyright: ignore[reportPrivateUsage]

FeatureType = TypeVar("FeatureType", bound=Locatable)
"""Type variable for features stored within the overlap detector."""

HitType = TypeVar("HitType")
"""Type variable for what a query of one span finds."""

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


@dataclass(slots=True)
class _QuerySpan:
    """A span of a query that is not its own only span, with the strand it is matched by."""

    refname: str
    start: int
    end: int
    strand: BedStrand | None


def _is_own_span(feature: Any) -> bool:
    """Return whether a feature is its own only span.

    An unregistered linear BED record is, and so is an unregistered record that is not BED-like
    and has a `refname`, `start` and `end`.
    """
    kind: type = type(feature)
    if _REGISTRY and kind not in _UNREGISTERED and _registered_for(kind) is not None:
        return False
    if getattr(kind, "spans", None) is SimpleBed.spans:
        return True
    if isinstance(feature, BedLike):
        return False
    try:
        _ = feature.refname, feature.start, feature.end
    except (AttributeError, KeyError):
        return False
    return True


Located: TypeAlias = tuple[tuple[ReferenceSpan, ...], BedStrand | None]
"""The spans of a feature that is not its own only span, and the strand those without one take."""


def _located(feature: Any) -> Located | None:
    """Return a feature's spans and the strand its spans without one take, or None if it has none.

    A feature that is its own only span has no spans of its own here, and only the spans of a BED
    record take its strand.
    """
    kind: type = type(feature)
    if _REGISTRY and kind not in _UNREGISTERED:
        to_bed = _registered_for(kind)
        if to_bed is not None:
            return tuple(to_bed(feature)), None
    if getattr(kind, "spans", None) is SimpleBed.spans:
        return None
    if isinstance(feature, BedLike):
        return tuple(feature.spans()), _strand(feature)
    try:
        _ = feature.refname, feature.start, feature.end
    except (AttributeError, KeyError):
        return tuple(_spans(feature)), None
    return None


def _registration(feature: Locatable) -> ToBed | None:
    """Return the function registered for a feature's class, if any."""
    kind: type = type(feature)
    return None if kind in _UNREGISTERED or not _REGISTRY else _registered_for(kind)


def _query_spans(feature: Any, to_bed: ToBed | None = None) -> tuple[ReferenceSpan, ...]:
    """Return the spans of a query that is not its own only span, each with its strand.

    The spans of a query registered with `to_bed` keep only their own strand, as a registered
    feature's do; otherwise a span without a strand takes its query's if the query is a BED record.
    A span is used as it is when it needs no strand from its query and its class is not registered.
    """
    if to_bed is None:
        fallback = _strand(feature) if isinstance(feature, BedLike) else None
        located = _spans(feature)
    else:
        fallback, located = None, to_bed(feature)
    spans: list[ReferenceSpan] = []
    for span in located:
        if (fallback is None or _strand(span) is not None) and _registration(span) is None:
            spans.append(span)
        else:
            strand = span_strand(span, fallback)
            spans.append(_QuerySpan(span.refname, span.start, span.end, strand))
    return tuple(spans)


def _as_spans(feature: Any) -> tuple[ReferenceSpan, ...]:
    """Return the spans of a query, which is its own only span if it has a refname, start and end.

    A BED12 query is matched by its whole extent, as an unregistered query with a `refname`,
    `start` and `end` always is.
    """
    to_bed = _registration(feature)
    if to_bed is None:
        try:
            _ = feature.refname, feature.start, feature.end
        except (AttributeError, KeyError):
            pass
        else:
            return (feature,)
    return _query_spans(feature, to_bed)


def _one_span(feature: Locatable, method: str) -> ReferenceSpan | None:
    """Return the only span of a query, or None if it has none.

    Raises:
        ValueError: If the query has more than one span.
    """
    spans = _as_spans(feature)
    if len(spans) > 1:
        raise ValueError(f"{method}() needs a query of one span, not {len(spans)}: {feature!r}")
    return spans[0] if spans else None


def _once_each(
    hits_of_spans: Iterable[Iterable[HitType]], key: Callable[[HitType], Hashable]
) -> Iterator[HitType]:
    """Yield what the spans of one query find, each once, but as often as any one span finds it."""
    yielded: Counter[Hashable] = Counter()
    for hits in hits_of_spans:
        found: Counter[Hashable] = Counter()
        for hit in hits:
            hit_key = key(hit)
            found[hit_key] += 1
            if found[hit_key] > yielded[hit_key]:
                yielded[hit_key] += 1
                yield hit


@dataclass(slots=True)
class _RefIndex:
    """Per-reference interval index state."""

    tree: IntervalTree = field(default_factory=IntervalMap)
    is_built: bool = False


class TreeDetector(Iterable[FeatureType], Generic[FeatureType]):
    """Detects and returns overlaps between a collection of reference features and query feature.

    The overlap detector may be built with any BED record, indexed by every one of its `spans()`,
    with any record registered with `register_locatable`, indexed by the spans its function
    returns, or with any feature-like Python object that has the following properties:

      * `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based half-open end position

    Queries may be any of these too. A query with a `refname`, `start` and `end` is matched by its
    whole extent unless its class is registered, so a BED12 query is matched by its whole extent.

    A feature overlaps or encloses a query when any of its spans on the query's reference does, so
    a BEDPE record is found by either end, but it is enclosed by a query only when all of its spans
    are. A BED12 record's spans are its blocks, so a query that touches only an intron does not
    find it; to find a BED12 record by its whole span, add it as a BED6 record. A query yields
    each matching feature once, in the order the index finds it, which is repeatable but is not
    the order the features were added. A query of several spans, like a registered read and its
    mate, finds the features any of its spans overlaps; it cannot be asked what encloses it or
    what it encloses.

    A feature with no spans, like an unmapped read, is held but never found, and a query with no
    spans finds nothing.

    A zero-length feature, such as an insertion, overlaps features holding either base beside it.

    Every query may be limited to features on the same strand as the query with `stranded=True`.
    The strand compared is that of each matching span, or, for a BED record, of the record if the
    span has none. A feature or query span without a strand never matches a stranded query.

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

        A feature with no spans, like an unmapped read, is held but never found.

        Raises:
            TypeError: If the spans of a feature cannot be found, in which case none are added.
        """
        if not self._holds_spans and all(map(_is_own_span, features)):
            self._add_features(features)
            return
        located = [_located(feature) for feature in features]
        for spans_and_strand in located:
            if spans_and_strand is not None:
                for span in spans_and_strand[0]:
                    _check_span(span.refname, span.start, span.end)
        if not self._holds_spans:
            self._hold_spans()
        self._add_spans(features, located)

    def _add_features(self, features: Iterable[Any]) -> None:
        """Add features that are each their own only span, holding the features in the trees."""
        added = self._features
        refname_to_index = self._refname_to_index
        for feature in features:
            start, end = _closed(feature.start, feature.end)
            index = refname_to_index.get(feature.refname)
            if index is None:
                index = refname_to_index[feature.refname] = _RefIndex()
            index.tree.add(start, end, feature)
            index.is_built = False  # mark that this tree needs re-indexing
            added.append(feature)

    def _add_spans(self, features: Iterable[Any], located: Iterable[Located | None]) -> None:
        """Add features by every span, holding the number of each span in the trees."""
        added = self._features
        span_refnames = self._span_refnames
        span_starts = self._span_starts
        span_ends = self._span_ends
        span_features = self._span_features
        span_strands = self._span_strands
        first_spans = self._first_spans
        refname_to_index = self._refname_to_index
        for feature, spans_and_strand in zip(features, located, strict=True):
            number = len(added)
            spans, strand = ((feature,), None) if spans_and_strand is None else spans_and_strand
            for span in spans:
                start, end = _closed(span.start, span.end)
                index = refname_to_index.get(span.refname)
                if index is None:
                    index = refname_to_index[span.refname] = _RefIndex()
                index.tree.add(start, end, len(span_refnames))
                index.is_built = False  # mark that this tree needs re-indexing
                span_refnames.append(span.refname)
                span_starts.append(span.start)
                span_ends.append(span.end)
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

    def _span_hits(self, query: ReferenceSpan, stranded: bool) -> list[int]:
        """Return the spans overlapping a query of one span, on its strand if stranded."""
        tree = self._tree_for(query.refname)
        if tree is None:
            return []

        hits: list[int] = tree.search_values(*_closed(query.start, query.end))
        if not stranded:
            return hits

        strand = _strand(query)
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
        query: Any = feature
        if _REGISTRY and (to_bed := _registration(feature)) is not None:
            yield from self._overlapping_any(_query_spans(query, to_bed), stranded)
            return
        try:
            refname: Refname = query.refname
        except (AttributeError, KeyError):
            yield from self._overlapping_any(_query_spans(query), stranded)
            return

        if self._holds_spans:
            try:
                hits = self._span_hits(query, stranded)
            except AttributeError:
                yield from self._overlapping_any(_query_spans(query), stranded)
                return
            yield from self._features_of(hits)
            return

        tree = self._tree_for(refname)
        if tree is None:
            return

        try:
            start, end = _closed(query.start, query.end)
        except AttributeError:
            yield from self._overlapping_any(_query_spans(query), stranded)
            return
        if not stranded:
            yield from tree.search_values(start, end)
            return

        strand = _strand(query)
        if strand is None:
            return
        for hit in tree.search_values(start, end):
            if _strand(hit) is strand:
                yield hit

    def _overlapping_any(
        self, spans: tuple[ReferenceSpan, ...], stranded: bool
    ) -> Iterator[FeatureType]:
        """Yield the features with a span overlapping any span of a query, each once."""
        if self._holds_spans:
            yield from self._features_of([
                hit for span in spans for hit in self._span_hits(span, stranded)
            ])
            return
        yield from _once_each((self.overlapping(span, stranded=stranded) for span in spans), id)

    def overlaps(self, feature: Locatable, *, stranded: bool = False) -> bool:
        """Determine if a query feature overlaps any other features."""
        query: Any = feature
        if _REGISTRY and (to_bed := _registration(feature)) is not None:
            return self._overlaps_any(_query_spans(query, to_bed), stranded)
        try:
            tree = self._tree_for(query.refname)
            start, end = _closed(query.start, query.end)
        except (AttributeError, KeyError):
            return self._overlaps_any(_query_spans(query), stranded)
        if tree is None:
            return False

        if not stranded:
            return tree.has_overlaps(start, end)

        strand = _strand(query)
        if strand is None:
            return False
        if not self._holds_spans:
            return any(_strand(hit) is strand for hit in tree.iter_values(start, end))
        span_strands = self._span_strands
        return any(span_strands[hit] is strand for hit in tree.iter_values(start, end))

    def _overlaps_any(self, spans: tuple[ReferenceSpan, ...], stranded: bool) -> bool:
        """Determine if any span of a query overlaps any other features."""
        if len(spans) == 1:
            return self.overlaps(spans[0], stranded=stranded)
        return any(self.overlaps(span, stranded=stranded) for span in spans)

    def enclosing(self, feature: Locatable, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with a span that completely encloses the given query feature.

        Raises:
            ValueError: If the query has more than one span.
        """
        query = _one_span(feature, "enclosing")
        if query is None:
            return

        if not self._holds_spans:
            hits: Iterator[Any] = self.overlapping(query, stranded=stranded)
            for hit in hits:
                if encloses(hit.start, hit.end, query):
                    yield hit
            return

        starts, ends = self._span_starts, self._span_ends
        yield from self._features_of([
            hit
            for hit in self._span_hits(query, stranded)
            if encloses(starts[hit], ends[hit], query)
        ])

    def enclosed_by(self, feature: Locatable, *, stranded: bool = False) -> Iterator[FeatureType]:
        """Yields all the features with every span enclosed by the given query feature.

        Raises:
            ValueError: If the query has more than one span.
        """
        query = _one_span(feature, "enclosed_by")
        if query is None:
            return

        if not self._holds_spans:
            hits: Iterator[Any] = self.overlapping(query, stranded=stranded)
            for hit in hits:
                if is_enclosed_by(hit.start, hit.end, query):
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
