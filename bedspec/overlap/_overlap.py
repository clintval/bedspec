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
class _Span(Hashable, Protocol):
    """A span with a start and an end. 0-based open-ended."""

    @property
    def start(self) -> int:
        """A 0-based start position."""

    @property
    def end(self) -> int:
        """A 0-based open-ended position."""


@runtime_checkable
class _GenomicSpanWithChrom(_Span, Protocol):
    """A genomic feature where reference sequence is accessed with `chrom`."""

    @property
    def chrom(self) -> str:
        """A reference sequence name."""


@runtime_checkable
class _GenomicSpanWithContig(_Span, Protocol):
    """A genomic feature where reference sequence is accessed with `contig`."""

    @property
    def contig(self) -> str:
        """A reference sequence name."""


@runtime_checkable
class _GenomicSpanWithRefName(_Span, Protocol):
    """A genomic feature where reference sequence is accessed with `refname`."""

    @property
    def refname(self) -> str:
        """A reference sequence name."""


GenomicSpanLike = TypeVar(
    "GenomicSpanLike",
    bound=_GenomicSpanWithChrom | _GenomicSpanWithContig | _GenomicSpanWithRefName,
)
"""
A 0-based exclusive genomic feature where the reference sequence name is accessed with any of the 3
most common property names ("chrom", "contig", "refname").
"""

_GenericGenomicSpanLike = TypeVar(
    "_GenericGenomicSpanLike",
    bound=_GenomicSpanWithChrom | _GenomicSpanWithContig | _GenomicSpanWithRefName,
)
"""
A generic 0-based exclusive genomic feature where the reference sequence name is accessed with any
of the most common property names ("chrom", "contig", "refname"). This type variable is used for
describing the generic type contained within the :class:`~bedspec.overlap.OverlapDetector`.
"""

Refname: TypeAlias = str
"""A type alias for a reference sequence name string."""


class OverlapDetector(Generic[_GenericGenomicSpanLike], Iterable[_GenericGenomicSpanLike]):
    """Detects and returns overlaps between a collection of genomic features and an interval.

    The overlap detector may be built with any genomic feature-like Python object that has the
    following properties:

      * `chrom` or `contig` or `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based exclusive end position

    This detector is most efficiently used when all features to be queried are added ahead of time.

    """

    def __init__(self, features: Iterable[_GenericGenomicSpanLike] | None = None) -> None:
        self._refname_to_features: dict[Refname, list[_GenericGenomicSpanLike]] = defaultdict(list)
        self._refname_to_tree: dict[Refname, cr.cgranges] = defaultdict(cr.cgranges)  # type: ignore[attr-defined,name-defined]
        self._refname_to_is_indexed: dict[Refname, bool] = defaultdict(lambda: False)
        if features is not None:
            self.add_all(features)

    def __iter__(self) -> Iterator[_GenericGenomicSpanLike]:
        """Iterate over the features in the overlap detector."""
        return chain(*self._refname_to_features.values())

    @staticmethod
    def _reference_sequence_name(feature: GenomicSpanLike) -> Refname:
        """Return the reference name of a given genomic feature."""
        if isinstance(feature, _GenomicSpanWithContig):
            return feature.contig
        if isinstance(feature, _GenomicSpanWithChrom):
            return feature.chrom
        elif isinstance(feature, _GenomicSpanWithRefName):
            return feature.refname
        else:
            raise ValueError(
                f"Genomic feature is missing a reference sequence name property: {feature}"
            )

    def add(self, feature: _GenericGenomicSpanLike) -> None:
        """Add a genomic feature to this overlap detector."""
        if not isinstance(feature, Hashable):
            raise ValueError(f"Genomic feature is not hashable but should be: {feature}")

        refname: Refname = self._reference_sequence_name(feature)
        feature_idx: int = len(self._refname_to_features[refname])

        self._refname_to_features[refname].append(feature)
        self._refname_to_tree[refname].add(refname, feature.start, feature.end, feature_idx)
        self._refname_to_is_indexed[refname] = False  # mark that this tree needs re-indexing

    def add_all(self, features: Iterable[_GenericGenomicSpanLike]) -> None:
        """Adds one or more genomic features to this overlap detector."""
        for feature in features:
            self.add(feature)

    def overlapping(self, feature: GenomicSpanLike) -> Iterator[_GenericGenomicSpanLike]:
        """Yields all the overlapping features for a given genomic span."""
        refname: Refname = self._reference_sequence_name(feature)

        if refname in self._refname_to_tree and not self._refname_to_is_indexed[refname]:
            self._refname_to_tree[refname].index()  # index the tree if we find it is not indexed

        for *_, idx in self._refname_to_tree[refname].overlap(refname, feature.start, feature.end):
            yield self._refname_to_features[refname][idx]

    def overlaps_any(self, feature: GenomicSpanLike) -> bool:
        """Determine if a given genomic span overlaps any features."""
        return next(self.overlapping(feature), None) is not None

    def those_enclosing(self, feature: GenomicSpanLike) -> Iterator[_GenericGenomicSpanLike]:
        """Yields all the overlapping features that completely enclose the given genomic span."""
        for overlap in self.overlapping(feature):
            if feature.start >= overlap.start and feature.end <= overlap.end:
                yield overlap

    def those_enclosed_by(self, feature: GenomicSpanLike) -> Iterator[_GenericGenomicSpanLike]:
        """Yields all the overlapping features that are enclosed by the given genomic span."""
        for overlap in self.overlapping(feature):
            if feature.start <= overlap.start and feature.end >= overlap.end:
                yield overlap
