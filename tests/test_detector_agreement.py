import math
import random
from collections import Counter
from collections.abc import Iterable
from collections.abc import Iterator
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import fields
from pathlib import Path
from typing import Any

import pytest
from pybgzf import IndexFormat
from typing_extensions import override

from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import ReferenceSpan
from bedspec.overlap import TabixDetector
from bedspec.overlap import TreeDetector

STRANDS = [BedStrand.Positive, BedStrand.Negative, None]
METHODS = ["overlapping", "enclosing", "enclosed_by"]

SHIFTS = [14, 17, 20, 23, 26, 29]
"""The widths, as powers of two, of the bins of a tabix or CSI index."""

LONGEST = 10_000_000
"""The longest feature, or query, made at random."""

LARGEST = 2**31 - 1
"""The largest position the in-memory detector stores, as a 32-bit integer."""

REFERENCES: dict[str, tuple[int, int]] = {
    "chr1": (2**29, 8000),
    "chr10": (2**27 + 5, 800),
    "chr2": (2**29, 20),
    "chrM": (16_569, 2),
}
"""The length of each reference, at most what tabix allows, and how many features it scatters."""

CSI_REFERENCES: dict[str, tuple[int, int]] = {**REFERENCES, "chrBig": (LARGEST - 1, 1500)}
"""The references of a CSI index, which also allows one longer than tabix allows."""

ABSENT = ["chr", "chrX"]
"""References no feature is written on."""


@dataclass(frozen=True)
class Blocks6(Bed6):
    """A BED6 whose territory is its first and last thirds and a point between them."""

    @override
    def territory(self) -> Iterator[ReferenceSpan]:
        """Yield the first and last thirds, and a zero-length middle on the other strand."""
        third = len(self) // 3
        if third == 0:
            yield self
            return
        middle = self.start + len(self) // 2
        strand = None if self.strand is None else self.strand.opposite()
        yield Bed3(self.refname, start=self.start, end=self.start + third)
        yield Bed6(self.refname, start=middle, end=middle, name=None, score=None, strand=strand)
        yield Bed3(self.refname, start=self.end - third, end=self.end)


def near_boundary(rng: random.Random, length: int) -> int:
    """Return a bin boundary inside a reference, or a position beside one."""
    shift = rng.choice([shift for shift in SHIFTS if 1 << shift < length])
    return (rng.randint(1, (length - 1) >> shift) << shift) + rng.choice([-1, 0, 1])


def log_length(rng: random.Random) -> int:
    """Return a length from 1 to about 10 Mb, log-uniformly, or zero one time in twenty."""
    return 0 if rng.random() < 0.05 else int(math.exp(rng.uniform(0, math.log(LONGEST))))


def feature_at(rng: random.Random, refname: str, start: int, end: int) -> Bed6:
    """Return a feature on a random strand, or none, with a long name to fill BGZF blocks."""
    name = rng.choice([None, "x" * rng.randint(100, 255)])
    score = rng.choice([None, 0, 1000])
    return Bed6(refname, start=start, end=end, name=name, score=score, strand=rng.choice(STRANDS))


def scattered(rng: random.Random, refname: str, length: int, count: int) -> Iterator[Bed6]:
    """Yield features starting anywhere or on bin boundaries, some spanning many bins."""
    for _ in range(count):
        start = near_boundary(rng, length) if rng.random() < 0.5 else rng.randint(1, length)
        yield feature_at(rng, refname, start, min(start + log_length(rng), length))


def nested(rng: random.Random, refname: str, length: int) -> Iterator[Bed6]:
    """Yield a long feature and a chain nested in it, sharing its start, its end, or neither."""
    span = int(math.exp(rng.uniform(math.log(1000), math.log(LONGEST))))
    start = rng.choice([near_boundary(rng, length - span), rng.randint(1, length - span)])
    end = start + span
    yield feature_at(rng, refname, start, end)
    depth = rng.randint(2, 100)
    step = span // (4 * depth)
    for level in range(1, depth):
        inset = level * step
        inner_start, inner_end = rng.choice([
            (start + inset, end - inset),
            (start, end - inset),
            (start + inset, end),
            (start + inset, start + inset + level),
        ])
        yield feature_at(rng, refname, inner_start, inner_end)


def touching(rng: random.Random, refname: str, length: int) -> Iterator[Bed6]:
    """Yield a run of features, each starting where the last ends, across a bin boundary."""
    start = near_boundary(rng, length) - rng.randint(0, 100)
    for _ in range(rng.randint(2, 30)):
        end = min(length, start + rng.choice([0, 1, 2, rng.randint(3, 1000)]))
        yield feature_at(rng, refname, start, end)
        start = end


def hostile_features(rng: random.Random, refname: str, length: int, count: int) -> list[Bed6]:
    """Return scattered, nested, touching, zero-length, and duplicate features, sorted by start."""
    insertions = [near_boundary(rng, length) for _ in range(count // 20)]
    features = [
        *scattered(rng, refname, length, count),
        *(inner for _ in range(count // 100) for inner in nested(rng, refname, length)),
        *(adjacent for _ in range(count // 100) for adjacent in touching(rng, refname, length)),
        *(feature_at(rng, refname, position, position) for position in insertions),
        *(feature_at(rng, refname, start, end) for start, end in [(0, 1), (1, 1), (0, length)]),
        *(feature_at(rng, refname, start, length) for start in [length - 1, length]),
    ]
    features.extend(rng.choices(features, k=len(features) // 20))
    return sorted(features, key=lambda feature: feature.start)


def query_at(rng: random.Random, refname: str, start: int, end: int) -> Bed6:
    """Return a query on a random strand, or none, ending by the largest position stored."""
    strand = rng.choice(STRANDS)
    return Bed6(refname, start=start, end=min(end, LARGEST), name=None, score=None, strand=strand)


def edge_queries(rng: random.Random, feature: ReferenceSpan) -> Iterator[Bed6]:
    """Yield zero-length or one-base queries at each edge of a feature, it, and one base wider."""
    edges = {feature.start - 1, feature.start, feature.end - 1, feature.end, feature.end + 1}
    for position in sorted(edges - {-1}):
        yield query_at(rng, feature.refname, position, position + rng.choice([0, 1]))
    yield query_at(rng, feature.refname, feature.start, feature.end)
    yield query_at(rng, feature.refname, max(feature.start - 1, 0), feature.end + 1)


def hostile_queries(
    rng: random.Random, features: Sequence[Bed6], lengths: dict[str, int]
) -> Iterator[Bed6]:
    """Yield queries at feature edges, inside long features, on bin boundaries, and at random."""
    for feature in rng.sample(features, 80):
        yield from edge_queries(rng, feature)
    for long in rng.sample([feature for feature in features if len(feature) > 1000], 80):
        position = rng.randint(long.start, long.end)
        yield query_at(rng, long.refname, position, position + rng.choice([0, 1, 50]))
    for refname, length in lengths.items():
        for _ in range(30):
            position = near_boundary(rng, length)
            yield query_at(rng, refname, position, position + rng.choice([0, 1, 2]))
        for _ in range(30):
            start = rng.randint(0, length)
            yield query_at(rng, refname, start, start + log_length(rng))
        ends = [(0, 0), (0, length), (length, length), (length - 1, length + 1), (length, LARGEST)]
        for start, end in ends:
            yield query_at(rng, refname, start, end)
    for refname in ABSENT:
        yield query_at(rng, refname, 0, 0)
        yield query_at(rng, refname, 0, LARGEST)


def bgzf_blocks(path: Path) -> int:
    """Return the number of BGZF blocks in a file, reading each block's size from its header."""
    data = path.read_bytes()
    offset = count = 0
    while offset < len(data):
        offset += int.from_bytes(data[offset + 16 : offset + 18], "little") + 1
        count += 1
    return count


def assert_detectors_agree(
    tree: TreeDetector[Any], tabix: TabixDetector[Any], queries: Iterable[Bed6], index: str
) -> None:
    """Assert that both detectors find the same features for every query, method, and strand."""
    for query in queries:
        for stranded in (False, True):
            for method in METHODS:
                expected: Counter[Any] = Counter(getattr(tree, method)(query, stranded=stranded))
                actual: Counter[Any] = Counter(getattr(tabix, method)(query, stranded=stranded))
                assert actual == expected, f"{method}({query}, {stranded=}) with {index}"
            found = bool(list(tree.overlapping(query, stranded=stranded)))
            for detector in (tree, tabix):
                name = type(detector).__name__
                assert detector.overlaps(query, stranded=stranded) is found, (
                    f"{name}.overlaps({query}, {stranded=}) with {index}"
                )


@pytest.mark.parametrize("index", list(IndexFormat))
def test_both_detectors_agree_on_hostile_features(index: IndexFormat, tmp_path: Path) -> None:
    """Test that both detectors agree on nested, duplicate, and touching features across bins."""
    rng = random.Random(f"hostile {index.name}")
    references = CSI_REFERENCES if index is IndexFormat.CSI else REFERENCES
    features = [
        feature
        for refname, (length, count) in references.items()
        for feature in hostile_features(rng, refname, length, count)
    ]
    path = tmp_path / "hostile.bed.gz"
    with BedWriter.from_path[Bed6](path, index=index, threads=4) as writer:
        for feature in features:
            writer.write(feature)
    assert bgzf_blocks(path) > 25

    lengths = {refname: length for refname, (length, _) in references.items()}
    queries = list(hostile_queries(rng, features, lengths))
    with TabixDetector[Bed6](path) as tabix:
        assert_detectors_agree(TreeDetector(features), tabix, queries, index.name)


@pytest.mark.parametrize("index", list(IndexFormat))
def test_both_detectors_agree_on_every_span_of_features_with_several(
    index: IndexFormat, tmp_path: Path
) -> None:
    """Test that both detectors agree on features whose territory has several spans."""
    rng = random.Random(f"blocks {index.name}")
    features = [
        Blocks6(**{field.name: getattr(feature, field.name) for field in fields(Bed6)})
        for refname, (length, count) in REFERENCES.items()
        for feature in hostile_features(rng, refname, length, count // 10)
    ]
    path = tmp_path / "blocks.bed.gz"
    with BedWriter.from_path[Blocks6](path, index=index) as writer:
        for feature in features:
            writer.write(feature)

    lengths = {refname: length for refname, (length, _) in REFERENCES.items()}
    queries = [
        *hostile_queries(rng, features, lengths),
        *(
            query
            for feature in rng.sample(features, 100)
            for span in feature.territory()
            for query in edge_queries(rng, span)
        ),
    ]
    with TabixDetector[Blocks6](path) as tabix:
        assert_detectors_agree(TreeDetector(features), tabix, queries, index.name)
