import random
from collections import Counter
from collections.abc import Callable
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pytest
from pybgzf import IndexFormat

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedLike
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import ReferenceSpan
from bedspec.overlap import IndexedOverlapDetector
from bedspec.overlap import OverlapDetector

REFNAMES = ["chr1", "chr2", "chr10"]
STRANDS = [BedStrand.Positive, BedStrand.Negative, None]
METHODS = ["overlapping", "enclosing", "enclosed_by"]


def random_span(rng: random.Random) -> tuple[int, int]:
    """Return a random start and end, zero-length about one time in ten, never at position 0."""
    start = rng.randint(1, 3000)
    return start, start + (0 if rng.random() < 0.1 else rng.randint(1, 200))


def random_bed3(rng: random.Random, refname: str) -> Bed3:
    """Return a random BED3 record."""
    start, end = random_span(rng)
    return Bed3(refname, start=start, end=end)


def random_bed6(rng: random.Random, refname: str) -> Bed6:
    """Return a random BED6 record on a random strand, or none."""
    start, end = random_span(rng)
    return Bed6(
        refname,
        start=start,
        end=end,
        name="x",
        score=rng.randint(0, 1000),
        strand=rng.choice(STRANDS),
    )


def random_bed12(rng: random.Random, refname: str) -> Bed12:
    """Return a random BED12 record, with one block when it has any length."""
    start, end = random_span(rng)
    blocks = end > start and rng.random() < 0.5
    return Bed12(
        refname,
        start=start,
        end=end,
        name=None,
        score=None,
        strand=rng.choice(STRANDS),
        thick_start=start,
        thick_end=end,
        item_rgb=None,
        block_count=1 if blocks else None,
        block_sizes=(end - start,) if blocks else None,
        block_starts=(0,) if blocks else None,
    )


def random_bed2(rng: random.Random, refname: str) -> Bed2:
    """Return a random BED2 record."""
    return Bed2(refname, start=rng.randint(0, 3000))


def random_query(rng: random.Random) -> Bed6:
    """Return a random query, sometimes zero-length or on a reference with no features."""
    start = rng.randint(0, 3200)
    end = start + (0 if rng.random() < 0.2 else rng.randint(1, 300))
    refname = rng.choice([*REFNAMES, "chrX"])
    return Bed6(refname, start=start, end=end, name=None, score=None, strand=rng.choice(STRANDS))


def span(record: BedLike) -> ReferenceSpan:
    """Return the one interval a BED record covers."""
    return next(record.territory())


def write_sorted(
    path: Path, record_type: type[Any], records: Iterable[Any], index: IndexFormat
) -> None:
    """Write records sorted by reference, in the order of `REFNAMES`, then start, with an index."""
    ordered = sorted(records, key=lambda record: (REFNAMES.index(record.refname), record.start))
    with BedWriter.from_path[record_type](path, index=index) as writer:
        writer.write_comment("features")
        for record in ordered:
            writer.write(record)


@pytest.mark.parametrize(
    "record_type,build",
    [(Bed3, random_bed3), (Bed6, random_bed6), (Bed12, random_bed12), (Bed2, random_bed2)],
)
@pytest.mark.parametrize("index", list(IndexFormat))
def test_queries_agree_with_the_overlap_detector(
    record_type: type[Any],
    build: Callable[[random.Random, str], Any],
    index: IndexFormat,
    tmp_path: Path,
) -> None:
    """Test that every query agrees with the in-memory overlap detector on random features."""
    rng = random.Random(42)
    records = [build(rng, rng.choice(REFNAMES)) for _ in range(4000)]
    path = tmp_path / "features.bed.gz"
    write_sorted(path, record_type, records, index)

    by_span: dict[Any, list[Any]] = {}
    for record in records:
        by_span.setdefault(span(record), []).append(record)
    oracle = OverlapDetector(by_span)

    def expected(method: str, query: Bed6, stranded: bool) -> Counter[Any]:
        found: Iterable[ReferenceSpan] = getattr(oracle, method)(query, stranded=False)
        matches = Counter(record for key in found for record in by_span[key])
        if not stranded:
            return matches
        return Counter({
            record: count
            for record, count in matches.items()
            if query.strand is not None and getattr(record, "strand", None) is query.strand
        })

    with IndexedOverlapDetector[record_type](path) as detector:  # type: ignore[valid-type]
        for _ in range(400):
            query = random_query(rng)
            for stranded in (False, True):
                for method in METHODS:
                    actual: Counter[Any] = Counter(
                        getattr(detector, method)(query, stranded=stranded)
                    )
                    assert actual == expected(method, query, stranded), (method, query, stranded)
                assert detector.overlaps(query, stranded=stranded) == bool(
                    expected("overlapping", query, stranded)
                )


def test_stranded_queries_agree_with_the_overlap_detector(tmp_path: Path) -> None:
    """Test that stranded queries find exactly what the in-memory detector finds."""
    rng = random.Random(7)
    records = [random_bed6(rng, rng.choice(REFNAMES)) for _ in range(2000)]
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed6, records, IndexFormat.TBI)
    oracle = OverlapDetector(records)

    with IndexedOverlapDetector[Bed6](path) as detector:
        for _ in range(400):
            query = random_query(rng)
            for method in METHODS:
                actual: Counter[Any] = Counter(getattr(detector, method)(query, stranded=True))
                assert actual == Counter(getattr(oracle, method)(query, stranded=True))


def test_zero_length_features_overlap_the_bases_beside_them(tmp_path: Path) -> None:
    """Test that a zero-length feature is found by queries of either base beside it."""
    insertion = Bed3("chr1", start=5, end=5)
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed3, [insertion], IndexFormat.TBI)

    with IndexedOverlapDetector[Bed3](path) as detector:
        assert list(detector.overlapping(Bed3("chr1", start=4, end=5))) == [insertion]
        assert list(detector.overlapping(Bed3("chr1", start=5, end=6))) == [insertion]
        assert list(detector.overlapping(Bed3("chr1", start=5, end=5))) == [insertion]
        assert list(detector.overlapping(Bed3("chr1", start=3, end=4))) == []
        assert list(detector.overlapping(Bed3("chr1", start=6, end=7))) == []


def test_a_zero_length_feature_at_the_start_of_a_reference_is_never_found(tmp_path: Path) -> None:
    """Test that a zero-length feature at position 0 is not found, as tabix never returns it."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed3, [Bed3("chr1", start=0, end=0)], IndexFormat.TBI)

    with IndexedOverlapDetector[Bed3](path) as detector:
        assert not detector.overlaps(Bed3("chr1", start=0, end=1))


def test_an_index_path_and_threads_may_be_given(tmp_path: Path) -> None:
    """Test that the index may be found at any path and read on several threads."""
    bed = Bed3("chr1", start=1, end=5)
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed3, [bed], IndexFormat.CSI)
    index_path = (tmp_path / "features.bed.gz.csi").rename(tmp_path / "elsewhere.csi")

    with IndexedOverlapDetector[Bed3](path, index_path=index_path, threads=2) as detector:
        assert list(detector.overlapping(bed)) == [bed]


def test_the_file_is_closed_when_the_context_ends(tmp_path: Path) -> None:
    """Test that the file is closed when the context ends."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed3, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    with IndexedOverlapDetector[Bed3](path) as detector:
        assert not detector.closed
    assert detector.closed


def test_a_record_type_is_required(tmp_path: Path) -> None:
    """Test that the detector must be subscripted with a BED type."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed3, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    with pytest.raises(TypeError, match=r"IndexedOverlapDetector\[Bed3\]"):
        with IndexedOverlapDetector(path):
            pass


def test_a_paired_bed_cannot_be_queried(tmp_path: Path) -> None:
    """Test that BEDPE, which has two intervals, is refused."""
    path = tmp_path / "features.bed.gz"
    write_sorted(path, Bed3, [Bed3("chr1", start=1, end=5)], IndexFormat.TBI)

    with pytest.raises(TypeError, match="BedPE"):
        with IndexedOverlapDetector[BedPE](path):  # type: ignore[type-var]  # pyright: ignore[reportInvalidTypeArguments]
            pass
