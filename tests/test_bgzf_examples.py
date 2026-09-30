from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import is_dataclass
from dataclasses import replace
from io import StringIO
from pathlib import Path
from typing import Any

import pytest
from pybgzf import IndexedReader
from pybgzf import IndexFormat

import bedspec
from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed3N
from bedspec import Bed4
from bedspec import Bed4N
from bedspec import Bed5
from bedspec import Bed5N
from bedspec import Bed6
from bedspec import Bed6N
from bedspec import Bed9
from bedspec import Bed9N
from bedspec import Bed12
from bedspec import Bed12N
from bedspec import BedColor
from bedspec import BedGraph
from bedspec import BedLike
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import BroadPeak
from bedspec import GappedPeak
from bedspec import NarrowPeak
from bedspec import PairBed
from bedspec import PointBed
from bedspec import SimpleBed
from bedspec.overlap import TabixDetector

PLUS = BedStrand.Positive
MINUS = BedStrand.Negative

Build = Callable[[str, int, int, BedStrand | None], Any]
"""Build a record of one BED type from its reference, start, end, and strand."""

RED = BedColor(255, 0, 0)
PEAK: dict[str, Any] = {"p_value": 3.5, "q_value": -1.0, "signal_value": 12.25}


def blocks(start: int, end: int) -> dict[str, Any]:
    """Return one block covering a feature with any length, or no blocks for a zero-length one."""
    if start == end:
        return {"block_count": None, "block_sizes": None, "block_starts": None}
    return {"block_count": 1, "block_sizes": (end - start,), "block_starts": (0,)}


def bed6(refname: str, start: int, end: int, strand: BedStrand | None) -> dict[str, Any]:
    """Return the fields of a BED6 record."""
    return {
        "end": end,
        "name": "f",
        "refname": refname,
        "score": 500,
        "start": start,
        "strand": strand,
    }


def bed9(refname: str, start: int, end: int, strand: BedStrand | None) -> dict[str, Any]:
    """Return the fields of a BED9 record."""
    return {
        **bed6(refname, start, end, strand),
        "item_rgb": RED,
        "thick_end": end,
        "thick_start": start,
    }


BUILDERS: dict[type[Any], Build] = {
    Bed3: lambda r, s, e, _: Bed3(r, start=s, end=e),
    Bed3N: lambda r, s, e, _: Bed3N(r, start=s, end=e, extra=("x",)),
    Bed4: lambda r, s, e, _: Bed4(r, start=s, end=e, name="f"),
    Bed4N: lambda r, s, e, _: Bed4N(r, start=s, end=e, name="f", extra=("x", "y")),
    Bed5: lambda r, s, e, _: Bed5(r, start=s, end=e, name="f", score=500),
    Bed5N: lambda r, s, e, _: Bed5N(r, start=s, end=e, name="f", score=500, extra=()),
    Bed6: lambda r, s, e, strand: Bed6(**bed6(r, s, e, strand)),
    Bed6N: lambda r, s, e, strand: Bed6N(**bed6(r, s, e, strand), extra=("x",)),
    Bed9: lambda r, s, e, strand: Bed9(**bed9(r, s, e, strand)),
    Bed9N: lambda r, s, e, strand: Bed9N(**bed9(r, s, e, strand), extra=("x",)),
    Bed12: lambda r, s, e, strand: Bed12(**bed9(r, s, e, strand), **blocks(s, e)),
    Bed12N: lambda r, s, e, strand: Bed12N(**bed9(r, s, e, strand), **blocks(s, e), extra=("x",)),
    BedGraph: lambda r, s, e, _: BedGraph(r, start=s, end=e, value=1.5),
    BroadPeak: lambda r, s, e, strand: BroadPeak(**bed6(r, s, e, strand), **PEAK),
    GappedPeak: lambda r, s, e, strand: GappedPeak(**bed9(r, s, e, strand), **blocks(s, e), **PEAK),
    NarrowPeak: lambda r, s, e, strand: NarrowPeak(**bed6(r, s, e, strand), **PEAK, peak=-1),
}
"""How to build a record of each BED type that describes an interval."""

INTERVALS: dict[str, tuple[str, int, int, BedStrand | None]] = {
    "A": ("chr1", 10, 20, PLUS),
    "B": ("chr1", 15, 30, MINUS),
    "C": ("chr1", 25, 25, PLUS),
    "D": ("chr1", 40, 50, None),
    "E": ("chr2", 5, 8, PLUS),
}
"""Features on two references, sorted by start, where C is zero-length and D has no strand."""


@dataclass(frozen=True)
class Query:
    """A query and the features every operation, and the index alone, should find for it.

    Attributes:
        feature: the refname, start, end, and strand of the query.
        overlapping: the features overlapping the query.
        enclosing: the features enclosing the query.
        enclosed_by: the features the query encloses.
        stranded: the overlapping features on the query's strand.
        lines: the features whose lines the index returns, with tabix's rule that a zero-length
            feature at `s` is in `[x, y)` only if `x < s < y`.
    """

    feature: tuple[str, int, int, BedStrand | None]
    overlapping: str
    enclosing: str
    enclosed_by: str
    stranded: str
    lines: str


INTERVAL_QUERIES = [
    Query(("chr1", 12, 14, PLUS), "A", "A", "", "A", "A"),
    Query(("chr1", 18, 22, MINUS), "AB", "B", "", "B", "AB"),
    Query(("chr1", 20, 30, PLUS), "BC", "B", "C", "C", "BC"),
    Query(("chr1", 24, 25, PLUS), "BC", "B", "C", "C", "B"),
    Query(("chr1", 25, 25, PLUS), "BC", "BC", "C", "C", ""),
    Query(("chr1", 26, 27, MINUS), "B", "B", "", "B", "B"),
    Query(("chr1", 30, 40, PLUS), "", "", "", "", ""),
    Query(("chr1", 35, 60, None), "D", "", "D", "", "D"),
    Query(("chr2", 0, 100, PLUS), "E", "", "E", "E", "E"),
    Query(("chr3", 0, 100, PLUS), "", "", "", "", ""),
]
"""Queries of `INTERVALS`: a zero-length feature covers the bases on either side of it."""

POINTS: dict[str, Bed2] = {
    "P": Bed2("chr1", start=10),
    "Q": Bed2("chr1", start=11),
    "R": Bed2("chr2", start=0),
}
"""BED2 points, where a point at `s` covers the one base `[s, s + 1)`."""

POINT_QUERIES = [
    Query(("chr1", 9, 10, PLUS), "", "", "", "", ""),
    Query(("chr1", 10, 11, PLUS), "P", "P", "P", "", "P"),
    Query(("chr1", 11, 12, PLUS), "Q", "Q", "Q", "", "Q"),
    Query(("chr1", 10, 12, PLUS), "PQ", "", "PQ", "", "PQ"),
    Query(("chr1", 11, 11, PLUS), "PQ", "PQ", "", "", ""),
    Query(("chr1", 12, 20, PLUS), "", "", "", "", ""),
    Query(("chr2", 0, 1, PLUS), "R", "R", "R", "", "R"),
]
"""Queries of `POINTS`, which have no strand, so a stranded query finds none of them."""

CASES: dict[type[Any], tuple[dict[str, Any], list[Query]]] = {
    Bed2: (POINTS, POINT_QUERIES),
    **{
        record_type: (
            {label: build(*interval) for label, interval in INTERVALS.items()},
            [
                query
                if hasattr(build("chr1", 0, 1, PLUS), "strand")
                else replace(query, stranded="")
                for query in INTERVAL_QUERIES
            ],
        )
        for record_type, build in BUILDERS.items()
    },
}
"""The records of each BED type that can be indexed, and queries of them.

A BED type without a strand never matches a stranded query.
"""


def test_every_indexable_bed_type_has_examples() -> None:
    """Test that every BED type bedspec ships, except BEDPE, has examples here."""
    shipped = {
        value
        for name in bedspec.__all__
        if isinstance(value := getattr(bedspec, name), type)
        and issubclass(value, BedLike)
        and is_dataclass(value)
        and value not in {BedLike, PairBed, PointBed, SimpleBed}
    }
    assert shipped == {*CASES, BedPE}


@pytest.mark.parametrize("index", [IndexFormat.CSI, IndexFormat.TBI], ids=["CSI", "TBI"])
@pytest.mark.parametrize("record_type", list(CASES), ids=lambda record_type: record_type.__name__)
def test_every_bed_type_is_indexed_and_queried(
    record_type: type[Any], index: IndexFormat, tmp_path: Path
) -> None:
    """Test that each BED type is written to indexed BGZF and every query finds what it should."""
    records, queries = CASES[record_type]
    path = tmp_path / "features.bed.gz"
    with BedWriter.from_path[record_type](path, index=index) as writer:
        writer.write_comment("examples")
        for record in records.values():
            writer.write(record)

    def found(results: Any) -> str:
        hits = list(results)
        return "".join(sorted(label for label, record in records.items() if record in hits))

    with TabixDetector[record_type](path) as detector:  # type: ignore[valid-type]
        for query in queries:
            refname, start, end, strand = query.feature
            feature = Bed6(refname, start=start, end=end, name=None, score=None, strand=strand)
            assert found(detector.overlapping(feature)) == query.overlapping, query
            assert detector.overlaps(feature) == bool(query.overlapping), query
            assert found(detector.enclosing(feature)) == query.enclosing, query
            assert found(detector.enclosed_by(feature)) == query.enclosed_by, query
            assert found(detector.overlapping(feature, stranded=True)) == query.stranded, query
            assert detector.overlaps(feature, stranded=True) == bool(query.stranded), query

    lines = text_lines(record_type, records)
    with IndexedReader(path) as reader:
        for query in queries:
            refname, start, end, _ = query.feature
            expected = [lines[label] for label in query.lines]
            assert list(reader.query(refname, start, end)) == expected, query


def text_lines(record_type: type[Any], records: dict[str, Any]) -> dict[str, str]:
    """Return the line a BED writer writes for each record, without its line break."""
    lines: dict[str, str] = {}
    for label, record in records.items():
        handle = StringIO()
        BedWriter[record_type](handle).write(record)  # type: ignore[valid-type]
        lines[label] = handle.getvalue().rstrip("\r\n")
    return lines
