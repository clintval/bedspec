import gzip
import random
import shutil
import subprocess
from collections import Counter
from collections.abc import Callable
from dataclasses import is_dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import pytest
from pybgzf import IndexedReader
from pybgzf import IndexFormat
from typeline import TsvWriter

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
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import BroadPeak
from bedspec import GappedPeak
from bedspec import NarrowPeak
from bedspec import PairBed
from bedspec import PointBed
from bedspec import SimpleBed
from bedspec.overlap import IndexedOverlapDetector
from bedspec.overlap import OverlapDetector

Fields = dict[str, Any]
"""The fields of a record, by name."""

REFNAMES = ["chr1", "chr2", "chr10"]
METHODS = ["overlapping", "enclosing", "enclosed_by"]
TABIX = shutil.which("tabix")
INDEX_SUFFIXES = {IndexFormat.CSI: ".csi", IndexFormat.TBI: ".tbi"}


def span_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random span, zero-length about one time in ten, never zero-length at 0."""
    start = rng.randint(1, 5000)
    end = start + (0 if rng.random() < 0.1 else rng.randint(1, 300))
    return {"refname": refname, "start": start, "end": end}


def named_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random span with a name, or none."""
    return {**span_fields(rng, refname), "name": rng.choice([None, "gene", "gene a", "peak-7"])}


def scored_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random named span with a score, or none."""
    return {**named_fields(rng, refname), "score": rng.choice([None, 0, rng.randint(0, 1000)])}


def stranded_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random scored span on a strand, or none."""
    strand = rng.choice([BedStrand.Positive, BedStrand.Negative, None])
    return {**scored_fields(rng, refname), "strand": strand}


def thick_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random stranded span with a thick part and a color, or neither."""
    fields = stranded_fields(rng, refname)
    thick = rng.random() < 0.5
    thick_start = rng.randint(fields["start"], fields["end"]) if thick else None
    color = BedColor(rng.randint(0, 255), rng.randint(0, 255), rng.randint(0, 255))
    return {
        **fields,
        "item_rgb": rng.choice([None, color]),
        "thick_end": rng.randint(thick_start, fields["end"]) if thick_start is not None else None,
        "thick_start": thick_start,
    }


def block_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random span with a thick part and one or two blocks, or none."""
    fields = thick_fields(rng, refname)
    length = fields["end"] - fields["start"]
    if length < 2 or rng.random() < 0.3:
        return {**fields, "block_count": None, "block_sizes": None, "block_starts": None}
    if rng.random() < 0.5:
        return {**fields, "block_count": 1, "block_sizes": (length,), "block_starts": (0,)}
    first = rng.randint(1, length - 1)
    second = rng.randint(first, length - 1)
    return {
        **fields,
        "block_count": 2,
        "block_sizes": (first, length - second),
        "block_starts": (0, second),
    }


def with_extra(
    build: Callable[[random.Random, str], Fields],
) -> Callable[[random.Random, str], Fields]:
    """Add zero to three extra columns to the fields a factory builds."""

    def build_with_extra(rng: random.Random, refname: str) -> Fields:
        extra = tuple(rng.choice(["1.5", "x", "a b"]) for _ in range(rng.randint(0, 3)))
        return {**build(rng, refname), "extra": extra}

    return build_with_extra


def with_scores(
    build: Callable[[random.Random, str], Fields],
) -> Callable[[random.Random, str], Fields]:
    """Add an ENCODE signal value, and a p-value and q-value that may be missing as -1."""

    def build_with_scores(rng: random.Random, refname: str) -> Fields:
        return {
            **build(rng, refname),
            "p_value": rng.choice([-1.0, round(rng.uniform(0, 50), 3)]),
            "q_value": rng.choice([-1.0, round(rng.uniform(0, 50), 3)]),
            "signal_value": round(rng.uniform(0, 100), 3),
        }

    return build_with_scores


def narrow_peak_fields(rng: random.Random, refname: str) -> Fields:
    """Return a random narrowPeak, with a summit when it has any length, or -1."""
    fields = with_scores(stranded_fields)(rng, refname)
    length = fields["end"] - fields["start"]
    return {**fields, "peak": rng.choice([-1, rng.randrange(length)]) if length else -1}


FACTORIES: dict[type[Any], Callable[[random.Random, str], Fields]] = {
    Bed2: lambda rng, refname: {"refname": refname, "start": rng.randint(0, 5000)},
    Bed3: span_fields,
    Bed3N: with_extra(span_fields),
    Bed4: named_fields,
    Bed4N: with_extra(named_fields),
    Bed5: scored_fields,
    Bed5N: with_extra(scored_fields),
    Bed6: stranded_fields,
    Bed6N: with_extra(stranded_fields),
    Bed9: thick_fields,
    Bed9N: with_extra(thick_fields),
    Bed12: block_fields,
    Bed12N: with_extra(block_fields),
    BedGraph: lambda rng, refname: {**span_fields(rng, refname), "value": rng.uniform(-5, 5)},
    BroadPeak: with_scores(stranded_fields),
    GappedPeak: with_scores(block_fields),
    NarrowPeak: narrow_peak_fields,
}
"""How to build random fields for each BED type that can be indexed."""


def random_records(record_type: type[Any], rng: random.Random, count: int = 2000) -> list[Any]:
    """Return random records sorted by reference, in the order of `REFNAMES`, then start."""
    build = FACTORIES[record_type]
    records = [record_type(**build(rng, rng.choice(REFNAMES))) for _ in range(count)]
    return sorted(records, key=lambda record: (REFNAMES.index(record.refname), record.start))


def random_query(rng: random.Random) -> Bed6:
    """Return a random query, sometimes zero-length or on a reference with no features."""
    start = rng.randint(0, 5400)
    end = start + (0 if rng.random() < 0.2 else rng.randint(1, 400))
    strand = rng.choice([BedStrand.Positive, BedStrand.Negative, None])
    return Bed6(
        rng.choice([*REFNAMES, "chrX"]), start=start, end=end, name=None, score=None, strand=strand
    )


def write_with_comments(writer: TsvWriter[Any], records: list[Any]) -> None:
    """Write records with a comment before them and another halfway through."""
    writer.write_comment("random features")
    for number, record in enumerate(records):
        if number == len(records) // 2:
            writer.write_comment("# halfway")
        writer.write(record)


def plain_text(record_type: type[Any], records: list[Any]) -> str:
    """Return the text a BED writer writes for records, with the same comments."""
    handle = StringIO()
    write_with_comments(BedWriter[record_type](handle), records)  # type: ignore[valid-type]
    return handle.getvalue()


def oracle(record_type: type[Any], records: list[Any]) -> Callable[[str, Bed6, bool], Counter[Any]]:
    """Return what the in-memory overlap detector finds for a method, query, and strandedness."""
    if issubclass(record_type, PointBed):
        detector = OverlapDetector([
            Bed3(r.refname, start=r.start, end=r.start + 1) for r in records
        ])

        def found(method: str, query: Bed6, stranded: bool) -> Counter[Any]:
            spans = getattr(detector, method)(query, stranded=stranded)
            return Counter(Bed2(span.refname, start=span.start) for span in spans)

        return found

    in_memory = OverlapDetector(records)
    return lambda method, query, stranded: Counter(
        getattr(in_memory, method)(query, stranded=stranded)
    )


def tabix_index(tabix: str, path: Path, record_type: type[Any], index: IndexFormat) -> bytes:
    """Index a copy of a BGZF file with htslib's tabix and return the decompressed index."""
    copy = path.parent / "tabix" / path.name
    copy.parent.mkdir()
    _ = shutil.copyfile(path, copy)
    preset = ["-0", "-s1", "-b2", "-e2"] if issubclass(record_type, PointBed) else ["-p", "bed"]
    csi = ["-C"] if index is IndexFormat.CSI else []
    _ = subprocess.run([tabix, *csi, *preset, str(copy)], check=True, capture_output=True)
    return gzip.decompress(Path(f"{copy}{INDEX_SUFFIXES[index]}").read_bytes())


def tabix_lines(tabix: str, path: Path, refname: str, start: int, end: int) -> list[str]:
    """Return the lines htslib's tabix finds in a 0-based, half-open region."""
    region = f"{refname}:{start + 1}-{end}"
    result = subprocess.run([tabix, str(path), region], check=True, capture_output=True, text=True)
    stdout: str = result.stdout
    return stdout.splitlines()


def test_every_shipped_bed_type_is_covered() -> None:
    """Test that every BED type bedspec ships has a factory here, or is BEDPE, never indexed."""
    bases = {BedLike, PairBed, PointBed, SimpleBed}
    shipped = {
        value
        for name in bedspec.__all__
        if isinstance(value := getattr(bedspec, name), type)
        and issubclass(value, BedLike)
        and is_dataclass(value)
        and value not in bases
    }
    assert shipped == {*FACTORIES, BedPE}


def write_indexed(
    record_type: type[Any], index: IndexFormat, threads: int, path: Path
) -> list[Any]:
    """Write random records of a BED type to an indexed BGZF file and return them."""
    records = random_records(record_type, random.Random(record_type.__name__))
    with BedWriter.from_path[record_type](path, index=index, threads=threads) as writer:
        write_with_comments(writer, records)
    return records


def assert_queries_agree(
    detector: IndexedOverlapDetector[Any], expected: Callable[[str, Bed6, bool], Counter[Any]]
) -> None:
    """Assert that random queries find what the in-memory overlap detector finds."""
    rng = random.Random(0)
    for _ in range(100):
        query = random_query(rng)
        for stranded in (False, True):
            for method in METHODS:
                actual: Counter[Any] = Counter(getattr(detector, method)(query, stranded=stranded))
                assert actual == expected(method, query, stranded), (method, query, stranded)
            found = expected("overlapping", query, stranded)
            assert detector.overlaps(query, stranded=stranded) == bool(found)


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("index", list(IndexFormat))
@pytest.mark.parametrize("record_type", list(FACTORIES), ids=lambda t: t.__name__)
def test_every_bed_type_round_trips_and_is_queried(
    record_type: type[Any], index: IndexFormat, threads: int, tmp_path: Path
) -> None:
    """Test that a BED type round-trips through indexed BGZF and is found by every query."""
    path = tmp_path / "features.bed.gz"
    records = write_indexed(record_type, index, threads, path)

    with BedReader.from_path[record_type](path) as reader:
        assert list(reader) == records
    assert gzip.decompress(path.read_bytes()).decode() == plain_text(record_type, records)
    with IndexedOverlapDetector[record_type](path, threads=threads) as detector:  # type: ignore[valid-type]
        assert_queries_agree(detector, oracle(record_type, records))


@pytest.mark.skipif(TABIX is None, reason="htslib's tabix is not installed")
@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("index", list(IndexFormat))
@pytest.mark.parametrize("record_type", list(FACTORIES), ids=lambda t: t.__name__)
def test_every_bed_type_is_indexed_as_tabix_indexes_it(
    record_type: type[Any], index: IndexFormat, threads: int, tmp_path: Path
) -> None:
    """Test that the index matches tabix's byte for byte, and both find the same lines."""
    assert TABIX is not None
    path = tmp_path / "features.bed.gz"
    records = write_indexed(record_type, index, threads, path)

    written = gzip.decompress(Path(f"{path}{INDEX_SUFFIXES[index]}").read_bytes())
    assert written == tabix_index(TABIX, path, record_type, index)

    rng = random.Random(1)
    with IndexedReader(path) as reader:
        for _ in range(20):
            query = random_query(rng)
            end = max(query.end, query.start + 1)
            lines = list(reader.query(query.refname, query.start, end))
            assert lines == tabix_lines(TABIX, path, query.refname, query.start, end), query

    reference = path.parent / "tabix" / f"{path.name}{INDEX_SUFFIXES[index]}"
    with IndexedOverlapDetector[record_type](path, index_path=reference) as detector:  # type: ignore[valid-type]
        assert_queries_agree(detector, oracle(record_type, records))


def test_bedpe_round_trips_through_bgzf_but_is_not_indexed(tmp_path: Path) -> None:
    """Test that BEDPE is written to BGZF and read back, but an index of it is refused."""
    rng = random.Random(0)
    records = [
        BedPE.from_bed6(
            Bed6(**stranded_fields(rng, rng.choice(REFNAMES))),
            Bed6(**stranded_fields(rng, rng.choice(REFNAMES))),
            name=rng.choice([None, "pair"]),
            score=rng.choice([None, rng.randint(0, 5000)]),
        )
        for _ in range(500)
    ]
    path = tmp_path / "pairs.bedpe.gz"
    with BedWriter.from_path[BedPE](path, threads=4) as writer:
        write_with_comments(writer, records)

    with BedReader.from_path[BedPE](path) as reader:
        assert list(reader) == records
    assert gzip.decompress(path.read_bytes()).decode() == plain_text(BedPE, records)
    with pytest.raises(ValueError, match="BedPE"):
        _ = BedWriter.from_path[BedPE](tmp_path / "indexed.bedpe.gz", index=IndexFormat.TBI)
