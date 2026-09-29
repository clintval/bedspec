from dataclasses import fields
from io import StringIO
from pathlib import Path
from typing import Any

import pytest

from bedspec import Bed6
from bedspec import Bed9
from bedspec import Bed9N
from bedspec import Bed12
from bedspec import BedColor
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import BroadPeak
from bedspec import GappedPeak
from bedspec import NarrowPeak


def roundtrip(record: Any, text: str, tmp_path: Path) -> None:
    """Assert that a record is written as the given text and read back as itself."""
    with BedWriter.from_path[type(record)](tmp_path / "test.bed") as writer:
        writer.write(record)
    assert (tmp_path / "test.bed").read_text() == text
    with BedReader.from_path[type(record)](tmp_path / "test.bed") as reader:
        assert list(reader) == [record]


def test_bed9_roundtrips(tmp_path: Path) -> None:
    """Test that a BED9 record is written and read back."""
    record = Bed9(
        refname="chr1",
        start=2,
        end=10,
        name="bed9",
        score=2,
        strand=BedStrand.Positive,
        thick_start=3,
        thick_end=4,
        item_rgb=BedColor(255, 0, 0),
    )
    roundtrip(record, "chr1\t2\t10\tbed9\t2\t+\t3\t4\t255,0,0\n", tmp_path)


def test_bed9_checks_its_thick_bounds() -> None:
    """Test that a BED9 record checks its thick part like a BED12 record does."""
    with pytest.raises(ValueError, match="thick_start and thick_end must satisfy"):
        Bed9(
            refname="chr1",
            start=2,
            end=10,
            name=None,
            score=None,
            strand=None,
            thick_start=1,
            thick_end=4,
            item_rgb=None,
        )


def test_bed9n_keeps_extra_columns() -> None:
    """Test that a BED9+N record keeps any further columns as text."""
    line = "chr1\t2\t10\tbed9\t2\t+\t3\t4\t0\tx\ty\n"
    with BedReader[Bed9N](StringIO(line)) as reader:
        (record,) = list(reader)
    assert record.item_rgb is None
    assert record.extra == ("x", "y")


def test_narrow_peak_roundtrips(tmp_path: Path) -> None:
    """Test that an ENCODE narrowPeak record is written and read back."""
    record = NarrowPeak(
        refname="chr1",
        start=9356548,
        end=9356648,
        name=None,
        score=0,
        strand=None,
        signal_value=182.0,
        p_value=5.0945,
        q_value=-1.0,
        peak=50,
    )
    roundtrip(record, "chr1\t9356548\t9356648\t.\t0\t.\t182.0\t5.0945\t-1.0\t50\n", tmp_path)


def test_narrow_peak_reads_encode_text() -> None:
    """Test that a narrowPeak line as ENCODE writes it is read, with -1 for values not given."""
    line = "chr1\t9356548\t9356648\t.\t0\t.\t182\t5.0945\t-1\t-1\n"
    with BedReader[NarrowPeak](StringIO(line)) as reader:
        (record,) = list(reader)
    assert (record.signal_value, record.p_value, record.q_value, record.peak) == (
        182.0,
        5.0945,
        -1.0,
        -1,
    )


@pytest.mark.parametrize("peak", [-2, 100])
def test_narrow_peak_must_sit_within_the_feature(peak: int) -> None:
    """Test that a narrowPeak's summit must be -1 or an offset within the feature."""
    with pytest.raises(ValueError, match="peak must be -1 or an offset within the feature!"):
        NarrowPeak(
            refname="chr1",
            start=0,
            end=100,
            name=None,
            score=0,
            strand=None,
            signal_value=1.0,
            p_value=-1.0,
            q_value=-1.0,
            peak=peak,
        )


def test_broad_peak_roundtrips(tmp_path: Path) -> None:
    """Test that an ENCODE broadPeak record is written and read back."""
    record = BroadPeak(
        refname="chr1",
        start=10,
        end=500,
        name="peak1",
        score=1000,
        strand=BedStrand.Negative,
        signal_value=3.5,
        p_value=-1.0,
        q_value=2.25,
    )
    roundtrip(record, "chr1\t10\t500\tpeak1\t1000\t-\t3.5\t-1.0\t2.25\n", tmp_path)


def test_gapped_peak_roundtrips(tmp_path: Path) -> None:
    """Test that an ENCODE gappedPeak record is written and read back."""
    record = GappedPeak(
        refname="chr1",
        start=2,
        end=10,
        name="peak1",
        score=500,
        strand=None,
        thick_start=2,
        thick_end=10,
        item_rgb=None,
        block_count=2,
        block_sizes=(2, 3),
        block_starts=(0, 5),
        signal_value=3.5,
        p_value=4.0,
        q_value=-1.0,
    )
    roundtrip(
        record, "chr1\t2\t10\tpeak1\t500\t.\t2\t10\t0\t2\t2,3\t0,5\t3.5\t4.0\t-1.0\n", tmp_path
    )


@pytest.mark.parametrize(
    "extended,base",
    [(NarrowPeak, Bed6), (BroadPeak, Bed6), (GappedPeak, Bed12), (Bed9N, Bed9)],
)
def test_extended_formats_add_fields_to_their_base(extended: type[Any], base: type[Any]) -> None:
    """Test that each extended format is its base BED type followed by its own fields."""
    assert issubclass(extended, base)
    names = [field.name for field in fields(extended)]
    assert names[: len(fields(base))] == [field.name for field in fields(base)]


@pytest.mark.parametrize("peak_type", [NarrowPeak, BroadPeak])
def test_peaks_check_their_score(peak_type: type[Any]) -> None:
    """Test that the ENCODE peak formats check their score like BED6 records do."""
    values: dict[str, Any] = {"signal_value": 1.0, "p_value": -1.0, "q_value": -1.0}
    if peak_type is NarrowPeak:
        values["peak"] = -1
    with pytest.raises(ValueError, match="score must be between 0 and 1000!"):
        peak_type(refname="chr1", start=0, end=5, name=None, score=1001, strand=None, **values)
