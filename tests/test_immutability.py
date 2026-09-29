from dataclasses import FrozenInstanceError
from dataclasses import fields
from dataclasses import replace
from io import StringIO

import pytest

from bedspec import Bed3
from bedspec import Bed12
from bedspec import Bed12N
from bedspec import BedColor
from bedspec import BedReader
from bedspec import BedStrand

BED12 = Bed12(
    refname="chr1",
    start=2,
    end=10,
    name="bed12",
    score=2,
    strand=BedStrand.Positive,
    thick_start=2,
    thick_end=10,
    item_rgb=BedColor(1, 2, 3),
    block_count=2,
    block_sizes=(2, 3),
    block_starts=(0, 5),
)


def test_records_cannot_be_changed() -> None:
    """Test that a record's fields cannot be changed after it is built."""
    record = Bed3(refname="chr1", start=1, end=2)
    with pytest.raises(FrozenInstanceError):
        record.start = 5  # type: ignore[misc]


def test_colors_cannot_be_changed() -> None:
    """Test that a color's values cannot be changed after it is built."""
    color = BedColor(1, 2, 3)
    with pytest.raises(FrozenInstanceError):
        color.r = 5  # type: ignore[misc]


def test_replace_builds_a_checked_copy() -> None:
    """Test that a changed copy of a record is checked like any other record."""
    record = Bed3(refname="chr1", start=1, end=2)
    assert replace(record, end=9) == Bed3(refname="chr1", start=1, end=9)
    with pytest.raises(ValueError, match="end must be greater than or equal to start!"):
        replace(record, end=0)


def test_bed12_records_can_be_hashed() -> None:
    """Test that BED12 records, blocks and all, can be hashed and kept in sets."""
    assert len({BED12, replace(BED12)}) == 1
    values = {field.name: getattr(BED12, field.name) for field in fields(BED12)}
    assert len({Bed12N(**values, extra=("x",)), Bed12N(**values, extra=("x",))}) == 1


def test_bed12_blocks_are_read_as_tuples() -> None:
    """Test that the blocks of a BED12 record are read as tuples."""
    line = "chr1\t2\t10\tbed12\t2\t+\t2\t10\t1,2,3\t2\t2,3\t0,5\n"
    with BedReader[Bed12](StringIO(line)) as reader:
        (record,) = list(reader)
    assert record == BED12
    assert hash(record) == hash(BED12)
