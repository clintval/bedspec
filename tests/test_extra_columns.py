from dataclasses import fields
from pathlib import Path
from typing import Any

import pytest
from typeline import ExtraColumns

from bedspec import Bed3
from bedspec import Bed3N
from bedspec import Bed4
from bedspec import Bed4N
from bedspec import Bed5
from bedspec import Bed5N
from bedspec import Bed6
from bedspec import Bed6N
from bedspec import Bed12
from bedspec import Bed12N
from bedspec import BedReader
from bedspec import BedStrand
from bedspec import BedWriter


@pytest.mark.parametrize(
    "extended,base",
    [(Bed3N, Bed3), (Bed4N, Bed4), (Bed5N, Bed5), (Bed6N, Bed6), (Bed12N, Bed12)],
)
def test_extended_bed_types_add_extra_columns_to_their_base(
    extended: type[Any], base: type[Any]
) -> None:
    """Test that each extended BED type is its base type with a last field for extra columns."""
    assert issubclass(extended, base)
    assert [field.name for field in fields(extended)] == [
        *(field.name for field in fields(base)),
        "extra",
    ]
    assert fields(extended)[-1].type == ExtraColumns


def test_bed_reader_keeps_extra_columns(tmp_path: Path) -> None:
    """Test that the BED reader keeps any columns past a BED6 record as text."""
    (tmp_path / "test.bed").write_text(
        "track name=peaks\nchr1\t1\t2\tpeak1\t0\t+\n\nchr1\t5\t9\tpeak2\t7\t-\t3.2\t.\n"
    )

    with BedReader.from_path[Bed6N](tmp_path / "test.bed") as reader:
        assert list(reader) == [
            Bed6N("chr1", start=1, end=2, name="peak1", score=0, strand=BedStrand.Positive),
            Bed6N(
                "chr1",
                start=5,
                end=9,
                name="peak2",
                score=7,
                strand=BedStrand.Negative,
                extra=("3.2", "."),
            ),
        ]


def test_extra_columns_round_trip_through_a_bed_file(tmp_path: Path) -> None:
    """Test that extra columns are written back after the BED columns, unchanged."""
    text = "chr1\t1\t2\tpeak1\t0\t+\nchr1\t5\t9\tpeak2\t7\t-\t3.2\t.\n"
    (tmp_path / "in.bed").write_text(text)

    with (
        BedReader.from_path[Bed6N](tmp_path / "in.bed") as reader,
        BedWriter.from_path[Bed6N](tmp_path / "out.bed") as writer,
    ):
        for record in reader:
            writer.write(record)

    assert (tmp_path / "out.bed").read_text() == text
