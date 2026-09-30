from io import StringIO
from pathlib import Path

import pytest

from bedspec import Bed4
from bedspec import BedReader
from bedspec import BedWriter


def test_quotes_in_a_name_are_ordinary_text(tmp_path: Path) -> None:
    """Test that BED has no quoting, so names holding quotes are written and read as they are."""
    records = [
        Bed4(refname="chr1", start=1, end=2, name='"peak'),
        Bed4(refname="chr1", start=3, end=4, name='say "hi"'),
        Bed4(refname="chr1", start=5, end=6, name="plain"),
    ]
    with BedWriter.from_path[Bed4](tmp_path / "test.bed") as writer:
        for record in records:
            writer.write(record)

    assert (tmp_path / "test.bed").read_text() == (
        'chr1\t1\t2\t"peak\nchr1\t3\t4\tsay "hi"\nchr1\t5\t6\tplain\n'
    )
    with BedReader.from_path[Bed4](tmp_path / "test.bed") as reader:
        assert list(reader) == records


def test_a_name_holding_a_tab_is_refused() -> None:
    """Test that a name holding a tab, which BED cannot hold, is refused rather than written."""
    with pytest.raises(ValueError, match=r"without quoting"):
        BedWriter[Bed4](StringIO()).write(Bed4(refname="chr1", start=1, end=2, name="a\tb"))
