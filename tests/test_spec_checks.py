from typing import Any

import pytest

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedPE
from bedspec import BedStrand
from bedspec.overlap import OverlapDetector


def bed12(**fields: Any) -> Bed12:
    """Build a valid BED12 record spanning 2 to 10, with any fields replaced."""
    defaults: dict[str, Any] = {
        "refname": "chr1",
        "start": 2,
        "end": 10,
        "name": "bed12",
        "score": 2,
        "strand": BedStrand.Positive,
        "thick_start": 2,
        "thick_end": 10,
        "item_rgb": None,
        "block_count": 2,
        "block_sizes": (2, 3),
        "block_starts": (0, 5),
    }
    return Bed12(**(defaults | fields))


def test_zero_length_features_are_allowed() -> None:
    """Test that a feature may start where it ends, as insertions do."""
    assert len(Bed3(refname="chr1", start=5, end=5)) == 0
    assert len(Bed3(refname="chr1", start=0, end=0)) == 0


def test_zero_length_pairs_are_allowed() -> None:
    """Test that either interval of a pair may start where it ends."""
    record = BedPE(
        refname1="chr1",
        start1=5,
        end1=5,
        refname2="chr2",
        start2=0,
        end2=0,
        name=None,
        score=None,
        strand1=None,
        strand2=None,
    )
    assert record.start1 == record.end1


def test_start_must_not_be_negative() -> None:
    """Test that every kind of record refuses a negative start."""
    with pytest.raises(ValueError, match="start must be greater than or equal to 0!"):
        _ = Bed2(refname="chr1", start=-1)
    with pytest.raises(ValueError, match="start must be greater than or equal to 0!"):
        _ = Bed3(refname="chr1", start=-1, end=5)


def test_end_must_not_come_before_start() -> None:
    """Test that a record refuses an end before its start."""
    with pytest.raises(ValueError, match="end must be greater than or equal to start!"):
        _ = Bed3(refname="chr1", start=5, end=4)


@pytest.mark.parametrize(
    "start1,end1,start2,end2,message",
    [
        (-1, 5, 1, 2, "start1 must be greater than or equal to 0!"),
        (5, 4, 1, 2, "end1 must be greater than or equal to start1!"),
        (1, 2, -1, 5, "start2 must be greater than or equal to 0!"),
        (1, 2, 5, 4, "end2 must be greater than or equal to start2!"),
    ],
)
def test_pairs_check_both_intervals(
    start1: int, end1: int, start2: int, end2: int, message: str
) -> None:
    """Test that a pair checks the start and end of both its intervals."""
    with pytest.raises(ValueError, match=message):
        _ = BedPE(
            refname1="chr1",
            start1=start1,
            end1=end1,
            refname2="chr1",
            start2=start2,
            end2=end2,
            name=None,
            score=None,
            strand1=None,
            strand2=None,
        )


def test_refname_must_not_be_empty() -> None:
    """Test that a record refuses an empty reference sequence name."""
    with pytest.raises(ValueError, match="refname must not be empty!"):
        _ = Bed3(refname="", start=1, end=2)
    with pytest.raises(ValueError, match="refname must not be empty!"):
        _ = Bed2(refname="", start=1)


@pytest.mark.parametrize("name", ["", "x" * 256])
def test_name_must_be_1_to_255_characters(name: str) -> None:
    """Test that a record refuses a name that is empty or longer than 255 characters."""
    with pytest.raises(ValueError, match="name must be 1 to 255 characters long!"):
        _ = Bed4(refname="chr1", start=1, end=2, name=name)
    with pytest.raises(ValueError, match="name must be 1 to 255 characters long!"):
        _ = Bed6(refname="chr1", start=1, end=2, name=name, score=None, strand=None)
    with pytest.raises(ValueError, match="name must be 1 to 255 characters long!"):
        _ = bed12(name=name)


def test_names_at_the_limits_are_allowed() -> None:
    """Test that a name of 1 or 255 characters, or no name, is allowed."""
    assert Bed4(refname="chr1", start=1, end=2, name="x").name == "x"
    assert Bed4(refname="chr1", start=1, end=2, name="x" * 255).name == "x" * 255
    assert Bed4(refname="chr1", start=1, end=2, name=None).name is None


@pytest.mark.parametrize("score", [-1, 1001])
def test_score_must_be_between_0_and_1000(score: int) -> None:
    """Test that a record refuses a score outside 0 to 1000."""
    with pytest.raises(ValueError, match="score must be between 0 and 1000!"):
        _ = Bed5(refname="chr1", start=1, end=2, name=None, score=score)
    with pytest.raises(ValueError, match="score must be between 0 and 1000!"):
        _ = Bed6(refname="chr1", start=1, end=2, name=None, score=score, strand=None)
    with pytest.raises(ValueError, match="score must be between 0 and 1000!"):
        _ = bed12(score=score)


@pytest.mark.parametrize("score", [0, 1000, None])
def test_scores_at_the_limits_are_allowed(score: int | None) -> None:
    """Test that a score of 0 or 1000, or no score, is allowed."""
    assert Bed5(refname="chr1", start=1, end=2, name=None, score=score).score == score


@pytest.mark.parametrize(
    "thick_start,thick_end",
    [(1, 10), (2, 11), (6, 5)],
)
def test_bed12_thick_bounds_must_sit_within_the_feature(thick_start: int, thick_end: int) -> None:
    """Test that the thick part of a BED12 record must sit within the record, in order."""
    with pytest.raises(
        ValueError,
        match="thick_start and thick_end must satisfy start <= thick_start <= thick_end <= end!",
    ):
        _ = bed12(thick_start=thick_start, thick_end=thick_end)


def test_bed12_may_have_no_thick_part() -> None:
    """Test that a BED12 record may have a thick part of zero length."""
    assert bed12(thick_start=2, thick_end=2).thick_start == 2


@pytest.mark.parametrize(
    "block_sizes,block_starts",
    [
        ((4, 3), (0, 3)),
        ((3, 3, 3), (0, 5, 2)),
        ((9, 3), (0, 5)),
    ],
)
def test_bed12_blocks_must_ascend_without_overlapping(
    block_sizes: tuple[int, ...], block_starts: tuple[int, ...]
) -> None:
    """Test that the blocks of a BED12 record must ascend and must not overlap."""
    with pytest.raises(ValueError, match="Blocks must be in ascending order and must not overlap!"):
        _ = bed12(block_count=len(block_sizes), block_sizes=block_sizes, block_starts=block_starts)


def test_bed12_blocks_may_touch() -> None:
    """Test that the blocks of a BED12 record may touch."""
    assert bed12(block_sizes=(5, 3), block_starts=(0, 5)).block_count == 2


def test_overlap_detector_finds_zero_length_features_by_their_flanking_bases() -> None:
    """Test that a zero-length feature overlaps the features on both sides of it."""
    left = Bed3(refname="chr1", start=0, end=5)
    right = Bed3(refname="chr1", start=5, end=9)
    away = Bed3(refname="chr1", start=6, end=9)
    insertion = Bed3(refname="chr1", start=5, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([left, right, away])
    assert set(detector.overlapping(insertion)) == {left, right}
    assert set(detector.enclosing(insertion)) == {left, right}


def test_overlap_detector_holds_zero_length_features() -> None:
    """Test that zero-length features can be added and found."""
    insertion = Bed3(refname="chr1", start=5, end=5)
    at_start = Bed3(refname="chr1", start=0, end=0)
    detector: OverlapDetector[Bed3] = OverlapDetector([insertion, at_start])
    assert set(detector.overlapping(Bed3(refname="chr1", start=4, end=5))) == {insertion}
    assert set(detector.overlapping(Bed3(refname="chr1", start=5, end=6))) == {insertion}
    assert set(detector.overlapping(Bed3(refname="chr1", start=6, end=9))) == set()
    assert set(detector.overlapping(Bed3(refname="chr1", start=0, end=1))) == {at_start}
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=5, end=9))) == {insertion}


def test_bedpe_scores_are_not_limited_to_the_bed_range() -> None:
    """Test that a BEDPE score outside 0 to 1000 is allowed, and its intervals can be taken."""
    record = BedPE(
        refname1="chr1",
        start1=1,
        end1=2,
        refname2="chr2",
        start2=3,
        end2=4,
        name="pair",
        score=1001,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    assert record.score == 1001
    assert record.bed1.score is None
    assert record.bed2.score is None
    assert list(record.territory()) == [record.bed1, record.bed2]
