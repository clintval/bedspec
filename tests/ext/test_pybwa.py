"""Tests of the function that reads pybwa's records as spans."""

from pybwa import AuxHit
from pybwa import to_xa_hits

from bedspec import Bed3
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import register_for_bedspec
from bedspec.ext.pybwa import aux_hit_span
from bedspec.overlap import TreeDetector

XA = "chr4,-97592047,24M,3;chr8,+32368749,24M,3;"


def test_aux_hit_span_converts_a_1_based_hit_on_its_strand() -> None:
    """Test that XA hits, 1-based and closed, are read 0-based and half-open on their strand."""
    minus, plus = to_xa_hits(XA)
    assert aux_hit_span(minus) == [
        Bed6("chr4", start=97592046, end=97592070, name=None, score=None, strand=BedStrand.Negative)
    ]
    assert aux_hit_span(plus) == [
        Bed6("chr8", start=32368748, end=32368772, name=None, score=None, strand=BedStrand.Positive)
    ]


def test_a_registered_hit_is_found_by_its_first_base() -> None:
    """Test that a registered hit is read through its function, not its 1-based attributes."""
    minus, _ = to_xa_hits(XA)
    first_base = TreeDetector([Bed3("chr4", start=97592046, end=97592047)])
    assert not first_base.overlaps(minus)
    with register_for_bedspec(AuxHit, aux_hit_span):
        assert first_base.overlaps(minus)
