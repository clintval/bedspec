import dataclasses

import pytest

import bedspec
from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedColor
from bedspec import BedGraph
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedType
from bedspec import GenomicSpan
from bedspec import PairBed
from bedspec import PointBed
from bedspec import SimpleBed
from bedspec import Stranded


def test_bed_strand() -> None:
    """Test that BED strands behave as string."""
    assert BedStrand("+") == BedStrand.Positive
    assert BedStrand("-") == BedStrand.Negative
    assert str(BedStrand.Positive) == "+"
    assert str(BedStrand.Negative) == "-"


def test_bed_strand_opposite() -> None:
    """Test that we return an opposite BED strand."""
    assert BedStrand.Positive.opposite() == BedStrand.Negative
    assert BedStrand.Negative.opposite() == BedStrand.Positive


def test_bed_color() -> None:
    """Test the small helper class for BED color."""
    assert str(BedColor(2, 3, 4)) == "2,3,4"


@pytest.mark.parametrize("bed_type", (PointBed, SimpleBed, PairBed))
def test_bed_type_class_hierarchy(bed_type: type[BedType]) -> None:
    """Test that all abstract base classes are subclasses of BedType."""
    assert issubclass(bed_type, BedType)


@pytest.mark.parametrize("bed_type", (Bed2, Bed3, Bed4, Bed5, Bed6, BedGraph, BedPE))
def test_all_bed_types_are_dataclasses(bed_type: type[BedType]) -> None:
    """Test that a simple BED record behaves as expected."""
    assert dataclasses.is_dataclass(bed_type)


def test_locatable_structural_type() -> None:
    """Test that the GenomicSpan structural type is set correctly."""
    span: GenomicSpan = Bed6(
        contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive
    )
    assert isinstance(span, GenomicSpan)


def test_stranded_structural_type() -> None:
    """Test that the Stranded structural type is set correctly."""
    stranded: Stranded = Bed6(
        contig="chr1", start=1, end=2, name="foo", score=3, strand=BedStrand.Positive
    )
    assert isinstance(stranded, Stranded)


def test_dataclass_protocol_structural_type() -> None:
    """Test that the dataclass structural type is set correctly."""
    bed: bedspec._types.DataclassInstance = Bed2(contig="chr1", start=1)
    assert isinstance(bed, bedspec._types.DataclassInstance)
