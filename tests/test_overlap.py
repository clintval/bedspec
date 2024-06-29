from dataclasses import dataclass
from typing import TypeAlias

import pytest

from bedspec import Bed3
from bedspec import Bed4
from bedspec.overlap import OverlapDetector


def test_overlap_detector_as_iterable() -> None:
    """Test we can iterate over all the intervals we put into the overlap detector."""
    bed1 = Bed3(contig="chr1", start=1, end=2)
    bed2 = Bed3(contig="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2])
    assert list(detector) == [bed1, bed2]


def test_we_can_mix_types_in_the_overlap_detector() -> None:
    """Test mix input types when building the overlap detector."""
    bed1 = Bed3(contig="chr1", start=1, end=2)
    bed2 = Bed4(contig="chr2", start=4, end=5, name="Clint Valentine")
    detector: OverlapDetector[Bed3 | Bed4] = OverlapDetector([bed1, bed2])
    assert list(detector) == [bed1, bed2]


def test_we_can_add_a_feature_to_the_overlap_detector() -> None:
    """Test we can add a feature to the overlap detector."""
    bed1 = Bed3(contig="chr1", start=1, end=2)
    bed2 = Bed4(contig="chr2", start=4, end=5, name="Clint Valentine")
    detector: OverlapDetector[Bed3 | Bed4] = OverlapDetector()
    detector.add(bed1)
    detector.add(bed2)
    assert list(detector) == [bed1, bed2]


def test_that_we_require_hashable_features_in_the_overlap_detector() -> None:
    """Test that we require hashable features in the overlap detector."""

    @dataclass
    class MissingHashFeature:
        contig: str
        start: int
        end: int

    feature: MissingHashFeature = MissingHashFeature("chr1", 2, 3)
    detector: OverlapDetector[MissingHashFeature] = OverlapDetector()

    with pytest.raises(ValueError):
        detector.add(feature)


def test_structural_type_reference_name_raises_if_not_found() -> None:

    @dataclass(eq=True, frozen=True)
    class BadInterval:
        chromosome_name: str
        start: int
        end: int

    feature: BadInterval = BadInterval("chr1", 1, 2)

    with pytest.raises(ValueError):
        OverlapDetector._reference_sequence_name(feature)  # type: ignore[type-var]

def test_we_can_add_all_features_to_the_overlap_detector() -> None:
    """Test we can add all features to the overlap detector."""
    bed1 = Bed3(contig="chr1", start=1, end=2)
    bed2 = Bed4(contig="chr2", start=4, end=5, name="Clint Valentine")
    detector: OverlapDetector[Bed3 | Bed4] = OverlapDetector()
    detector.add_all([bed1, bed2])
    assert list(detector) == [bed1, bed2]


def test_we_can_query_with_different_type_in_the_overlap_detector() -> None:
    """Test we can query with a different type in the overlap detector."""
    bed1 = Bed3(contig="chr1", start=1, end=2)
    bed2 = Bed4(contig="chr1", start=1, end=2, name="Clint Valentine")
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1])
    assert list(detector.get_overlapping(bed2)) == [bed1]


def test_we_can_query_for_overlapping_features() -> None:
    """Test we can query for features that overlap using the overlap detector."""
    bed1 = Bed3(contig="chr1", start=2, end=5)
    bed2 = Bed3(contig="chr1", start=4, end=10)
    bed3 = Bed3(contig="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2, bed3])

    assert list(detector) == [bed1, bed2, bed3]

    assert list(detector.get_overlapping(Bed3("chr1", 0, 1))) == []
    assert list(detector.get_overlapping(Bed3("chr1", 2, 3))) == [bed1]
    assert list(detector.get_overlapping(Bed3("chr1", 4, 5))) == [bed1, bed2]
    assert list(detector.get_overlapping(Bed3("chr1", 5, 6))) == [bed2]
    assert list(detector.get_overlapping(Bed3("chr2", 0, 1))) == []
    assert list(detector.get_overlapping(Bed3("chr2", 4, 5))) == [bed3]


def test_we_can_query_if_at_least_one_feature_overlaps() -> None:
    """Test we can query if at least one feature overlaps using the overlap detector."""
    bed1 = Bed3(contig="chr1", start=2, end=5)
    bed2 = Bed3(contig="chr1", start=4, end=10)
    bed3 = Bed3(contig="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2, bed3])

    assert list(detector) == [bed1, bed2, bed3]

    assert not detector.overlaps_any(Bed3("chr1", 0, 1))
    assert detector.overlaps_any(Bed3("chr1", 2, 3))
    assert detector.overlaps_any(Bed3("chr1", 4, 5))
    assert detector.overlaps_any(Bed3("chr1", 5, 6))
    assert not detector.overlaps_any(Bed3("chr2", 0, 1))
    assert detector.overlaps_any(Bed3("chr2", 4, 5))


def test_we_support_features_with_all_three_common_reference_sequence_name_properties() -> None:
    """Test that we can store features with either of 3 reference sequence name properties."""

    @dataclass(eq=True, frozen=True)
    class FeatureWithChrom:
        chrom: str
        start: int
        end: int

    @dataclass(eq=True, frozen=True)
    class FeatureWithContig:
        contig: str
        start: int
        end: int

    @dataclass(eq=True, frozen=True)
    class FeatureWithRefname:
        refname: str
        start: int
        end: int

    feature_with_chrom: FeatureWithChrom = FeatureWithChrom("chr1", 1, 3)
    feature_with_contig: FeatureWithContig = FeatureWithContig("chr1", 1, 3)
    feature_with_refname: FeatureWithRefname = FeatureWithRefname("chr1", 1, 3)

    AllKinds: TypeAlias = FeatureWithChrom | FeatureWithContig | FeatureWithRefname
    features: list[AllKinds] = [feature_with_chrom, feature_with_contig, feature_with_refname]
    detector: OverlapDetector[AllKinds] = OverlapDetector(features)

    assert list(detector) == features
