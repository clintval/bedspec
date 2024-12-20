from bedspec import Bed3
from bedspec import Bed4
from bedspec.overlap import OverlapDetector


def test_overlap_detector_as_iterable() -> None:
    """Test we can iterate over all the intervals we put into the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed3(refname="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2])
    assert set(detector) == {bed1, bed2}


def test_we_can_mix_types_in_the_overlap_detector() -> None:
    """Test mix input types when building the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: OverlapDetector[Bed3 | Bed4] = OverlapDetector([bed1, bed2])
    assert set(detector) == {bed1, bed2}


def test_we_can_add_a_feature_to_the_overlap_detector() -> None:
    """Test we can add a feature to the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: OverlapDetector[Bed3 | Bed4] = OverlapDetector()
    detector.add(bed1)
    detector.add(bed2)
    assert set(detector) == {bed1, bed2}


def test_we_can_add_all_features_to_the_overlap_detector() -> None:
    """Test we can add all features to the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: OverlapDetector[Bed3 | Bed4] = OverlapDetector()
    beds: list[Bed3 | Bed4] = [bed1, bed2]
    detector.add(*beds)
    assert set(detector) == {bed1, bed2}


def test_we_can_query_with_different_type_in_the_overlap_detector() -> None:
    """Test we can query with a different type in the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr1", start=1, end=2, name="Clint Valentine")
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1])
    assert set(detector.overlapping(bed2)) == {bed1}


def test_we_can_those_enclosing_intervals() -> None:
    """Test that we can get intervals enclosing a given query feature."""
    bed1 = Bed3(refname="chr1", start=1, end=5)
    bed2 = Bed3(refname="chr1", start=3, end=9)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2])
    assert set(detector.enclosing(Bed3(refname="chr1", start=2, end=5))) == {bed1}
    assert set(detector.enclosing(Bed3(refname="chr1", start=3, end=8))) == {bed2}
    assert set(detector.enclosing(Bed3(refname="chr1", start=4, end=9))) == {bed2}
    assert set(detector.enclosing(Bed3(refname="chr1", start=3, end=9))) == {bed2}
    assert set(detector.enclosing(Bed3(refname="chr1", start=2, end=10))) == set()
    assert set(detector.enclosing(Bed3(refname="chr1", start=1, end=10))) == set()


def test_we_can_those_enclosed_by_intervals() -> None:
    """Test that we can get intervals enclosed by a given query feature."""
    bed1 = Bed3(refname="chr1", start=1, end=5)
    bed2 = Bed3(refname="chr1", start=3, end=9)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2])
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=2, end=5))) == set()
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=3, end=8))) == set()
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=4, end=9))) == set()
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=3, end=9))) == {bed2}
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=2, end=10))) == {bed2}
    assert set(detector.enclosed_by(Bed3(refname="chr1", start=1, end=10))) == {bed1, bed2}


def test_we_can_query_for_overlapping_features() -> None:
    """Test we can query for features that overlap using the overlap detector."""
    bed1 = Bed3(refname="chr1", start=2, end=5)
    bed2 = Bed3(refname="chr1", start=4, end=10)
    bed3 = Bed3(refname="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert set(detector.overlapping(Bed3("chr1", start=0, end=1))) == set()
    assert set(detector.overlapping(Bed3("chr1", start=2, end=3))) == {bed1}
    assert set(detector.overlapping(Bed3("chr1", start=4, end=5))) == {bed1, bed2}
    assert set(detector.overlapping(Bed3("chr1", start=5, end=6))) == {bed2}
    assert set(detector.overlapping(Bed3("chr2", start=0, end=1))) == set()
    assert set(detector.overlapping(Bed3("chr2", start=4, end=5))) == {bed3}


def test_we_can_query_if_at_least_one_feature_overlaps() -> None:
    """Test we can query if at least one feature overlaps using the overlap detector."""
    bed1 = Bed3(refname="chr1", start=2, end=5)
    bed2 = Bed3(refname="chr1", start=4, end=10)
    bed3 = Bed3(refname="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert not detector.overlaps(Bed3("chr1", start=0, end=1))
    assert detector.overlaps(Bed3("chr1", start=2, end=3))
    assert detector.overlaps(Bed3("chr1", start=4, end=5))
    assert detector.overlaps(Bed3("chr1", start=5, end=6))
    assert not detector.overlaps(Bed3("chr2", start=0, end=1))
    assert detector.overlaps(Bed3("chr2", start=4, end=5))


def test_we_can_add_to_the_overlap_detector_after_and_before_queries() -> None:
    """Test we can add to the overlap detector after and before queries."""
    bed1 = Bed3(refname="chr1", start=2, end=5)
    bed2 = Bed3(refname="chr1", start=6, end=10)
    bed3 = Bed3(refname="chr2", start=4, end=5)
    detector: OverlapDetector[Bed3] = OverlapDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert not detector.overlaps(Bed3("chr1", start=5, end=6))

    detector.add(Bed3("chr1", start=5, end=6))

    assert detector.overlaps(Bed3("chr1", start=5, end=6))
