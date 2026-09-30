import pytest

from bedspec import Bed3
from bedspec import Bed4
from bedspec.overlap import TreeDetector


def test_overlap_detector_as_iterable() -> None:
    """Test we can iterate over all the intervals we put into the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed3(refname="chr2", start=4, end=5)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2])
    assert set(detector) == {bed1, bed2}


def test_we_can_mix_types_in_the_overlap_detector() -> None:
    """Test mix input types when building the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: TreeDetector[Bed3 | Bed4] = TreeDetector([bed1, bed2])
    assert set(detector) == {bed1, bed2}


def test_we_can_add_a_feature_to_the_overlap_detector() -> None:
    """Test we can add a feature to the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: TreeDetector[Bed3 | Bed4] = TreeDetector()
    detector.add(bed1)
    detector.add(bed2)
    assert set(detector) == {bed1, bed2}


def test_we_can_add_all_features_to_the_overlap_detector() -> None:
    """Test we can add all features to the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr2", start=4, end=5, name="Clint Valentine")
    detector: TreeDetector[Bed3 | Bed4] = TreeDetector()
    beds: list[Bed3 | Bed4] = [bed1, bed2]
    detector.add(*beds)
    assert set(detector) == {bed1, bed2}


def test_we_can_query_with_different_type_in_the_overlap_detector() -> None:
    """Test we can query with a different type in the overlap detector."""
    bed1 = Bed3(refname="chr1", start=1, end=2)
    bed2 = Bed4(refname="chr1", start=1, end=2, name="Clint Valentine")
    detector: TreeDetector[Bed3] = TreeDetector([bed1])
    assert set(detector.overlapping(bed2)) == {bed1}


def test_we_can_those_enclosing_intervals() -> None:
    """Test that we can get intervals enclosing a given query feature."""
    bed1 = Bed3(refname="chr1", start=1, end=5)
    bed2 = Bed3(refname="chr1", start=3, end=9)
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2])
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
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2])
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
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2, bed3])

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
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2, bed3])

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
    detector: TreeDetector[Bed3] = TreeDetector([bed1, bed2, bed3])

    assert set(detector) == {bed1, bed2, bed3}

    assert not detector.overlaps(Bed3("chr1", start=5, end=6))

    detector.add(Bed3("chr1", start=5, end=6))

    assert detector.overlaps(Bed3("chr1", start=5, end=6))


def test_half_open_features_that_abut_do_not_overlap() -> None:
    """Test that half-open features which share only an endpoint do not overlap."""
    bed = Bed3(refname="chr1", start=10, end=20)
    detector: TreeDetector[Bed3] = TreeDetector([bed])
    assert not detector.overlaps(Bed3(refname="chr1", start=5, end=10))
    assert not detector.overlaps(Bed3(refname="chr1", start=20, end=25))
    assert list(detector.overlapping(Bed3(refname="chr1", start=9, end=11))) == [bed]
    assert list(detector.overlapping(Bed3(refname="chr1", start=19, end=21))) == [bed]


def test_querying_an_unknown_reference_finds_nothing() -> None:
    """Test that querying a reference sequence with no features returns no overlaps."""
    detector: TreeDetector[Bed3] = TreeDetector([Bed3(refname="chr1", start=10, end=20)])
    assert list(detector.overlapping(Bed3(refname="chr2", start=10, end=20))) == []


def features_added_out_of_order() -> list[Bed3]:
    """Return features on two references, added out of order by start."""
    return [
        Bed3(refname="chr1", start=30, end=40),
        Bed3(refname="chr2", start=5, end=9),
        Bed3(refname="chr1", start=20, end=25),
        Bed3(refname="chr1", start=10, end=15),
        Bed3(refname="chr2", start=1, end=3),
    ]


def test_querying_while_iterating_visits_every_feature_once() -> None:
    """Test that a query inside a loop over the detector neither skips nor repeats features."""
    features = features_added_out_of_order()
    detector: TreeDetector[Bed3] = TreeDetector(features)

    seen: list[Bed3] = []
    for feature in detector:
        seen.append(feature)
        _ = list(detector.overlapping(feature))

    assert sorted(seen, key=lambda f: (f.refname, f.start)) == sorted(
        features, key=lambda f: (f.refname, f.start)
    )


def test_iteration_order_does_not_depend_on_queries() -> None:
    """Test that features are iterated by reference, then by start, whether or not a query ran."""
    queried: TreeDetector[Bed3] = TreeDetector(features_added_out_of_order())
    _ = queried.overlaps(Bed3(refname="chr1", start=0, end=1))
    fresh: TreeDetector[Bed3] = TreeDetector(features_added_out_of_order())

    expected = [(f.refname, f.start) for f in features_added_out_of_order()]
    expected.sort(key=lambda key: (key[0] != "chr1", key[1]))
    assert [(f.refname, f.start) for f in fresh] == expected
    assert [(f.refname, f.start) for f in queried] == expected


@pytest.mark.parametrize("refname", ["chr1", "chr3"])
def test_adding_while_iterating_raises(refname: str) -> None:
    """Test that adding a feature while iterating raises rather than skipping or repeating any."""
    detector: TreeDetector[Bed3] = TreeDetector(features_added_out_of_order())

    with pytest.raises(RuntimeError, match="changed during iteration"):
        for _ in detector:
            detector.add(Bed3(refname=refname, start=0, end=1))  # noqa: B909
