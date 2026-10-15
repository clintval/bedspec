"""Tests of registering how detectors and territories find the spans of other classes of records."""

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Protocol

import pytest
from pybgzf import IndexFormat
from typing_extensions import override

from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedPE
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import Locatable
from bedspec import ReferenceSpan
from bedspec import Territory
from bedspec import register_locatable
from bedspec.overlap import TabixDetector
from bedspec.overlap import TreeDetector


@dataclass(frozen=True)
class Read:
    """A read aligned in blocks, placed on a reference the way no BED record is."""

    contig: str | None
    blocks: tuple[tuple[int, int], ...]
    reverse: bool = False


@dataclass(frozen=True)
class Spliced(Read):
    """A kind of read, which uses the registration of its parent class."""


@dataclass(frozen=True)
class OneBased:
    """A record with a refname, start and end in 1-based closed coordinates, like pybwa's AuxHit."""

    refname: str
    start: int
    end: int


@dataclass(frozen=True)
class Plain:
    """A record that is not a BED record but has a refname, start and end."""

    refname: str
    start: int
    end: int


@dataclass(frozen=True)
class Located:
    """A record that is not a BED record but says what its spans are."""

    parts: tuple[Bed3, ...]

    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield the parts of this record."""
        yield from self.parts


@dataclass(frozen=True)
class Hollow(Bed3):
    """A BED record that has no spans."""

    @override
    def spans(self) -> Iterator[ReferenceSpan]:
        """Yield no spans."""
        yield from ()


class Placed(Protocol):
    """A protocol, which cannot be registered."""

    contig: str


def read_blocks(read: Read) -> list[Bed6]:
    """Return the blocks of a read on its strand, or none if it is placed on no reference."""
    if read.contig is None:
        return []
    strand = BedStrand.Negative if read.reverse else BedStrand.Positive
    return [
        Bed6(read.contig, start=start, end=end, name=None, score=None, strand=strand)
        for start, end in read.blocks
    ]


def read_extent(read: Read) -> list[Bed3]:
    """Return the whole extent of a read, from its first block to its last."""
    assert read.contig is not None
    return [Bed3(read.contig, start=read.blocks[0][0], end=read.blocks[-1][1])]


def one_based(record: OneBased) -> list[Bed3]:
    """Return the span of a 1-based closed record, converted to 0-based half-open."""
    return [Bed3(record.refname, start=record.start - 1, end=record.end)]


SPLICED = Read("chr1", ((100, 110), (200, 210)))
"""A read whose two blocks flank an intron from 110 to 200."""

EXON1, INTRON, EXON2 = (
    Bed3("chr1", start=105, end=106),
    Bed3("chr1", start=150, end=160),
    Bed3("chr1", start=205, end=206),
)


def test_a_registered_record_is_found_by_its_spans_and_returned_as_itself() -> None:
    """Test that a registered record is indexed by the spans its function returns."""
    elsewhere = Read("chr2", ((0, 10),))
    with register_locatable(Read, read_blocks):
        detector = TreeDetector([SPLICED, elsewhere])
        assert list(detector) == [SPLICED, elsewhere]
        assert list(detector.overlapping(EXON1)) == [SPLICED]
        assert list(detector.overlapping(EXON2)) == [SPLICED]
        assert not detector.overlaps(INTRON)
        assert list(detector.enclosing(EXON1)) == [SPLICED]
        assert list(detector.enclosed_by(Bed3("chr1", start=0, end=1000))) == [SPLICED]


def test_a_registered_query_finds_what_any_of_its_spans_overlaps_once() -> None:
    """Test that a query of several spans finds each feature any span overlaps, once."""
    across = Bed3("chr1", start=100, end=210)
    targets = [EXON1, INTRON, EXON2, across]
    hollow = Hollow("chr9", start=0, end=1)
    with register_locatable(Read, read_blocks):
        for detector in (TreeDetector(targets), TreeDetector([*targets, hollow])):
            assert detector.overlaps(SPLICED)
            assert sorted(detector.overlapping(SPLICED), key=str) == sorted(
                [EXON1, EXON2, across], key=str
            )
            assert not detector.overlaps(Read("chr1", ((300, 310),)))


def test_a_feature_added_twice_is_found_twice_by_a_query_of_several_spans() -> None:
    """Test that a query of several spans yields a repeated feature as often as one span does."""
    with register_locatable(Read, read_blocks):
        detector = TreeDetector([EXON1, EXON1])
        assert list(detector.overlapping(SPLICED)) == [EXON1, EXON1]


def test_a_subclass_uses_the_registration_of_its_parent() -> None:
    """Test that a subclass of a registered class is read by its parent's function."""
    spliced = Spliced("chr1", ((100, 110), (200, 210)))
    with register_locatable(Read, read_blocks):
        assert not TreeDetector([spliced]).overlaps(INTRON)
        with register_locatable(Spliced, read_extent):
            assert TreeDetector([spliced]).overlaps(INTRON)
            assert not TreeDetector([SPLICED]).overlaps(INTRON)
        assert not TreeDetector([spliced]).overlaps(INTRON)


def test_a_registration_in_a_with_block_is_undone_when_the_block_exits() -> None:
    """Test that a registration made for a block holds only inside it."""
    with register_locatable(Read, read_blocks):
        assert TreeDetector([SPLICED]).overlaps(EXON1)
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([SPLICED])


def test_a_registration_in_a_with_block_restores_the_one_before_it() -> None:
    """Test that leaving nested blocks restores each earlier registration in turn."""
    registration = register_locatable(Read, read_blocks)
    try:
        with register_locatable(Read, read_extent):
            assert TreeDetector([SPLICED]).overlaps(INTRON)
            with register_locatable(Read, read_blocks):
                assert not TreeDetector([SPLICED]).overlaps(INTRON)
            assert TreeDetector([SPLICED]).overlaps(INTRON)
        assert not TreeDetector([SPLICED]).overlaps(INTRON)
    finally:
        _ = registration.__exit__(None, None, None)
    with pytest.raises(TypeError):
        _ = TreeDetector([SPLICED])


def test_a_registration_in_a_with_block_is_undone_when_the_block_raises() -> None:
    """Test that a registration made for a block is undone even when the block raises."""
    with pytest.raises(RuntimeError), register_locatable(Read, read_blocks):
        raise RuntimeError("the block failed")
    with pytest.raises(TypeError):
        _ = TreeDetector([SPLICED])


def test_a_registration_replaced_before_its_block_holds_for_the_block() -> None:
    """Test that entering a replaced registration installs it for the block, then yields back."""

    @dataclass(frozen=True)
    class Local(Read):
        """A kind of read registered only by this test."""

    spliced = Local("chr1", SPLICED.blocks)
    with register_locatable(Local, read_blocks):
        stale = register_locatable(Local, read_extent)
        _ = register_locatable(Local, read_blocks)
        with stale:
            assert TreeDetector([spliced]).overlaps(INTRON)
        assert not TreeDetector([spliced]).overlaps(INTRON)
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([spliced])


def test_entering_a_registration_keeps_its_own_function_registered_again_before() -> None:
    """Test that entering a registration restores what was registered before the block."""

    @dataclass(frozen=True)
    class Local(Read):
        """A kind of read registered only by this test."""

    spliced = Local("chr1", SPLICED.blocks)
    with register_locatable(Local, read_extent):
        scoped = register_locatable(Local, read_blocks)
        _ = register_locatable(Local, read_extent)
        _ = register_locatable(Local, read_blocks)
        with scoped:
            assert not TreeDetector([spliced]).overlaps(INTRON)
        assert not TreeDetector([spliced]).overlaps(INTRON)
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([spliced])


def test_entering_a_registration_twice_restores_each_level() -> None:
    """Test that entering the same registration twice restores at each block's exit."""

    @dataclass(frozen=True)
    class Local(Read):
        """A kind of read registered only by this test."""

    spliced = Local("chr1", SPLICED.blocks)
    with register_locatable(Local, read_extent):
        registration = register_locatable(Local, read_blocks)
        with registration:
            with registration:
                assert not TreeDetector([spliced]).overlaps(INTRON)
            assert not TreeDetector([spliced]).overlaps(INTRON)
        assert TreeDetector([spliced]).overlaps(INTRON)
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([spliced])


def test_a_permanent_registration_holds_until_replaced() -> None:
    """Test that a registration made without a with statement holds until replaced."""
    first = register_locatable(Read, read_blocks)
    try:
        assert not TreeDetector([SPLICED]).overlaps(INTRON)
        second = register_locatable(Read, read_extent)
        assert TreeDetector([SPLICED]).overlaps(INTRON)
        _ = second.__exit__(None, None, None)
        assert not TreeDetector([SPLICED]).overlaps(INTRON)
    finally:
        _ = first.__exit__(None, None, None)


def test_a_registration_beats_a_refname_start_and_end(tmp_path: Path) -> None:
    """Test that a registered class with native-looking attributes is read through its function.

    Read natively the hit would be [11, 20); its registration reads it as [10, 20), so every
    assertion here fails if its own attributes win.
    """
    hit = OneBased("chr1", start=11, end=20)
    first_base = Bed3("chr1", start=10, end=11)
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(first_base)
    with register_locatable(OneBased, one_based):
        assert TreeDetector([hit]).overlaps(first_base)
        assert TreeDetector([first_base]).overlaps(hit)
        assert list(TreeDetector([first_base]).enclosed_by(hit)) == [first_base]
        assert list(TreeDetector([Bed3("chr1", start=11, end=30)]).enclosing(hit)) == []
        assert list(Territory([hit])) == [Bed3("chr1", start=10, end=20)]
        assert hit not in Territory([Bed3("chr1", start=11, end=20)])
        with TabixDetector[Bed3](path) as detector:
            assert detector.overlaps(hit)
            assert list(detector.enclosed_by(hit)) == [first_base]
    assert not TreeDetector([hit]).overlaps(first_base)


@pytest.mark.parametrize("kind", [SPLICED, object, Placed, "Read"])
def test_only_a_class_of_records_can_be_registered(kind: Any) -> None:
    """Test that registering an instance, object itself, or a protocol raises."""
    with pytest.raises(TypeError, match="Only a class of records can be registered"):
        _ = register_locatable(kind, read_blocks)


def test_a_record_whose_spans_cannot_be_found_raises(tmp_path: Path) -> None:
    """Test that an unregistered record without spans raises, pointing at what to do."""
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(Bed3("chr1", start=0, end=10))
    message = r"register_locatable\(\), or use one from bedspec\.ext"
    with pytest.raises(TypeError, match=message):
        _ = TreeDetector([SPLICED])
    with pytest.raises(TypeError, match=message):
        _ = TreeDetector([EXON1]).overlaps(SPLICED)
    with pytest.raises(TypeError, match=message):
        _ = Territory([SPLICED])
    with TabixDetector[Bed3](path) as detector, pytest.raises(TypeError, match=message):
        _ = detector.overlaps(SPLICED)
    assert SPLICED not in Territory([Bed3("chr1", start=0, end=1000)])


def test_adding_features_adds_none_when_a_registered_function_raises() -> None:
    """Test that adding features adds none of them when one feature's function raises."""
    detector: TreeDetector[Bed3 | Read] = TreeDetector([EXON1])
    with register_locatable(Read, read_extent), pytest.raises(IndexError):
        detector.add(SPLICED, Read("chr1", ()))
    assert list(detector) == [EXON1]


def test_a_record_with_a_spans_method_is_found_by_its_spans() -> None:
    """Test that an unregistered record with a spans() method is found by them."""
    located = Located((Bed3("chr1", start=100, end=110), Bed3("chr1", start=200, end=210)))
    assert list(TreeDetector([located]).overlapping(EXON2)) == [located]
    assert TreeDetector([EXON1, EXON2]).overlaps(located)
    assert not TreeDetector([located]).overlaps(INTRON)
    assert list(Territory([located])) == list(located.parts)


def test_a_feature_with_no_spans_is_held_but_never_found() -> None:
    """Test that a feature with no spans, like an unmapped read, is held and iterated."""
    unmapped, hollow = Read(None, ()), Hollow("chr1", start=0, end=1000)
    everything = Bed3("chr1", start=0, end=1000)
    with register_locatable(Read, read_blocks):
        detector = TreeDetector([unmapped, hollow, SPLICED])
        assert list(detector) == [unmapped, hollow, SPLICED]
        assert list(detector.overlapping(everything)) == [SPLICED]
        assert list(detector.enclosed_by(everything)) == [SPLICED]


def test_a_query_with_no_spans_finds_nothing(tmp_path: Path) -> None:
    """Test that every method of both detectors finds nothing for a query with no spans."""
    unmapped, hollow = Read(None, ()), Hollow("chr1", start=0, end=1000)
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(Bed3("chr1", start=0, end=1000))
    with register_locatable(Read, read_blocks), TabixDetector[Bed3](path) as tabix:
        tree = TreeDetector([Bed3("chr1", start=0, end=1000)])
        for detector in (tree, tabix):
            for query in (unmapped, Located(())):
                assert not detector.overlaps(query)
                assert list(detector.overlapping(query)) == []
                assert list(detector.enclosing(query)) == []
                assert list(detector.enclosed_by(query)) == []
        assert tree.overlaps(hollow)


def test_a_territory_adds_no_bases_for_and_holds_no_feature_with_no_spans() -> None:
    """Test that a feature with no spans adds nothing to a territory and is never in one."""
    unmapped, hollow = Read(None, ()), Hollow("chr1", start=0, end=10)
    whole = Territory([Bed3("chr1", start=0, end=1000)])
    with register_locatable(Read, read_blocks):
        assert list(Territory([unmapped, hollow, SPLICED])) == [
            Bed3("chr1", start=100, end=110),
            Bed3("chr1", start=200, end=210),
        ]
        assert unmapped not in whole
        assert SPLICED in whole
    assert hollow not in whole


def test_a_stranded_query_matches_each_of_its_spans_by_its_own_strand() -> None:
    """Test that a stranded query by a registered record compares each span's strand."""
    plus = Bed6("chr1", start=100, end=110, name=None, score=None, strand=BedStrand.Positive)
    minus = Bed6("chr1", start=200, end=210, name=None, score=None, strand=BedStrand.Negative)
    detector = TreeDetector([plus, minus])
    with register_locatable(Read, read_blocks):
        assert list(detector.overlapping(SPLICED, stranded=True)) == [plus]
        reverse = Read("chr1", SPLICED.blocks, reverse=True)
        assert list(detector.overlapping(reverse, stranded=True)) == [minus]
    with register_locatable(Read, read_extent):
        assert not detector.overlaps(SPLICED, stranded=True)


def test_enclosing_and_enclosed_by_refuse_a_query_of_several_spans(tmp_path: Path) -> None:
    """Test that both detectors refuse to say what encloses, or is enclosed by, several spans."""
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(Bed3("chr1", start=0, end=1000))
    one_block = Read("chr1", ((100, 110),))
    with register_locatable(Read, read_blocks), TabixDetector[Bed3](path) as tabix:
        tree = TreeDetector([Bed3("chr1", start=0, end=1000)])
        for detector in (tree, tabix):
            for method in ("enclosing", "enclosed_by"):
                with pytest.raises(ValueError, match=rf"{method}\(\) needs a query of one span"):
                    _ = list(getattr(detector, method)(SPLICED))
            assert list(detector.enclosing(one_block)) == [Bed3("chr1", start=0, end=1000)]
            assert list(detector.enclosed_by(one_block)) == []


def test_a_tabix_detector_finds_each_record_any_span_of_a_query_overlaps(tmp_path: Path) -> None:
    """Test that a query of several spans yields each record once, but each repeated line too."""
    across = Bed3("chr1", start=100, end=210)
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        for record in (across, EXON1, EXON1, INTRON, EXON2):
            writer.write(record)
    with register_locatable(Read, read_blocks), TabixDetector[Bed3](path) as detector:
        assert list(detector.overlapping(SPLICED)) == [across, EXON1, EXON1, EXON2]
        assert detector.overlaps(SPLICED)
        assert not detector.overlaps(Read("chr1", ((300, 310),)))


def test_a_bed12_query_is_still_matched_by_its_whole_extent() -> None:
    """Test that a query with a refname, start and end is matched by its whole extent."""
    transcript = Bed12(
        "chr1",
        start=100,
        end=210,
        name=None,
        score=None,
        strand=None,
        thick_start=None,
        thick_end=None,
        item_rgb=None,
        block_count=2,
        block_sizes=(10, 10),
        block_starts=(0, 100),
    )
    detector = TreeDetector([INTRON])
    assert detector.overlaps(transcript)
    assert list(detector.enclosed_by(transcript)) == [INTRON]


def test_a_query_without_a_refname_is_matched_by_its_spans() -> None:
    """Test that a BED query without a refname, like a BEDPE record, is matched by its spans."""
    pair = BedPE(
        refname1="chr1",
        start1=100,
        end1=110,
        refname2="chr1",
        start2=200,
        end2=210,
        name=None,
        score=None,
        strand1=None,
        strand2=None,
    )
    detector = TreeDetector([EXON1, INTRON, EXON2])
    assert list(detector.overlapping(pair)) == [EXON1, EXON2]


def test_features_keep_the_spans_they_were_added_by() -> None:
    """Test that features stay indexed by the spans registered when they were added."""
    with register_locatable(Read, read_blocks):
        detector = TreeDetector([SPLICED])
    assert detector.overlaps(EXON1)
    assert not detector.overlaps(INTRON)
    with register_locatable(Read, read_extent):
        assert not detector.overlaps(INTRON)
    with pytest.raises(TypeError, match="register_locatable"):
        _ = detector.overlaps(SPLICED)


def test_a_lazy_query_uses_the_registration_when_it_is_consumed(tmp_path: Path) -> None:
    """Test that both detectors read a query's spans only when the query is consumed."""
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(EXON1)
    with TabixDetector[Bed3](path) as tabix:
        for detector in (TreeDetector([EXON1]), tabix):
            with register_locatable(Read, read_blocks):
                hits = detector.overlapping(SPLICED)
            with pytest.raises(TypeError, match="register_locatable"):
                _ = list(hits)
            unregistered = (detector.overlapping(SPLICED), detector.enclosing(SPLICED))
            for query in unregistered:
                with pytest.raises(TypeError, match="register_locatable"):
                    _ = list(query)


@dataclass(frozen=True)
class Mated:
    """A read and its mate, which may be on different references."""

    read: Bed3
    mate: Bed3


def mated_spans(pair: Mated) -> tuple[Bed3, Bed3]:
    """Return the spans of a read and its mate."""
    return pair.read, pair.mate


def test_a_query_with_spans_on_two_references_finds_targets_on_both(tmp_path: Path) -> None:
    """Test that both detectors find targets on each reference a query's spans are on."""
    targets = [Bed3("chr1", start=100, end=110), Bed3("chr2", start=50, end=60)]
    pair = Mated(Bed3("chr1", start=105, end=106), Bed3("chr2", start=55, end=56))
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        for target in targets:
            writer.write(target)
    with register_locatable(Mated, mated_spans), TabixDetector[Bed3](path) as tabix:
        assert list(TreeDetector(targets).overlapping(pair)) == targets
        assert list(tabix.overlapping(pair)) == targets


def test_a_bed2_query_is_matched_by_its_base(tmp_path: Path) -> None:
    """Test that both detectors match a point query, which has no end, by its one base."""
    target = Bed3("chr1", start=100, end=110)
    path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed3](path, index=IndexFormat.TBI) as writer:
        writer.write(target)
    held_by_spans = TreeDetector([target, Hollow("chr9", start=0, end=1)])
    with TabixDetector[Bed3](path) as tabix:
        for detector in (TreeDetector([target]), held_by_spans, tabix):
            assert detector.overlaps(Bed2("chr1", start=105))
            assert [hit for hit in detector.overlapping(Bed2("chr1", start=105))] == [target]
            assert list(detector.enclosing(Bed2("chr1", start=105))) == [target]
            assert not detector.overlaps(Bed2("chr1", start=110))
            assert not detector.overlaps(Bed2("chr2", start=105))


def test_records_of_every_locatable_kind_go_in_together() -> None:
    """Test that BED records, records with spans() and registered records mix as Locatable."""
    records: list[Locatable] = [EXON1, Located((EXON2,)), SPLICED]
    with register_locatable(Read, read_blocks):
        detector: TreeDetector[Locatable] = TreeDetector(records)
        assert set(detector.overlapping(Bed3("chr1", start=0, end=1000))) == set(records)
        assert list(Territory(records)) == [
            Bed3("chr1", start=100, end=110),
            Bed3("chr1", start=200, end=210),
        ]


def test_locatable_is_a_type_alias_and_not_a_class() -> None:
    """Test that Locatable cannot be mistaken for a class: isinstance and subclassing refuse it."""
    alias: Any = Locatable
    with pytest.raises(TypeError):
        _ = isinstance(SPLICED, alias)
    with pytest.raises(TypeError):
        _ = type("Mine", (alias,), {})


def test_a_plain_record_beside_a_registered_one_is_found_by_its_own_span() -> None:
    """Test that a record with a refname, start and end is its own span beside registered ones."""
    plain = Plain("chr1", start=300, end=310)
    with register_locatable(Read, read_blocks):
        detector: TreeDetector[Locatable] = TreeDetector([SPLICED, plain])
        assert list(detector.overlapping(Bed3("chr1", start=305, end=306))) == [plain]
        assert list(detector.overlapping(EXON1)) == [SPLICED]


@dataclass(frozen=True)
class Tagged(Bed6):
    """A kind of BED6 record, which is found by its own spans."""


@pytest.mark.parametrize("kind", [Bed3, Bed12, Tagged])
def test_a_bed_class_cannot_be_registered(kind: Any) -> None:
    """Test that registering a BED class raises, since a BED record is found by its own spans."""
    with pytest.raises(TypeError, match=r"override spans\(\) in a subclass"):
        _ = register_locatable(kind, read_extent)


class Mixin:
    """A class of records registered for its subclasses, one of them a BED record."""


@dataclass(frozen=True)
class MixedBed(Mixin, Bed3):
    """A BED record whose parent class is registered."""


def test_a_bed_record_is_found_by_its_own_spans_beside_a_registered_parent() -> None:
    """Test that a BED record is found by its own spans even when a parent class is registered."""
    record = MixedBed("chr1", start=100, end=110)
    with register_locatable(Mixin, lambda _: [Bed3("chr2", start=0, end=10)]):
        assert TreeDetector([record]).overlaps(EXON1)
        assert TreeDetector([EXON1]).overlaps(record)
        assert list(Territory([record])) == [Bed3("chr1", start=100, end=110)]


def test_only_a_function_can_be_registered() -> None:
    """Test that registering something that is not a function raises."""
    not_a_function: Any = None
    with pytest.raises(TypeError, match="must be a function returning spans"):
        _ = register_locatable(Read, not_a_function)
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([SPLICED])


@pytest.mark.parametrize(
    ("bad", "error"), [(("chr1", 20, 30), AttributeError), (Plain("chr1", 30, 20), ValueError)]
)
def test_adding_features_adds_none_when_a_registered_function_returns_a_bad_span(
    bad: object, error: type[Exception]
) -> None:
    """Test that a bad span adds none of the features, leaving the detector as it was."""
    detector: TreeDetector[Bed3 | Read] = TreeDetector([EXON1])

    def spans(_: Read) -> Any:
        return [Bed3("chr1", start=0, end=10), bad]

    with register_locatable(Read, spans), pytest.raises(error):
        detector.add(INTRON, SPLICED)
    assert list(detector) == [EXON1]
    assert list(detector.overlapping(Bed3("chr1", start=0, end=300))) == [EXON1]
    detector.add(EXON2)
    assert set(detector.overlapping(Bed3("chr1", start=0, end=300))) == {EXON1, EXON2}


class Disguised(Plain):
    """A record with a refname, start and end whose `__class__` claims it is a `Read`."""

    @override
    def __getattribute__(self, name: str) -> Any:
        """Claim to be a `Read`."""
        return Read if name == "__class__" else super().__getattribute__(name)


def test_a_record_is_located_by_its_type_when_added_and_when_queried() -> None:
    """Test that adding and querying both locate a record by its type, not by its `__class__`."""
    disguised = Disguised("chr1", start=100, end=110)
    with register_locatable(Read, read_extent):
        assert TreeDetector([disguised]).overlaps(EXON1)
        assert TreeDetector([EXON1]).overlaps(disguised)
        assert list(Territory([disguised])) == [Bed3("chr1", start=100, end=110)]


def test_a_span_of_a_registered_class_is_used_as_it_is() -> None:
    """Test that a span a registered function returns is not located again by its own class."""
    read = Read("chr1", ((100, 110),))
    base_before = Bed3("chr1", start=99, end=100)

    def as_one_based(read: Read) -> Any:
        return [OneBased("chr1", *read.blocks[0])]

    with register_locatable(OneBased, one_based), register_locatable(Read, as_one_based):
        assert TreeDetector([base_before]).overlaps(OneBased("chr1", 100, 110))
        assert not TreeDetector([base_before]).overlaps(read)
        assert not TreeDetector([read]).overlaps(base_before)
