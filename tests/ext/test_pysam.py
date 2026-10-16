"""Tests of the functions that read pysam's records as spans."""

# pyright: reportUnknownMemberType=false

from io import StringIO
from pathlib import Path
from typing import Any

import pysam
import pytest
from pybgzf import IndexFormat
from pysam.libctabixproxies import GTFProxy

import bedspec.ext.pysam
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec import ReferenceSpan
from bedspec import Territory
from bedspec import register_locatable
from bedspec.ext.pysam import as_blocks
from bedspec.ext.pysam import as_interval
from bedspec.ext.pysam import as_paired_blocks
from bedspec.ext.pysam import as_paired_intervals
from bedspec.overlap import TabixDetector
from bedspec.overlap import TreeDetector

HEADER = pysam.AlignmentHeader.from_dict({
    "SQ": [{"SN": "chr1", "LN": 10_000}, {"SN": "chr2", "LN": 10_000}]
})


def read(start: int | None, cigar: str = "10M", reverse: bool = False) -> pysam.AlignedSegment:
    """Return a read aligned to chr1 at a 0-based start, or unmapped if the start is None."""
    segment = pysam.AlignedSegment(HEADER)
    segment.query_name = "read"
    if start is None:
        segment.is_unmapped = True
        return segment
    segment.reference_name = "chr1"
    segment.reference_start = start
    segment.cigarstring = cigar
    segment.is_reverse = reverse
    return segment


def mated(
    segment: pysam.AlignedSegment,
    refname: str,
    start: int,
    cigar: str | None,
    reverse: bool = False,
) -> None:
    """Pair a read with a mate at a 0-based start, its CIGAR in the MC tag, or unmapped if None."""
    segment.is_paired = True
    segment.next_reference_name = refname
    segment.next_reference_start = start
    segment.mate_is_reverse = reverse
    if cigar is None:
        segment.mate_is_unmapped = True
    else:
        segment.set_tag("MC", cigar)


def stranded(start: int, end: int, strand: BedStrand, refname: str = "chr1") -> Bed6:
    """Return a span on a reference, chr1 by default, on a strand."""
    return Bed6(refname, start=start, end=end, name=None, score=None, strand=strand)


def tabix_records(lines: list[str], parser: Any) -> list[Any]:
    """Return the records pysam's tabix parser reads from lines of text."""
    generic: Any = pysam.tabix_generic_iterator
    return list(generic(StringIO("".join(f"{line}\n" for line in lines)), parser=parser))


def test_as_blocks_splits_an_alignment_at_its_introns_and_deletions() -> None:
    """Test that a read's blocks skip its introns and deletions, on its strand."""
    assert as_blocks(read(100, "10M500N10M2D5M")) == [
        stranded(100, 110, BedStrand.Positive),
        stranded(610, 620, BedStrand.Positive),
        stranded(622, 627, BedStrand.Positive),
    ]
    assert as_blocks(read(100, "5S10M", reverse=True)) == [stranded(100, 110, BedStrand.Negative)]


def test_as_interval_of_an_alignment_is_its_whole_alignment() -> None:
    """Test that a read's span runs from its first aligned base to its last."""
    assert as_interval(read(100, "10M500N10M")) == [stranded(100, 620, BedStrand.Positive)]
    assert as_interval(read(100, reverse=True)) == [stranded(100, 110, BedStrand.Negative)]


def test_an_unmapped_read_has_no_spans() -> None:
    """Test that an unmapped read, placed beside its mate or not, has no spans."""
    placed = read(None)
    placed.reference_name, placed.reference_start = "chr1", 50
    for unmapped in (read(None), placed):
        assert as_blocks(unmapped) == []
        assert as_interval(unmapped) == []
        assert as_paired_blocks(unmapped) == []
        assert as_paired_intervals(unmapped) == []


def test_as_paired_blocks_adds_the_blocks_of_the_mate_from_its_mate_cigar() -> None:
    """Test that a read's pair blocks are its own and its mate's, each on its own strand."""
    first = read(100, "10M500N10M")
    mated(first, "chr2", 900, "5S10M2D3I12M", reverse=True)
    assert as_paired_blocks(first) == [
        stranded(100, 110, BedStrand.Positive),
        stranded(610, 620, BedStrand.Positive),
        stranded(900, 910, BedStrand.Negative, refname="chr2"),
        stranded(912, 924, BedStrand.Negative, refname="chr2"),
    ]


def test_the_paired_spans_of_a_read_without_a_mapped_mate_are_its_own() -> None:
    """Test that an unpaired read, or one whose mate is unmapped, has only its own blocks."""
    alone = read(100)
    orphan = read(100)
    mated(orphan, "chr1", 100, None)
    for segment in (alone, orphan):
        assert as_paired_blocks(segment) == [stranded(100, 110, BedStrand.Positive)]
        assert as_paired_intervals(segment) == [stranded(100, 110, BedStrand.Positive)]


def test_the_paired_spans_of_an_unmapped_read_are_its_mates() -> None:
    """Test that an unmapped read with a mapped mate is found by its mate's blocks."""
    unmapped = read(None)
    mated(unmapped, "chr1", 300, "20M", reverse=True)
    assert as_paired_blocks(unmapped) == [stranded(300, 320, BedStrand.Negative)]
    assert as_paired_intervals(unmapped) == [stranded(300, 320, BedStrand.Negative)]


def test_the_paired_spans_without_a_mate_cigar_raise() -> None:
    """Test that a read with a mapped mate and no MC tag raises, naming the fix."""
    segment = read(100)
    mated(segment, "chr1", 300, "20M")
    segment.set_tag("MC", None)
    for paired in (as_paired_blocks, as_paired_intervals):
        with pytest.raises(ValueError, match="samtools fixmate -m"):
            _ = paired(segment)


def test_as_paired_intervals_adds_the_interval_of_the_mate_from_its_mate_cigar() -> None:
    """Test that a read's paired intervals are its own and its mate's, each on its own strand."""
    first = read(100, "10M500N10M", reverse=True)
    mated(first, "chr2", 900, "5S10M2D3I12M")
    assert as_paired_intervals(first) == [
        stranded(100, 620, BedStrand.Negative),
        stranded(900, 924, BedStrand.Positive, refname="chr2"),
    ]


def test_as_paired_intervals_never_joins_mates_on_one_reference() -> None:
    """Test that mates on one reference keep an interval each, the insert between them unheld."""
    first = read(100, "20M")
    mated(first, "chr1", 300, "20M", reverse=True)
    assert as_paired_intervals(first) == [
        stranded(100, 120, BedStrand.Positive),
        stranded(300, 320, BedStrand.Negative),
    ]


def test_as_interval_of_a_variant_covers_its_reference_allele() -> None:
    """Test that a variant spans the bases of its reference allele, from 0-based coordinates."""
    header = pysam.VariantHeader()
    header.add_line("##contig=<ID=chr1,length=1000>")
    deletion = header.new_record(contig="chr1", start=99, alleles=("CAG", "C"))
    assert as_interval(deletion) == [Bed3("chr1", start=99, end=102)]


@pytest.mark.parametrize(
    ("line", "parser", "strand"),
    [
        ('chr1\tsrc\texon\t11\t20\t.\t-\t.\tgene_id "g";', pysam.asGTF(), BedStrand.Negative),
        ("chr1\tsrc\tgene\t11\t20\t.\t+\t.\tID=g", pysam.asGFF3(), BedStrand.Positive),
        ("chr1\tsrc\tgene\t11\t20\t.\t.\t.\tID=g", pysam.asGFF3(), None),
        ("chr1\t10\t20\tname\t0\t-", pysam.asBed(), BedStrand.Negative),
        ("chr1\t10\t20", pysam.asBed(), None),
    ],
)
def test_as_interval_reads_gtf_gff3_and_bed_lines(
    line: str, parser: Any, strand: BedStrand | None
) -> None:
    """Test that a tabix line's span is 0-based, on its strand, or on none if it has none."""
    (record,) = tabix_records([line], parser)
    assert as_interval(record) == [
        Bed6("chr1", start=10, end=20, name=None, score=None, strand=strand)
    ]


def test_registered_reads_are_found_by_their_blocks_and_not_their_introns() -> None:
    """Test that reads registered by their blocks are found by them in a tree detector."""
    spliced, other = read(100, "10M500N10M"), read(1000)
    with register_locatable(pysam.AlignedSegment, as_blocks):
        detector = TreeDetector([spliced, other])
        assert list(detector.overlapping(Bed3("chr1", start=615, end=616))) == [spliced]
        assert not detector.overlaps(Bed3("chr1", start=300, end=310))
        assert list(detector) == [spliced, other]


def test_registered_reads_query_targets_by_their_blocks_and_strand() -> None:
    """Test that reads registered by their blocks find targets by them, on their strand."""
    plus, minus = stranded(100, 200, BedStrand.Positive), stranded(600, 700, BedStrand.Negative)
    targets = TreeDetector([plus, minus])
    with register_locatable(pysam.AlignedSegment, as_blocks):
        assert list(targets.overlapping(read(150, "10M500N10M"))) == [plus, minus]
        assert list(targets.overlapping(read(150, "10M500N10M"), stranded=True)) == [plus]
        assert list(targets.overlapping(read(None))) == []
    in_intron = TreeDetector([Bed3("chr1", start=300, end=310)])
    with register_locatable(pysam.AlignedSegment, as_interval):
        assert in_intron.overlaps(read(250, "10M300N10M"))
    with register_locatable(pysam.AlignedSegment, as_blocks):
        assert not in_intron.overlaps(read(250, "10M300N10M"))


def test_a_territory_of_registered_reads_holds_their_blocks() -> None:
    """Test that a territory of reads registered by their blocks holds their blocks only."""
    with register_locatable(pysam.AlignedSegment, as_blocks):
        territory = Territory([read(100, "10M500N10M"), read(None)])
        assert list(territory) == [
            Bed3("chr1", start=100, end=110),
            Bed3("chr1", start=610, end=620),
        ]
        assert read(None) not in territory
        assert read(102, "5M") in territory


def test_registered_tabix_lines_and_variants_are_found() -> None:
    """Test that tabix lines and variants registered by their spans are found."""
    (exon,) = tabix_records(['chr1\tsrc\texon\t11\t20\t.\t-\t.\tgene_id "g";'], pysam.asGTF())
    header = pysam.VariantHeader()
    header.add_line("##contig=<ID=chr1,length=1000>")
    snv = header.new_record(contig="chr1", start=14, alleles=("C", "T"))
    with (
        register_locatable(GTFProxy, as_interval),
        register_locatable(pysam.VariantRecord, as_interval),
    ):
        detector = TreeDetector([exon, snv])
        assert list(detector.overlapping(stranded(14, 15, BedStrand.Negative), stranded=True)) == [
            exon
        ]
        assert {id(hit) for hit in detector.overlapping(Bed3("chr1", start=14, end=15))} == {
            id(exon),
            id(snv),
        }


def test_the_readme_example_keeps_and_tags_both_mates_of_pairs_on_target(tmp_path: Path) -> None:
    """Test that the README's example keeps both mates of a pair that is on target."""
    targets_path = tmp_path / "targets.bed.gz"
    with BedWriter.from_path[Bed4](targets_path, index=IndexFormat.TBI) as writer:
        writer.write(Bed4("chr1", start=1000, end=1100, name="exon2"))
        writer.write(Bed4("chr1", start=1050, end=1060, name="exon1"))

    reads: list[pysam.AlignedSegment] = []
    for name, (start1, start2) in {"on": (100, 1040), "off": (200, 500)}.items():
        for first, (start, mate_start) in enumerate([(start1, start2), (start2, start1)]):
            segment = read(start, "20M", reverse=bool(first))
            segment.query_name = name
            segment.is_read1, segment.is_read2 = not first, bool(first)
            mated(segment, "chr1", mate_start, "20M")
            reads.append(segment)
    sorted_path, kept_path = tmp_path / "sorted.bam", tmp_path / "on-target.bam"
    with pysam.AlignmentFile(str(sorted_path), "wb", header=HEADER) as out:
        for segment in sorted(reads, key=lambda segment: segment.reference_start):
            _ = out.write(segment)

    with (
        register_locatable(pysam.AlignedSegment, as_paired_intervals),
        TabixDetector[Bed4](targets_path) as targets,
        pysam.AlignmentFile(str(sorted_path)) as bam,
        pysam.AlignmentFile(str(kept_path), "wb", template=bam) as kept,
    ):
        for segment in bam:
            if names := {str(target.name) for target in targets.overlapping(segment)}:
                segment.set_tag("tg", ",".join(sorted(names)))
                _ = kept.write(segment)

    with pysam.AlignmentFile(str(kept_path)) as bam:
        assert [(segment.query_name, segment.get_tag("tg")) for segment in bam] == [
            ("on", "exon1,exon2"),
            ("on", "exon1,exon2"),
        ]


def test_an_unregistered_tabix_line_raises_a_type_error_not_a_key_error() -> None:
    """Test that a tabix line, whose missing attributes raise KeyError, is refused clearly."""
    (exon,) = tabix_records(['chr1\tsrc\texon\t11\t20\t.\t-\t.\tgene_id "g";'], pysam.asGTF())
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([exon])
    with pytest.raises(TypeError, match="register_locatable"):
        _ = TreeDetector([Bed3("chr1", start=0, end=10)]).overlaps(exon)
    assert exon not in Territory([Bed3("chr1", start=0, end=100)])


def test_bedspec_ext_pysam_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.pysam holds its functions and none of the names it imports."""
    module = bedspec.ext.pysam
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert (
        public
        == sorted(module.__all__)
        == ["as_blocks", "as_interval", "as_paired_blocks", "as_paired_intervals"]
    )


@pytest.mark.parametrize("cigar", ["50M2I50M", "50=1X49=", "10S50M2I50M10S"])
def test_as_blocks_joins_blocks_that_abut_across_an_insertion_or_mismatch(cigar: str) -> None:
    """Test that an insertion or a mismatch operator does not split an alignment's blocks."""
    assert as_blocks(read(100, cigar)) == [stranded(100, 200, BedStrand.Positive)]


def test_as_paired_blocks_joins_the_blocks_of_a_mate_across_an_insertion() -> None:
    """Test that an insertion in the mate's CIGAR does not split the mate's blocks."""
    first = read(100, "20M")
    mated(first, "chr1", 300, "50M2I50=1X49=")
    assert as_paired_blocks(first) == [
        stranded(100, 120, BedStrand.Positive),
        stranded(300, 450, BedStrand.Positive),
    ]


def test_an_alignment_with_an_insertion_is_one_span_as_a_query() -> None:
    """Test that an alignment with an insertion can be asked what encloses it."""
    target = Bed3("chr1", start=0, end=1_000)
    with register_locatable(pysam.AlignedSegment, as_blocks):
        assert list(TreeDetector([target]).enclosing(read(100, "50M2I50M"))) == [target]


@pytest.mark.parametrize("cigar", ["10S", "5I"])
def test_an_alignment_on_no_reference_bases_has_no_interval_and_no_blocks(cigar: str) -> None:
    """Test that a mapped read whose CIGAR consumes no reference bases is found by nothing."""
    assert as_interval(read(100, cigar)) == []
    assert as_blocks(read(100, cigar)) == []


def test_a_mate_without_a_cigar_adds_nothing() -> None:
    """Test that a mate CIGAR of * adds no interval and no blocks for the mate."""
    segment = read(100)
    mated(segment, "chr1", 300, "*")
    assert as_paired_intervals(segment) == [stranded(100, 110, BedStrand.Positive)]
    assert as_paired_blocks(segment) == [stranded(100, 110, BedStrand.Positive)]


def test_a_mate_at_no_position_adds_nothing() -> None:
    """Test that a mate at no position adds no interval and no blocks, rather than raising."""
    segment = read(100)
    mated(segment, "chr1", 300, "20M")
    segment.next_reference_start = -1
    assert as_paired_intervals(segment) == [stranded(100, 110, BedStrand.Positive)]
    assert as_paired_blocks(segment) == [stranded(100, 110, BedStrand.Positive)]


def test_pairs_touching_the_panel_are_kept_and_the_rest_written_as_a_territory(
    tmp_path: Path,
) -> None:
    """Test that pairs touching a panel are kept, and the rest written as a sorted, indexed BED."""
    panel_path = tmp_path / "panel.bed.gz"
    with BedWriter.from_path[Bed4](panel_path, index=IndexFormat.TBI) as writer:
        writer.write(Bed4("chr1", start=1_000, end=1_100, name="BRCA1_ex2"))
        writer.write(Bed4("chr1", start=5_000, end=5_100, name="BRCA1_ex3"))

    on1, on2 = read(1_050, "50M"), read(1_300, "50M", reverse=True)
    off1, off2 = read(20_000, "50M"), read(20_300, "50M", reverse=True)
    half1, half2 = read(5_020, "50M"), read(None)
    half2.reference_name, half2.reference_start = "chr1", 5_020
    lost1, lost2 = read(None), read(None)
    mated(on1, "chr1", 1_300, "50M")
    mated(on2, "chr1", 1_050, "50M")
    mated(off1, "chr1", 20_300, "50M")
    mated(off2, "chr1", 20_000, "50M")
    mated(half1, "chr1", 5_020, None)
    mated(half2, "chr1", 5_020, "50M")
    for name, segments in {"on": (on1, on2), "off": (off1, off2), "half": (half1, half2)}.items():
        for segment in segments:
            segment.query_name = name
    for segment in (lost1, lost2):
        segment.query_name, segment.is_paired, segment.mate_is_unmapped = "lost", True, True
    bam_path = tmp_path / "sorted.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=HEADER) as out:
        for segment in (on1, on2, half1, half2, off1, off2, lost1, lost2):
            _ = out.write(segment)

    with (
        register_locatable(pysam.AlignedSegment, as_paired_intervals),
        TabixDetector[Bed4](panel_path) as panel,
        pysam.AlignmentFile(str(bam_path)) as bam,
        pysam.AlignmentFile(str(tmp_path / "on-target.bam"), "wb", template=bam) as on_target,
        BedWriter.from_path[Bed3](tmp_path / "off-target.bed.gz", index=IndexFormat.TBI) as writer,
    ):
        off_target: list[ReferenceSpan] = []
        for segment in bam:
            if names := {str(target.name) for target in panel.overlapping(segment)}:
                segment.set_tag("tg", ",".join(sorted(names)))
                _ = on_target.write(segment)
            else:
                off_target.extend(as_interval(segment))
        writer.write_all(Territory(off_target))

    with pysam.AlignmentFile(str(tmp_path / "on-target.bam")) as kept:
        assert [(segment.query_name, segment.get_tag("tg")) for segment in kept] == [
            ("on", "BRCA1_ex2"),
            ("on", "BRCA1_ex2"),
            ("half", "BRCA1_ex3"),
            ("half", "BRCA1_ex3"),
        ]
    with TabixDetector[Bed3](tmp_path / "off-target.bed.gz") as written:
        assert list(written.overlapping(Bed3("chr1", start=0, end=10_000))) == []
        assert list(written.overlapping(Bed3("chr1", start=10_000, end=30_000))) == [
            Bed3("chr1", start=20_000, end=20_050),
            Bed3("chr1", start=20_300, end=20_350),
        ]


@pytest.mark.parametrize("paired", [as_paired_intervals, as_paired_blocks])
def test_a_stranded_query_finds_a_pair_by_each_mate_on_its_own_strand(paired: Any) -> None:
    """Test that a registered pair is found by whichever mate is on a target's strand."""
    first = read(100, "20M")
    mated(first, "chr1", 300, "20M", reverse=True)
    on_their_strands = [
        stranded(100, 120, BedStrand.Positive),
        stranded(300, 320, BedStrand.Negative),
    ]
    on_the_other_strands = [
        stranded(100, 120, BedStrand.Negative),
        stranded(300, 320, BedStrand.Positive),
    ]
    with register_locatable(pysam.AlignedSegment, paired):
        assert (
            list(TreeDetector(on_their_strands).overlapping(first, stranded=True))
            == on_their_strands
        )
        assert not TreeDetector(on_the_other_strands).overlaps(first, stranded=True)


@pytest.mark.parametrize("mate_cigar", ["50S", "garbage"])
def test_a_mate_whose_cigar_holds_no_reference_base_adds_nothing(mate_cigar: str) -> None:
    """Test that a mate CIGAR consuming no reference base adds no span, as for the read itself."""
    first = read(100, "20M")
    mated(first, "chr1", 500, mate_cigar)
    assert as_paired_intervals(first) == [stranded(100, 120, BedStrand.Positive)]
    assert as_paired_blocks(first) == [stranded(100, 120, BedStrand.Positive)]


def test_the_intervals_of_reads_are_typed_as_stranded_bed_records() -> None:
    """Test that a read's intervals type as BED6 records, so they can go where BED6 is expected."""
    first = read(100, "20M")
    mated(first, "chr1", 300, "20M", reverse=True)
    interval: list[Bed6] = as_interval(first)
    paired: list[Bed6] = as_paired_intervals(first)
    assert interval == [stranded(100, 120, BedStrand.Positive)]
    assert paired == [*interval, stranded(300, 320, BedStrand.Negative)]
