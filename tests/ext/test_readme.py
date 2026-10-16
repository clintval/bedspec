"""Run every README example, including those that read files, against files made for them."""

# pyright: reportUnknownMemberType=false

import gzip
import os
from collections import Counter
from pathlib import Path
from typing import Any

import pysam
from pybgzf import IndexFormat

from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed6
from bedspec import BedStrand
from bedspec import BedWriter
from bedspec.overlap import TabixDetector

README = Path(__file__).parents[2] / "README.md"

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6"},
    "SQ": [{"SN": "chr2", "LN": 10_000}, {"SN": "chr1", "LN": 10_000}],
})
"""A header that lists chr2 before chr1, so sorting by its order differs from sorting by name."""


def mate(
    name: str, start: int, cigar: str, mate_start: int, mate_cigar: str, first: bool
) -> pysam.AlignedSegment:
    """Return one mate of a pair on chr1, with its mate's CIGAR in its MC tag."""
    segment = pysam.AlignedSegment(HEADER)
    segment.query_name = name
    segment.reference_name = "chr1" if name != "off2" else "chr2"
    segment.reference_start = start
    segment.cigarstring = cigar
    length = segment.infer_query_length() or 0
    segment.query_sequence = "A" * length
    segment.query_qualities = pysam.qualitystring_to_array("I" * length)
    segment.mapping_quality = 60
    segment.is_paired = True
    segment.is_read1, segment.is_read2 = first, not first
    segment.is_reverse, segment.mate_is_reverse = not first, first
    segment.next_reference_name = segment.reference_name
    segment.next_reference_start = mate_start
    segment.set_tag("MC", mate_cigar)
    return segment


def pair(name: str, start1: int, start2: int, cigar1: str = "50M") -> list[pysam.AlignedSegment]:
    """Return a pair whose second mate aligns 50 bases at `start2`."""
    return [
        mate(name, start1, cigar1, start2, "50M", first=True),
        mate(name, start2, "50M", start1, cigar1, first=False),
    ]


def unmapped(name: str) -> list[pysam.AlignedSegment]:
    """Return a pair whose mates are both unmapped."""
    segments: list[pysam.AlignedSegment] = []
    for first in (True, False):
        segment = pysam.AlignedSegment(HEADER)
        segment.query_name = name
        segment.query_sequence = "A" * 50
        segment.query_qualities = pysam.qualitystring_to_array("I" * 50)
        segment.is_paired = segment.is_unmapped = segment.mate_is_unmapped = True
        segment.is_read1, segment.is_read2 = first, not first
        segments.append(segment)
    return segments


def write_files() -> None:
    """Write every file the README examples read into the current directory."""
    grouped = [
        *pair("on1", 1_020, 1_300),
        *pair("half", 1_090, 5_000),
        *pair("off1", 3_000, 3_200),
        *pair("off2", 7_000, 7_400),
        *pair("spliced", 900, 2_000, cigar1="20M500N30M"),
        *unmapped("lost"),
    ]
    with pysam.AlignmentFile("grouped.bam", "wb", header=HEADER) as out:
        for segment in grouped:
            _ = out.write(segment)
    with pysam.AlignmentFile("unsorted.bam", "wb", header=HEADER) as out:
        for segment in reversed(grouped):
            _ = out.write(segment)
    _ = pysam.sort("-o", "sorted.bam", "unsorted.bam")
    _ = pysam.sort("-o", "rna.bam", "unsorted.bam")
    for name in ("panel.bed.gz", "targets.bed.gz"):
        with BedWriter.from_path[Bed4](name, index=IndexFormat.TBI) as writer:
            writer.write(Bed4("chr2", start=500, end=600, name="T2"))
            writer.write(Bed4("chr1", start=1_000, end=1_100, name="T1"))
    with BedWriter.from_path[Bed6]("exons.bed.gz", index=IndexFormat.TBI) as writer:
        for start, end in ((900, 920), (1_420, 1_450)):
            writer.write(
                Bed6("chr1", start=start, end=end, name=None, score=None, strand=BedStrand.Positive)
            )
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=10_000)
    with pysam.VariantFile("calls.vcf.gz", "wz", header=header) as out:
        for start in (1_050, 1_150, 4_000):
            _ = out.write(
                out.new_record(contig="chr1", start=start, stop=start + 1, alleles=("A", "T"))
            )
    _ = Path("genes.gff").write_text("##gff-version 3\nchr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n")


def readme_blocks() -> list[tuple[int, str]]:
    """Return the line and code of every Python block in the README, skipped by sybil or not."""
    lines = README.read_text().split("\n")
    starts = [index for index, line in enumerate(lines) if line.startswith("```python")]
    return [
        (start + 1, "\n".join(lines[start + 1 : lines.index("```", start + 1)])) for start in starts
    ]


def test_every_readme_example_runs_against_the_files_it_reads(
    tmp_path: Path, monkeypatch: Any
) -> None:
    """Test that every README example runs in order, as sybil runs them, reading real files."""
    monkeypatch.chdir(tmp_path)
    write_files()
    namespace: dict[str, object] = {"__name__": "__readme__"}
    for line, code in readme_blocks():
        exec(compile(code, f"{README}:{line}", "exec"), namespace)  # noqa: S102
    with pysam.AlignmentFile("on-target.bam") as kept:
        assert Counter(segment.query_name for segment in kept) == {
            "on1": 2,
            "half": 2,
            "spliced": 2,
        }
    with pysam.AlignmentFile("kept.bam") as kept:
        assert Counter(segment.query_name for segment in kept) == {"on1": 2, "half": 2}
    with gzip.open("off-target.bed.gz", "rt") as written:
        assert written.read().splitlines() == [
            "chr2\t7000\t7050",
            "chr2\t7400\t7450",
            "chr1\t3000\t3050",
            "chr1\t3200\t3250",
        ]
    with TabixDetector[Bed3]("off-target.bed.gz") as detector:
        assert detector.overlaps(Bed3("chr1", start=3_010, end=3_020))
    assert os.path.exists("off-target.bed.gz.tbi")
