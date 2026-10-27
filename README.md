# bedspec

[![PyPi Release](https://badge.fury.io/py/bedspec.svg)](https://badge.fury.io/py/bedspec)
[![CI](https://github.com/clintval/bedspec/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/bedspec/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12_|_3.13_|_3.14-blue)](https://github.com/clintval/typeline)
[![basedpyright](https://img.shields.io/badge/basedpyright-checked-42b983)](https://docs.basedpyright.com/latest/)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://docs.astral.sh/uv/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

An HTS-specs compliant BED toolkit.

## Installation

The package can be installed with `pip`:

```console
pip install bedspec
```

To sort more records than fit in memory with spillsort, install the `sort` extra:

```console
pip install 'bedspec[sort]'
```

## At a Glance

From a BAM in any order, write every read whose pair misses a capture panel as a BED line, sorted on disk in the BAM's reference order and indexed with tabix:

<!-- skip: next -->
```python
from pybgzf import IndexFormat
from pysam import AlignedSegment, AlignmentFile
from spillsort import Sorter

from bedspec import Bed3, Bed4, BedWriter, ReferenceOrder, register_locatable
from bedspec.ext.pysam import as_paired_intervals
from bedspec.overlap import TabixDetector

with (
    register_locatable(AlignedSegment, as_paired_intervals),
    AlignmentFile("unsorted.bam") as reader,
    TabixDetector[Bed4]("panel.bed.gz") as panel,
    Sorter[Bed3](ReferenceOrder(reader.references)) as off_target_sorter,
    BedWriter.from_path[Bed3]("off-target.bed.gz", index=IndexFormat.TBI) as writer,
):
    for alignment in reader:
        if alignment.is_mapped and not panel.overlaps(alignment):
            bed = Bed3(
                alignment.reference_name,
                start=alignment.reference_start,
                end=alignment.reference_end,
            )
            off_target_sorter.add(bed)

    writer.write_all(off_target_sorter.sort())
```

An alignment pair does not overlap the panel when neither mate's interval overlaps any region in the panel.
An unsorted alignment file can be queried this way only if `samtools fixmate -m` has pre-added the standard `MC` tag upon the records.
This example also shows how, with the `spillsort` library installed, you can write the non-overlapping alignments as BED3 intervals in coordinate sorted order to a block-compressed and tabix indexed output file.

## Quickstart

### Building a BED Feature

```python
from bedspec import Bed3

bed = Bed3("chr1", start=2, end=8)
```

Records are checked against the BED spec when they are built.
BED has no quoting, so text is written and read as it is, and a value holding a tab is refused.
A feature may start where it ends, as an insertion does.

Records are immutable and hashable.
Use `dataclasses.replace` to build a changed copy, which is checked like any other record:

```python
from dataclasses import replace

assert replace(bed, end=10) == Bed3("chr1", start=2, end=10)
```

### Writing

```python
from tempfile import NamedTemporaryFile

from bedspec import BedWriter

temp_file = NamedTemporaryFile(mode="w+t", suffix=".txt")

with BedWriter.from_path[Bed3](temp_file.name) as writer:
    writer.write(bed)
```

### Reading

```python
from bedspec import BedReader

with BedReader.from_path[Bed3](temp_file.name) as reader:
    beds = list(reader)

assert beds == [Bed3("chr1", start=2, end=8)]
```

### Compressed and Indexed BED

A path ending in `.gz`, `.bgz`, or `.bgzf` is written as BGZF, which any gzip reader can read, and a compressed file is read by its contents.
Ask for a tabix or CSI index to have one written beside the file, with features sorted by reference and start.

```python
from pybgzf import IndexFormat

with BedWriter.from_path[Bed3](f"{temp_file.name}.gz", index=IndexFormat.TBI, threads=4) as writer:
    writer.write(Bed3("chr1", start=2, end=8))
    writer.write(Bed3("chr1", start=6, end=9))

with BedReader.from_path[Bed3](f"{temp_file.name}.gz") as reader:
    beds = list(reader)

assert beds == [Bed3("chr1", start=2, end=8), Bed3("chr1", start=6, end=9)]
```

Query an indexed file on disk with the same operations as the [overlap detector](#overlap-detection):

```python
from bedspec.overlap import TabixDetector

with TabixDetector[Bed3](f"{temp_file.name}.gz") as detector:
    enclosing = list(detector.enclosing(Bed3("chr1", start=7, end=8)))

assert enclosing == [Bed3("chr1", start=2, end=8), Bed3("chr1", start=6, end=9)]
```

No index query returns a zero-length feature at the start of a reference, so writing one to an indexed file warns.

### BED Types

This package provides builtin classes for the following BED formats:

```python
from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import Bed9
from bedspec import Bed12
from bedspec import BedGraph
from bedspec import BedPE
```

It also provides the ENCODE peak formats:

```python
from bedspec import BroadPeak
from bedspec import GappedPeak
from bedspec import NarrowPeak
```

ENCODE writes -1 for a p-value, q-value, or summit that is not given, and so do these types.

For BED files with extra columns (BEDn+m), use `Bed3N`, `Bed4N`, `Bed5N`, `Bed6N`, `Bed9N`, or `Bed12N`.
Each is its BED type plus an `extra` field that keeps any further columns as text.

```python
from pathlib import Path

from bedspec import Bed6N

Path(temp_file.name).write_text("chr1\t5\t9\tpeak\t7\t-\t3.2\t0.01\n")

with BedReader.from_path[Bed6N](temp_file.name) as reader:
    beds = list(reader)

assert beds[0].name == "peak"
assert beds[0].extra == ("3.2", "0.01")
```

### Spans and Territory

Every BED record yields its pieces on reference sequences with `spans()`, each with its strand.
A `Bed2` yields its one base, a `Bed3` to `Bed9` yields itself, a `BedPE` yields both of its ends, and a `Bed12` yields its blocks, so its introns are not part of it:

```python
from bedspec import BedStrand

transcript = Bed12(
    "chr1",
    start=100,
    end=200,
    name="tx",
    score=None,
    strand=BedStrand.Positive,
    thick_start=None,
    thick_end=None,
    item_rgb=None,
    block_count=2,
    block_sizes=(10, 20),
    block_starts=(0, 80),
)

pieces = [(span.refname, span.start, span.end, span.strand) for span in transcript.spans()]

assert pieces == [
    ("chr1", 100, 110, BedStrand.Positive),
    ("chr1", 180, 200, BedStrand.Positive),
]
```

The `territory()` method returns the bases those spans hold, without their strands, as a [`Territory`](#merging-intersecting-and-subtracting):

```python
from bedspec import Territory

assert transcript.territory() == Territory([
    Bed3("chr1", start=100, end=110),
    Bed3("chr1", start=180, end=200),
])
assert transcript.territory().length == 30
```

The `thick_territory()` method returns the bases of a BED9 or BED12 record within its thick part, which for a transcript is its coding sequence, and no bases when the thick part is unset:

```python
from dataclasses import replace

coding = replace(transcript, thick_start=105, thick_end=190)

assert coding.thick_territory() == Territory([
    Bed3("chr1", start=105, end=110),
    Bed3("chr1", start=180, end=190),
])
assert transcript.thick_territory() == Territory([])
```

### Overlap Detection

Use a fast overlap detector for any collection of interval types, including third-party:

```python
from bedspec import Bed3, Bed4
from bedspec.overlap import TreeDetector

bed1 = Bed3("chr1", start=1, end=4)
bed2 = Bed3("chr1", start=5, end=9)

detector = TreeDetector[Bed3]([bed1, bed2])

my_feature = Bed4("chr1", start=2, end=3, name="hi-mom")

assert detector.overlaps(my_feature)
```

The overlap detector supports the following operations:

- `overlapping`: return all overlapping features
- `overlaps`: test if any overlapping features exist
- `enclosed_by`: return those enclosed by the input feature
- `enclosing`: return those enclosing the input feature

A zero-length feature overlaps the features that hold either base beside it.

A BED record is found by any of its `spans()`, so `Bed2` points and `BedPE` pairs are supported, a `BedPE` is found by either end, and a `Bed12` is found by its blocks.
A query that touches only an intron of a `Bed12` does not find it, so to find a `Bed12` by its whole span, add it as a `Bed6`:

```python
intron = Bed3("chr1", start=150, end=160)

assert list(TreeDetector[Bed12]([transcript]).overlapping(intron)) == []

whole = Bed6(
    transcript.refname,
    start=transcript.start,
    end=transcript.end,
    name=transcript.name,
    score=transcript.score,
    strand=transcript.strand,
)

assert list(TreeDetector[Bed6]([whole]).overlapping(intron)) == [whole]
```

Each matching feature is returned once, even when several of its spans match.
A feature encloses the input feature when any one of its spans does.
A feature is enclosed by the input feature only when all of its spans are, so a `BedPE` needs both ends inside.
A query with a `refname`, `start`, and `end` is matched by its whole extent, so a `Bed12` query is matched across its introns, and any other query by its spans, so a `Bed2` query is matched by its base and a `BedPE` query by either end.
A `TabixDetector` answers every query the same way: its index finds each record by its whole extent, and each record found is then tested by its spans, so a query in an intron finds nothing there either.

Each operation takes `stranded=True` to find only features on the same strand as the query.
For a `BedPE`, each end is compared by its own strand, and each block of a `Bed12` by the record's strand.
For the opposite strand, flip the query's strand with `dataclasses.replace`:

```python
from dataclasses import replace

from bedspec import Bed6, BedStrand

plus = Bed6("chr1", start=1, end=4, name=None, score=None, strand=BedStrand.Positive)
minus = Bed6("chr1", start=1, end=4, name=None, score=None, strand=BedStrand.Negative)
stranded = TreeDetector[Bed6]([plus, minus])

opposite = replace(plus, strand=plus.strand.opposite())

assert list(stranded.overlapping(plus, stranded=True)) == [plus]
assert list(stranded.overlapping(opposite, stranded=True)) == [minus]
```

### Other Libraries' Records

Detectors and territories take records from other libraries once `register_locatable()` is given a function that returns each record's 0-based half-open spans, and they return those records as themselves.
In a `with` statement, the registration lasts only for the block.

With the functions in `bedspec.ext.pysam`, keep every pair of a coordinate-sorted BAM that touches a target, tagging both mates with the targets' names:

<!-- skip: next -->
```python
import pysam

from bedspec import Bed4, register_locatable
from bedspec.ext.pysam import as_paired_intervals
from bedspec.overlap import TabixDetector

with (
    register_locatable(pysam.AlignedSegment, as_paired_intervals),
    TabixDetector[Bed4]("targets.bed.gz") as targets,
    pysam.AlignmentFile("sorted.bam") as bam,
    pysam.AlignmentFile("on-target.bam", "wb", template=bam) as kept,
):
    for read in bam:
        if names := {target.name for target in targets.overlapping(read) if target.name}:
            read.set_tag("tg", ",".join(sorted(names)))
            kept.write(read)
```

Find reads by their aligned blocks, so a read spliced right over an exon does not count, and on the read's strand against stranded exons:

<!-- skip: next -->
```python
from bedspec import Bed6
from bedspec.ext.pysam import as_blocks

with (
    register_locatable(pysam.AlignedSegment, as_blocks),
    TabixDetector[Bed6]("exons.bed.gz") as exons,
    pysam.AlignmentFile("rna.bam") as bam,
):
    on_their_exons = [read for read in bam if exons.overlaps(read, stranded=True)]
```

Hold variants in a detector and get the `VariantRecord`s back, with no wrapper class:

<!-- skip: next -->
```python
from bedspec import Bed3
from bedspec.ext.pysam import as_interval
from bedspec.overlap import TreeDetector

with (
    register_locatable(pysam.VariantRecord, as_interval),
    pysam.VariantFile("calls.vcf.gz") as vcf,
):
    variants = TreeDetector(vcf)
    in_exon = list(variants.overlapping(Bed3("chr1", start=1_000, end=1_200)))
```

Build the territory a BAM's reads cover, keeping only the merged bases and not the reads, and ask what it holds:

<!-- skip: next -->
```python
from bedspec import Territory

with (
    register_locatable(pysam.AlignedSegment, as_blocks),
    pysam.AlignmentFile("sorted.bam") as bam,
):
    covered = Territory(bam)

assert Bed3("chr1", start=1_050, end=1_060) in covered
```

Keep whole templates, supplementary alignments included, from a BAM grouped by read name:

<!-- skip: next -->
```python
from fgpyo.sam import Template

from bedspec.ext.fgpyo import as_template_blocks

with (
    register_locatable(Template, as_template_blocks),
    TabixDetector[Bed4]("targets.bed.gz") as targets,
    pysam.AlignmentFile("grouped.bam") as bam,
    pysam.AlignmentFile("kept.bam", "wb", template=bam) as kept,
):
    for template in Template.iterator(bam):
        if targets.overlaps(template):
            template.write_to(kept)
```

Find the targets in the per-target coverage that Picard's `CollectHsMetrics` writes, read with [captainslog](https://github.com/clintval/captainslog), by their 1-based coordinates made 0-based half-open:

<!-- skip: next -->
```python
from captainslog.picard import PerTargetCoverage, PicardCoverageReader

from bedspec.ext.captainslog import as_interval

with (
    register_locatable(PerTargetCoverage, as_interval),
    PicardCoverageReader.from_path[PerTargetCoverage]("hs.per_target_coverage.txt") as rows,
):
    targets = TreeDetector(rows)
    [target] = targets.overlapping(Bed3("chr1", start=1_000, end=1_001))

assert (target.name, target.start) == ("T1", 1_001)
```

Register any other type with a function of your own, which can be a `lambda`; here pybedtools, which already reads a GFF file's coordinates as 0-based half-open, though `bedspec.ext.pybedtools.as_interval` does the same and keeps strand:

<!-- skip: next -->
```python
from pybedtools import BedTool, Interval

from bedspec import Bed3, register_locatable
from bedspec.overlap import TreeDetector

register_locatable(Interval, lambda i: [Bed3(i.chrom, start=i.start, end=i.end)])

genes = TreeDetector(BedTool("genes.gff"))
genes.overlaps(Interval("chr1", 100, 200))
```

The `bedspec.ext` modules have functions for these libraries, each installed with its extra, like `pip install 'bedspec[pysam]'`.
Each module names them the same way: `as_interval` finds a record by its interval, from its first base to its last, `as_blocks` by the aligned blocks within it, `as_paired_intervals` and `as_paired_blocks` add those of its mate, and `as_template_intervals` and `as_template_blocks` those of every record of a template, each on its own strand and never joined into one.

| Module | Function | Finds |
| --- | --- | --- |
| `bedspec.ext.pysam` | `as_interval` | an `AlignedSegment` by its whole alignment, on its strand; a `VariantRecord` from its start to its stop, so a symbolic allele to its `END`; a tabix GTF, GFF3, or BED line, on its strand |
| | `as_blocks` | an `AlignedSegment` by its aligned blocks, on its strand |
| | `as_paired_intervals` | an `AlignedSegment` by its interval and its mate's, from its `MC` tag, each on its strand |
| | `as_paired_blocks` | an `AlignedSegment` by its aligned blocks and its mate's, from its `MC` tag, each on its strand |
| `bedspec.ext.pybedtools` | `as_interval` | an `Interval`, on its strand |
| `bedspec.ext.fgpyo` | `as_interval` | a `SupplementaryAlignment` from an `SA` tag, on its strand |
| | `as_template_intervals` | a `Template` by the intervals of its primary and supplementary records, each on its strand |
| | `as_template_blocks` | a `Template` by the aligned blocks of its primary and supplementary records, each on its strand |
| `bedspec.ext.cyvcf2` | `as_interval` | a `Variant` from its start to its end, so a symbolic allele to its `END` |
| `bedspec.ext.mappy` | `as_interval` | an `Alignment` by its whole alignment, on its strand |
| | `as_blocks` | an `Alignment` by its aligned blocks, on its strand |
| `bedspec.ext.pybwa` | `as_interval` | an `AuxHit` from an `XA` tag, on its strand, converted from the 1-based coordinates it holds |
| `bedspec.ext.biopython` | `as_interval` | a `Bio.Align.Alignment`, like a BED12 line, by its first to last aligned base, on its strand |
| | `as_blocks` | a `Bio.Align.Alignment` by its blocks, on its strand |
| `bedspec.ext.captainslog` | `as_interval` | a Picard or Riker per-target or per-base coverage row, or a `FingerprintingDetailMetrics` SNP, converted from the 1-based coordinates it holds; a `RrbsCpgDetailMetrics` by the C of its CpG on the forward strand |

Every module names its functions the same way, so rename at import to use two together:

<!-- skip: next -->
```python
from bedspec.ext.cyvcf2 import as_interval as variant_interval
from bedspec.ext.pysam import as_blocks as read_blocks
```

### Merging, Intersecting, and Subtracting

A `Territory` is the set of bases held by any number of spans and BED records, kept as the fewest spans that hold them.
Spans that overlap or abut are joined, and a zero-length span holds no bases:

```python
from bedspec import Territory

targets = Territory([
    Bed3("chr1", start=0, end=10),
    Bed3("chr1", start=5, end=20),
    Bed3("chr1", start=20, end=30),
    Bed3("chr2", start=0, end=5),
])

assert list(targets) == [Bed3("chr1", start=0, end=30), Bed3("chr2", start=0, end=5)]
assert targets.length == 35
```

Territories combine like sets of bases, with `|` for the union, `&` for the intersection, and `-` for the difference.
For example, the bases of the targets that a gene holds, less a mask of bases to ignore:

```python
gene = Territory([Bed3("chr1", start=25, end=40)])
mask = Territory([Bed3("chr1", start=27, end=28)])

share = (targets & gene) - mask

assert list(share) == [Bed3("chr1", start=25, end=27), Bed3("chr1", start=28, end=30)]
assert share.length == 4
```

Test one base with `contains`, or whether every base of a span or BED record is held with `in`:

```python
assert share.contains("chr1", 26)
assert not share.contains("chr1", 27)
assert Bed3("chr1", start=25, end=27) in share
```

A BED record adds every one of its spans, so a `BedPE` adds both of its ends and a `Bed12` only its blocks, the same bases as its own `territory()`.
Spans are yielded as `Bed3` records, by start within each reference, with references in the order they were first added.
Territories are immutable, hashable, and equal when they hold the same bases.

### Sorting by Reference

Sort records in the order of a genome's reference sequences with a `ReferenceOrder`, a key function built from their names, or read from a `.fai`, a `.dict`, or a SAM or BAM header with `ReferenceOrder.from_path()`:

```python
from bedspec import Bed3, ReferenceOrder

order = ReferenceOrder(["chr2", "chr1"])
records = [Bed3("chr1", start=2, end=3), Bed3("chr2", start=7, end=8)]

assert sorted(records, key=order) == [Bed3("chr2", start=7, end=8), Bed3("chr1", start=2, end=3)]
```

### Custom BED Types

To create a custom BED record, inherit from the relevant BED-type (`PointBed`, `SimpleBed`, `PairBed`).
Custom BED records must be frozen dataclasses too.
A record of several pieces overrides `spans()` to yield them, and its `territory()` and overlap detection follow.

For example, to create a custom BED3+1 class:

```python
from dataclasses import dataclass

from bedspec import SimpleBed

@dataclass(frozen=True)
class Bed3Plus1(SimpleBed):
    refname: str
    start: int
    end: int
    my_custom_field: float | None
```

You can also inherit and extend a pre-existing BED class:

```python
from dataclasses import dataclass

from bedspec import Bed3

@dataclass(frozen=True)
class Bed3Plus1(Bed3):
    my_custom_field: float | None

record = Bed3Plus1(refname="chr1", start=2, end=3, my_custom_field=0.1)

assert record.my_custom_field == 0.1
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
