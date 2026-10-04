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

## Quickstart

### Building a BED Feature

```pycon
>>> from bedspec import Bed3
>>> 
>>> bed = Bed3("chr1", start=2, end=8)

```

Records are checked against the BED spec when they are built.
BED has no quoting, so text is written and read as it is, and a value holding a tab is refused.
A feature may start where it ends, as an insertion does.

Records are immutable and hashable.
Use `dataclasses.replace` to build a changed copy, which is checked like any other record:

```pycon
>>> from dataclasses import replace
>>>
>>> replace(bed, end=10)
Bed3(refname='chr1', start=2, end=10)

```

### Writing

```pycon
>>> from bedspec import BedWriter
>>> from tempfile import NamedTemporaryFile
>>> 
>>> temp_file = NamedTemporaryFile(mode="w+t", suffix=".txt")
>>>
>>> with BedWriter.from_path[Bed3](temp_file.name) as writer:
...     writer.write(bed)

```

### Reading

```pycon
>>> from bedspec import BedReader
>>> 
>>> with BedReader.from_path[Bed3](temp_file.name) as reader:
...     for bed in reader:
...         print(bed)
Bed3(refname='chr1', start=2, end=8)

```

### Compressed and Indexed BED

A path ending in `.gz`, `.bgz`, or `.bgzf` is written as BGZF, which any gzip reader can read, and a compressed file is read by its contents.
Ask for a tabix or CSI index to have one written beside the file, with features sorted by reference and start.

```pycon
>>> from pybgzf import IndexFormat
>>>
>>> with BedWriter.from_path[Bed3](f"{temp_file.name}.gz", index=IndexFormat.TBI, threads=4) as writer:
...     writer.write(Bed3("chr1", start=2, end=8))
...     writer.write(Bed3("chr1", start=6, end=9))
>>>
>>> with BedReader.from_path[Bed3](f"{temp_file.name}.gz") as reader:
...     print(list(reader))
[Bed3(refname='chr1', start=2, end=8), Bed3(refname='chr1', start=6, end=9)]

```

Query an indexed file on disk with the same operations as the [overlap detector](#overlap-detection):

```pycon
>>> from bedspec.overlap import TabixDetector
>>>
>>> with TabixDetector[Bed3](f"{temp_file.name}.gz") as detector:
...     print(list(detector.enclosing(Bed3("chr1", start=7, end=8))))
[Bed3(refname='chr1', start=2, end=8), Bed3(refname='chr1', start=6, end=9)]

```

No index query returns a zero-length feature at the start of a reference, so writing one to an indexed file warns.

### BED Types

This package provides builtin classes for the following BED formats:

```pycon
>>> from bedspec import Bed2
>>> from bedspec import Bed3
>>> from bedspec import Bed4
>>> from bedspec import Bed5
>>> from bedspec import Bed6
>>> from bedspec import Bed9
>>> from bedspec import Bed12
>>> from bedspec import BedGraph
>>> from bedspec import BedPE

```

It also provides the ENCODE peak formats:

```pycon
>>> from bedspec import BroadPeak
>>> from bedspec import GappedPeak
>>> from bedspec import NarrowPeak

```

ENCODE writes -1 for a p-value, q-value, or summit that is not given, and so do these types.

For BED files with extra columns (BEDn+m), use `Bed3N`, `Bed4N`, `Bed5N`, `Bed6N`, `Bed9N`, or `Bed12N`.
Each is its BED type plus an `extra` field that keeps any further columns as text.

```pycon
>>> from bedspec import Bed6N
>>>
>>> _ = open(temp_file.name, "w").write("chr1\t5\t9\tpeak\t7\t-\t3.2\t0.01\n")
>>>
>>> with BedReader.from_path[Bed6N](temp_file.name) as reader:
...     for bed in reader:
...         print(bed.name, bed.extra)
peak ('3.2', '0.01')

```

### Spans and Territory

Every BED record yields its pieces on reference sequences with `spans()`, each with its strand.
A `Bed2` yields its one base, a `Bed3` to `Bed9` yields itself, a `BedPE` yields both of its ends, and a `Bed12` yields its blocks, so its introns are not part of it:

```pycon
>>> from bedspec import BedStrand
>>>
>>> transcript = Bed12(
...     "chr1",
...     start=100,
...     end=200,
...     name="tx",
...     score=None,
...     strand=BedStrand.Positive,
...     thick_start=None,
...     thick_end=None,
...     item_rgb=None,
...     block_count=2,
...     block_sizes=(10, 20),
...     block_starts=(0, 80),
... )
>>> for span in transcript.spans():
...     print(span.refname, span.start, span.end, span.strand)
chr1 100 110 +
chr1 180 200 +

```

`territory()` returns the bases those spans hold, without their strands, as a [`Territory`](#merging-intersecting-and-subtracting):

```pycon
>>> transcript.territory()
Territory([Bed3(refname='chr1', start=100, end=110), Bed3(refname='chr1', start=180, end=200)])
>>> transcript.territory().length
30

```

### Overlap Detection

Use a fast overlap detector for any collection of interval types, including third-party:

```pycon
>>> from bedspec import Bed3, Bed4
>>> from bedspec.overlap import TreeDetector
>>>
>>> bed1 = Bed3("chr1", start=1, end=4)
>>> bed2 = Bed3("chr1", start=5, end=9)
>>> 
>>> detector = TreeDetector[Bed3]([bed1, bed2])
>>> 
>>> my_feature = Bed4("chr1", start=2, end=3, name="hi-mom")
>>> detector.overlaps(my_feature)
True

```

The overlap detector supports the following operations:

- `overlapping`: return all overlapping features
- `overlaps`: test if any overlapping features exist
- `enclosed_by`: return those enclosed by the input feature
- `enclosing`: return those enclosing the input feature

A zero-length feature overlaps the features that hold either base beside it.

A BED record is found by any of its `spans()`, so `Bed2` points and `BedPE` pairs are supported, a `BedPE` is found by either end, and a `Bed12` is found by its blocks.
A query that touches only an intron of a `Bed12` does not find it, so to find a `Bed12` by its whole span, add it as a `Bed6`:

```pycon
>>> intron = Bed3("chr1", start=150, end=160)
>>>
>>> list(TreeDetector[Bed12]([transcript]).overlapping(intron))
[]
>>> whole = Bed6(
...     transcript.refname,
...     start=transcript.start,
...     end=transcript.end,
...     name=transcript.name,
...     score=transcript.score,
...     strand=transcript.strand,
... )
>>> list(TreeDetector[Bed6]([whole]).overlapping(intron)) == [whole]
True

```

Each matching feature is returned once, even when several of its spans match.
A feature encloses the input feature when any one of its spans does.
A feature is enclosed by the input feature only when all of its spans are, so a `BedPE` needs both ends inside.
Queries must be spans with an `end`, so a `Bed2` can be added but cannot be used as a query.
A `TabixDetector` answers every query the same way: its index finds each record by its whole extent, and each record found is then tested by its spans, so a query in an intron finds nothing there either.

Each operation takes `stranded=True` to find only features on the same strand as the query.
For a `BedPE`, each end is compared by its own strand, and each block of a `Bed12` by the record's strand.
For the opposite strand, flip the query's strand with `dataclasses.replace`:

```pycon
>>> from dataclasses import replace
>>> from bedspec import Bed6, BedStrand
>>>
>>> plus = Bed6("chr1", start=1, end=4, name=None, score=None, strand=BedStrand.Positive)
>>> minus = Bed6("chr1", start=1, end=4, name=None, score=None, strand=BedStrand.Negative)
>>> stranded = TreeDetector[Bed6]([plus, minus])
>>>
>>> list(stranded.overlapping(plus, stranded=True)) == [plus]
True
>>> list(stranded.overlapping(replace(plus, strand=plus.strand.opposite()), stranded=True)) == [minus]
True

```

### Merging, Intersecting, and Subtracting

A `Territory` is the set of bases held by any number of spans and BED records, kept as the fewest spans that hold them.
Spans that overlap or abut are joined, and a zero-length span holds no bases:

```pycon
>>> from bedspec import Territory
>>>
>>> targets = Territory([
...     Bed3("chr1", start=0, end=10),
...     Bed3("chr1", start=5, end=20),
...     Bed3("chr1", start=20, end=30),
...     Bed3("chr2", start=0, end=5),
... ])
>>> list(targets)
[Bed3(refname='chr1', start=0, end=30), Bed3(refname='chr2', start=0, end=5)]
>>> targets.length
35

```

Territories combine like sets of bases, with `|` for the union, `&` for the intersection, and `-` for the difference.
For example, the bases of the targets that a gene holds, less a mask of bases to ignore:

```pycon
>>> gene = Territory([Bed3("chr1", start=25, end=40)])
>>> mask = Territory([Bed3("chr1", start=27, end=28)])
>>>
>>> share = (targets & gene) - mask
>>> share
Territory([Bed3(refname='chr1', start=25, end=27), Bed3(refname='chr1', start=28, end=30)])
>>> share.length
4

```

Test one base with `contains`, or whether every base of a span or BED record is held with `in`:

```pycon
>>> share.contains("chr1", 26), share.contains("chr1", 27)
(True, False)
>>> Bed3("chr1", start=25, end=27) in share
True

```

A BED record adds every one of its spans, so a `BedPE` adds both of its ends and a `Bed12` only its blocks, the same bases as its own `territory()`.
Spans are yielded as `Bed3` records, by start within each reference, with references in the order they were first added.
Territories are immutable, hashable, and equal when they hold the same bases.

### Custom BED Types

To create a custom BED record, inherit from the relevant BED-type (`PointBed`, `SimpleBed`, `PairBed`).
Custom BED records must be frozen dataclasses too.
A record of several pieces overrides `spans()` to yield them, and its `territory()` and overlap detection follow.

For example, to create a custom BED3+1 class:

```pycon
>>> from dataclasses import dataclass
>>> 
>>> from bedspec import SimpleBed
>>> 
>>> @dataclass(frozen=True)
... class Bed3Plus1(SimpleBed):
...     refname: str
...     start: int
...     end: int
...     my_custom_field: float | None

```

You can also inherit and extend a pre-existing BED class:

```pycon
>>> from dataclasses import dataclass
>>>
>>> from bedspec import Bed3
>>>
>>> @dataclass(frozen=True)
... class Bed3Plus1(Bed3):
...     my_custom_field: float | None
>>>
>>> Bed3Plus1(refname="chr1", start=2, end=3, my_custom_field=0.1)
Bed3Plus1(refname='chr1', start=2, end=3, my_custom_field=0.1)

```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
