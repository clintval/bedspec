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
A feature may start where it ends, as an insertion does.

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

### BED Types

This package provides builtin classes for the following BED formats:

```pycon
>>> from bedspec import Bed2
>>> from bedspec import Bed3
>>> from bedspec import Bed4
>>> from bedspec import Bed5
>>> from bedspec import Bed6
>>> from bedspec import Bed12
>>> from bedspec import BedGraph
>>> from bedspec import BedPE

```

For BED files with extra columns (BEDn+m), use `Bed3N`, `Bed4N`, `Bed5N`, `Bed6N`, or `Bed12N`.
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

### Overlap Detection

Use a fast overlap detector for any collection of interval types, including third-party:

```pycon
>>> from bedspec import Bed3, Bed4
>>> from bedspec.overlap import OverlapDetector
>>>
>>> bed1 = Bed3("chr1", start=1, end=4)
>>> bed2 = Bed3("chr1", start=5, end=9)
>>> 
>>> detector = OverlapDetector[Bed3]([bed1, bed2])
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

### Custom BED Types

To create a custom BED record, inherit from the relevant BED-type (`PointBed`, `SimpleBed`, `PairBed`).

For example, to create a custom BED3+1 class:

```pycon
>>> from dataclasses import dataclass
>>> 
>>> from bedspec import SimpleBed
>>> 
>>> @dataclass
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
>>> @dataclass
... class Bed3Plus1(Bed3):
...     my_custom_field: float | None
>>>
>>> Bed3Plus1(refname="chr1", start=2, end=3, my_custom_field=0.1)
Bed3Plus1(refname='chr1', start=2, end=3, my_custom_field=0.1)

```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
