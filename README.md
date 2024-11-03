# bedspec

[![PyPi Release](https://badge.fury.io/py/bedspec.svg)](https://badge.fury.io/py/bedspec)
[![CI](https://github.com/clintval/bedspec/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/bedspec/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.11_|_3.12-blue)](https://github.com/clintval/bedspec)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

An HTS-specs compliant BED toolkit.

## Installation

The package can be installed with `pip`:

```console
pip install bedspec
```

## Quickstart

### Building a BED Feature

```python
from bedspec import Bed3

bed = Bed3("chr1", start=2, end=8)
```

### Writing

```python
from bedspec import BedWriter

with BedWriter.from_path("test.bed") as writer:
    writer.write(bed)
```

### Reading

```python
from bedspec import BedReader

with BedReader.from_path("test.bed", Bed3) as reader:
    for bed in reader:
        print(bed)
```
```console
Bed3(refname="chr1", start=2, end=8)
```

### BED Types

This package provides builtin classes for the following BED formats:

```python
from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import Bed12
from bedspec import BedGraph
from bedspec import BedPE
```

### Overlap Detection

Use a fast overlap detector for any collection of interval types, including third-party:

```python
from bedspec import Bed3, Bed4
from bedspec.overlap import OverlapDetector

bed1 = Bed3("chr1", start=1, end=4)
bed2 = Bed3("chr1", start=5, end=9)

detector = OverlapDetector[Bed3]([bed1, bed2])

my_feature = Bed4("chr1", start=2, end=3, name="hi-mom")

assert detector.overlaps(my_feature) is True
```

The overlap detector supports the following operations:

- `overlapping`: return all overlapping features
- `overlaps`: test if any overlapping features exist
- `enclosed_by`: return those enclosed by the input feature
- `enclosing`: return those enclosing the input feature

### Custom BED Types

To create a custom BED record, inherit from the relevant BED-type (`PointBed`, `SimpleBed`, `PairBed`).

For example, to create a custom BED3+1 class:

```python
from dataclasses import dataclass

from bedspec import SimpleBed

@dataclass(eq=True)
class Bed3Plus1(SimpleBed):
    refname: str
    start: int
    end: int
    my_custom_field: float | None
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
