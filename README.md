# bedspec

[![PyPi Release](https://badge.fury.io/py/bedspec.svg)](https://badge.fury.io/py/bedspec)
[![CI](https://github.com/clintval/bedspec/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/clintval/bedspec/actions/workflows/tests.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12-blue)](https://github.com/clintval/bedspec)
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

### Writing

```python
from bedspec import Bed3, BedWriter

with BedWriter[Bed3].from_path("test.bed") as writer:
    writer.write(Bed3("chr1", start=2, end=8))
```

### Reading

```python
from bedspec import Bed3, BedReader

with BedReader[Bed3].from_path("test.bed") as reader:
    print(list(reader))
```
```console
[Bed3(contig="chr1", start=2, start=8)]
```

### BED Types

This package provides builtin classes for the following BED formats:

```python
from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedGraph
from bedspec import BedPE
```

### Overlap Detection

Use a fast overlap detector for any collection of interval types, including third-party:

```python
from bedspec import Bed3, Bed4
from bedspec.overlap import OverlapDetector

bed1: Bed4 = Bed4(contig="chr1", start=1, end=4, name="bed1")
bed2: Bed4 = Bed4(contig="chr1", start=5, end=9, name="bed2")

detector: OverlapDetector[Bed4] = OverlapDetector([bed1, bed2])

assert detector.get_overlapping(Bed3(contig="chr1", start=2, 3)) == bed1
```

### Create Custom BED Types

To create a custom BED record, inherit from the relevant BED-type (`PointBed`, `SimpleBed`, `PairBed`).

For example, to create a custom BED3+1 class:

```python
from dataclasses import dataclass

from bedspec import SimpleBed

@dataclass(eq=True, frozen=True)
class MyCustomBed(SimpleBed):
    contig: str
    start: int
    end: int
    my_custom_field: float
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
