# bedspec

[![CI](https://github.com/clintval/bedspec/actions/workflows/pythonpackage.yml/badge.svg?branch=main)](https://github.com/clintval/bedspec/actions/workflows/pythonpackage.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12-blue)](https://github.com/clintval/bedspec)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
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
from bedspec import BedWriter, Bed3

bed = Bed3("chr1", start=2, end=8)

with BedWriter(open("test.bed", "w")) as writer:
    writer.write(bed)
```

### Reading

```python
from bedspec import BedReader, Bed3

with BedReader[Bed3](open("test.bed")) as reader:
    for bed in reader:
        print(bed)
```
```console
Bed3(contig="chr1", start=2, start=8)
```

### BED Types

This package provides pre-defined classes for the following BED formats:

```python
from bedspec import Bed2
from bedspec import Bed3
from bedspec import Bed4
from bedspec import Bed5
from bedspec import Bed6
from bedspec import BedPE
```

### Custom BED Types

Creating custom records is as simple as inheriting from the relevent BED-type:

| Type        | Description                                         |
| ---         | ---                                                 |
| `PointBed`  | BED ecords that are a single point (1-length) only. |
| `SimpleBed` | BED ecords that are a single interval.              |
| `PairBed`   | BED ecords that are a pair of intervals.            |

For example, to create a custom BED3+1 class:

```python
from dataclasses import dataclass

from bedspec import SimpleBed

@dataclass
class MyCustomBed(SimpleBed):
    contig: str
    start: int
    end: int
    my_custom_field: float
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
