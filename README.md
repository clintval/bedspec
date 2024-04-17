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
from bedspec import Bed3, BedWriter

bed = Bed3("chr1", start=3, end=9)

with BedWriter(open("test.bed", "w")) as handle:
    handle.write(bed)
```

### Reading

```python
from bedspec import Bed3, BedReader

with BedReader(open("test.bed")) as handle:
    for bed in handle:
        print(bed)
```
```console
Bed3(contig="chr1", start=3, start=9)
```

### BED Types

This package provides pre-defined classes for the following BED formats:

- `Bed2`
- `Bed3`
- `Bed4`
- `Bed5`
- `Bed6`
- `BedPE`

### Custom BED Types

Creating custom records is as simple as inheriting from the relevent BED-type:

| Type        | Description                                      |
| ---         | ---                                              |
| `PointBed`  | Records that are a single point (1-length) only. |
| `SimpleBed` | Records that are a single interval.              |
| `PairBed`   | Records that are a pair of intervals.            |

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
