# bedspec

[![CI](https://github.com/clintval/bedspec/actions/workflows/pythonpackage.yml/badge.svg?branch=main)](https://github.com/clintval/bedspec/actions/workflows/pythonpackage.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12-blue)](https://github.com/clintval/bedspec)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

An HTS-specs compliant BED toolkit.

## Installation

`bedspec` may be installed with `pip`:

```console
pip install bedspec
```

## Quickstart

### Writing BED

```python
from bedspec import Bed3
from bedspec.io import BedWriter

with BedWriter(open("test.bed", "w")) as handle:
    handle.write_comment("browser position chr7:127471196-127495720")
    handle.write(Bed3(contig="chr7", start=127_471_196, start=127_495_720))
```


### Reading BED

```python
from bedspec import Bed3
from bedspec.io import BedReader

with BedReader(open("test.bed")) as handle:
    for bed in handle:
        print(bed)
```
```console
Bed3(contig="chr7", start=127_471_196, start=127_495_720)
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
