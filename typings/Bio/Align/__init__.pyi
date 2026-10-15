from collections.abc import Iterator
from os import PathLike
from typing import IO

import numpy as np
import numpy.typing as npt
from Bio.SeqRecord import SeqRecord

class Alignment:
    @property
    def target(self) -> SeqRecord: ...
    @property
    def query(self) -> SeqRecord: ...
    @property
    def coordinates(self) -> npt.NDArray[np.int64]: ...
    @property
    def aligned(self) -> npt.NDArray[np.int64]: ...

def parse(source: str | PathLike[str] | IO[str], fmt: str) -> Iterator[Alignment]: ...
