from collections.abc import Iterator
from collections.abc import Sequence
from os import PathLike
from typing import IO

import numpy as np
import numpy.typing as npt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Alignment:
    def __init__(
        self,
        sequences: Sequence[SeqRecord | Seq | str],
        coordinates: npt.NDArray[np.int64] | None = None,
    ) -> None: ...
    @property
    def target(self) -> SeqRecord | Seq | str: ...
    @property
    def query(self) -> SeqRecord | Seq | str: ...
    @property
    def coordinates(self) -> npt.NDArray[np.int64]: ...
    @property
    def aligned(self) -> npt.NDArray[np.int64]: ...

def parse(source: str | PathLike[str] | IO[str], fmt: str) -> Iterator[Alignment]: ...
