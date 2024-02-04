from dataclasses import dataclass
from types import TracebackType

from typing import ContextManager
from typing import Iterable
from typing import Protocol
from typing import TypeVar


@dataclass
class Interval(Protocol):
    """A 0-based half-closed genomic interval.

    Attributes:
        contig: The reference sequence name of this interval.
        start: The 0-based start position of this interval (BED format).
        end: The half-open end position of this interval (BED format).

    """

    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        """Validate that the the coordinates of this interval are sensible."""
        if self.start < 0 or self.start >= self.end:
            raise ValueError("Start cannot be <0! End must be greater than start!")


BedType = TypeVar("BedType", bound=Interval)


class Bed3(BedType):
    num_static_fields: int = 3


class Bed4(BedType):
    num_static_fields: int = 4


class Bed5(BedType):
    num_static_fields: int = 5


class Bed6(BedType):
    num_static_fields: int = 6


class Bed8(BedType):
    num_static_fields: int = 8


class Bed9(BedType):
    num_static_fields: int = 9


class Bed12(BedType):
    num_static_fields: int = 12


class BedReader[BedType](ContextManager, Iterable[BedType]):

    def __enter__(self) -> "BedReader":
        return super().__enter__()

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        return super().__exit__(__exc_type, __exc_value, __traceback)
