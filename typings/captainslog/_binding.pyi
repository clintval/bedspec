from collections.abc import Callable
from collections.abc import Iterator
from pathlib import Path
from typing import Generic
from typing import TypeVar

from typing_extensions import Self

_RecordType = TypeVar("_RecordType")

class _Reader(Generic[_RecordType]):
    def __enter__(self) -> Self: ...
    def __exit__(self, *args: object) -> None: ...
    def __iter__(self) -> Iterator[_RecordType]: ...

class FromPath:
    def __getitem__(
        self, record_type: type[_RecordType], /
    ) -> Callable[[Path | str], _Reader[_RecordType]]: ...
