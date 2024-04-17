import inspect
import io
import typing
from csv import DictWriter
from dataclasses import asdict as as_dict
from types import FrameType
from types import TracebackType
from typing import Any
from typing import ContextManager
from typing import Generic
from typing import Iterator
from typing import TypeVar
from typing import cast

from bedspec import BedType

T = TypeVar("T", bound=BedType)


class BedWriter(Generic[T], ContextManager):
    bed_kind: type[T] | None

    def __class_getitem__(cls, key: Any) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return typing._GenericAlias(cls, key)  # type: ignore[attr-defined,no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedWriter[T]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        argvalues = inspect.getargvalues(signature)
        typelevel = argvalues.locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __enter__(self) -> "BedWriter[T]":
        """Enter this context."""
        return self

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED writer wihout knowing yet what BED types we will write."""
        self._handle = handle
        self._writer: DictWriter | None = None

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return super().__exit__(__exc_type, __exc_value, __traceback)

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()

    def write(self, bed: T) -> None:
        """Write a BED record to the BED output."""
        if self.bed_kind is not None:
            if type(bed) is not self.bed_kind:
                raise TypeError(
                    f"BedWriter can only continue to write features of the same type."
                    f" Will not write a {type(bed).__name__} after a {self.bed_kind.__name__}."
                )
        else:
            self.bed_kind = type(bed)

        if self._writer is None:
            self._writer = DictWriter(
                self._handle, delimiter="\t", fieldnames=self.bed_kind.fieldnames()
            )

        self._writer.writerow(as_dict(bed))

    def write_comment(self, comment: str) -> None:
        """Write a comment to the BED output."""
        for line in comment.splitlines():
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                self._handle.write(f"{comment}\n")
            else:
                self._handle.write(f"# {comment}\n")


class BedReader(Generic[T], ContextManager):
    bed_kind: type[T] | None

    def __class_getitem__(cls, key: Any) -> type:
        """Wrap all objects of this class to become generic aliases."""
        return typing._GenericAlias(cls, key)  # type: ignore[attr-defined,no-any-return]

    def __new__(cls, handle: io.TextIOWrapper) -> "BedReader[T]":
        """Bind the kind of BED type to this class for later introspection."""
        signature = cast(FrameType, cast(FrameType, inspect.currentframe()).f_back)
        argvalues = inspect.getargvalues(signature)
        typelevel = argvalues.locals.get("self", None)
        bed_kind = None if typelevel is None else typelevel.__args__[0]
        instance = super().__new__(cls)
        instance.bed_kind = bed_kind
        return instance

    def __init__(self, handle: io.TextIOWrapper) -> None:
        """Initialize a BED reader wihout knowing yet what BED types we will write."""
        self._handle = handle

    def __enter__(self) -> "BedReader[T]":
        """Enter this context."""
        return self

    def __iter__(self) -> Iterator[T]:
        """Iterate through the BED records of this IO handle."""
        for line in self._handle:
            if line.strip() == "":
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            yield self._decode(line)

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> bool | None:
        """Close and exit this context."""
        self.close()
        return super().__exit__(__exc_type, __exc_value, __traceback)

    def _decode(self, line: str) -> T:
        if self.bed_kind is None:
            raise NotImplementedError("Untyped reading is not yet supported!")
        return cast(T, self.bed_kind.decode(line))

    def close(self) -> None:
        """Close the underlying IO handle."""
        self._handle.close()
