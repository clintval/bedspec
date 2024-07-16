import dataclasses
from abc import ABC
from abc import abstractmethod
from dataclasses import Field
from dataclasses import fields
from enum import StrEnum
from enum import unique
from typing import Any
from typing import ClassVar
from typing import Iterator
from typing import Protocol
from typing import Type
from typing import cast
from typing import get_type_hints
from typing import runtime_checkable

import bedspec


@runtime_checkable
class DataclassInstance(Protocol):
    """A protocol for objects that are dataclass instances."""

    __dataclass_fields__: ClassVar[dict[str, Field[Any]]]


@unique
class BedStrand(StrEnum):
    """BED strands for forward and reverse orientations."""

    Positive = "+"
    """The positive BED strand."""

    Negative = "-"
    """The negative BED strand."""

    def opposite(self) -> "BedStrand":
        """Return the opposite BED strand."""
        if self is BedStrand.Positive:
            return BedStrand.Negative
        else:
            return BedStrand.Positive


@runtime_checkable
class GenomicSpan(Protocol):
    """A structural protocol for 0-based half-open objects located on a reference sequence."""

    contig: str
    start: int
    end: int


@runtime_checkable
class Named(Protocol):
    """A structural protocol for a named BED type."""

    name: str | None


@runtime_checkable
class Stranded(Protocol):
    """A structural protocol for stranded BED types."""

    strand: BedStrand | None


class BedType(ABC, DataclassInstance):
    """An abstract base class for all types of BED records."""

    def __new__(cls, *args: object, **kwargs: object) -> "BedType":
        if not dataclasses.is_dataclass(cls):
            raise TypeError("You must annotate custom BED class definitions with @dataclass!")
        instance: BedType = cast(BedType, object.__new__(cls))
        return instance

    @classmethod
    def decode(cls, line: str) -> "BedType":
        """Decode a line of text into a BED record."""
        row: list[str] = line.strip().split()
        coerced: dict[str, object] = {}

        try:
            zipped = list(zip(fields(cls), row, strict=True))
        except ValueError as exception:
            raise ValueError(
                f"Expected {len(fields(cls))} fields but found {len(row)} in record:"
                f" '{' '.join(row)}'"
            ) from exception

        hints: dict[str, Type] = get_type_hints(cls)

        for field, value in zipped:
            try:
                if (
                    bedspec._typing.is_optional(hints[field.name])
                    and value == bedspec._io.MISSING_FIELD
                ):
                    coerced[field.name] = None
                else:
                    coerced[field.name] = bedspec._typing.singular_non_optional_type(field.type)(
                        value
                    )
            except ValueError as exception:
                raise TypeError(
                    f"Tried to build the BED field '{field.name}' (of type '{field.type.__name__}')"
                    f" from the value '{value}' but couldn't for record '{' '.join(row)}'"
                ) from exception

        return cls(**coerced)

    @abstractmethod
    def territory(self) -> Iterator[GenomicSpan]:
        """Return intervals that describe the territory of this BED record."""
