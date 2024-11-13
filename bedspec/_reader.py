from io import TextIOWrapper
from pathlib import Path
from types import NoneType
from types import UnionType
from typing import Any
from typing import get_args
from typing import get_origin

from typeline import TsvRecordReader
from typing_extensions import Self
from typing_extensions import override

from bedspec._bedspec import COMMENT_PREFIXES
from bedspec._bedspec import MISSING_FIELD
from bedspec._bedspec import BedColor
from bedspec._bedspec import BedStrand
from bedspec._bedspec import BedType


class BedReader(TsvRecordReader[BedType]):
    """A reader of BED records."""

    @override
    def __init__(
        self,
        handle: TextIOWrapper,
        record_type: type[BedType],
        /,
        header: bool = False,
        comment_prefixes: set[str] = COMMENT_PREFIXES,
    ):
        """Instantiate a new BED reader.

        Args:
            handle: a file-like object to read delimited data from.
            record_type: the type of BED record we will be writing.
            header: whether we expect the first line to be a header or not.
            comment_prefixes: skip lines that have any of these string prefixes.
        """
        super().__init__(handle, record_type, header=header, comment_prefixes=comment_prefixes)

    @override
    def _decode(self, field_type: type[Any] | str | Any, item: str) -> str:
        """A callback for overriding the string formatting of builtin and custom types."""
        if field_type is BedStrand:
            return f'"{item}"'
        elif field_type is BedColor:
            color = BedColor.from_string(item)
            return f'{{"r":{color.r},"g":{color.g},"b":{color.b}}}'

        is_union: bool = isinstance(field_type, UnionType)
        type_args: tuple[type, ...] = get_args(field_type)
        is_optional: bool = is_union and NoneType in type_args

        if is_optional:
            if item == MISSING_FIELD:
                return "null"
            elif type_args.count(BedStrand) == 1:
                return f'"{item}"'
            elif type_args.count(BedColor) == 1:
                if item == "0":
                    return "null"
                else:
                    color = BedColor.from_string(item)
                    return f'{{"r":{color.r},"g":{color.g},"b":{color.b}}}'

        type_origin: type | None = get_origin(field_type)

        if type_origin in (frozenset, list, set, tuple):
            return f"[{item.rstrip(',')}]"

        return super()._decode(field_type, item=item)

    @classmethod
    @override
    def from_path(
        cls,
        path: Path | str,
        record_type: type[BedType],
        /,
        header: bool = False,
        comment_prefixes: set[str] = COMMENT_PREFIXES,
    ) -> Self:
        """Construct a BED reader from a file path.

        Args:
            path: the path to the file to read delimited data from.
            record_type: the type of the object we will be writing.
            header: whether we expect the first line to be a header or not.
            comment_prefixes: skip lines that have any of these string prefixes.
        """
        handle = Path(path).open("r")
        reader = cls(handle, record_type, header=header, comment_prefixes=comment_prefixes)
        return reader
