import json
from dataclasses import asdict as as_dict
from io import TextIOWrapper
from pathlib import Path
from types import NoneType
from types import UnionType
from typing import Any
from typing import cast
from typing import get_args
from typing import get_origin

from typeline import TsvStructReader
from typing_extensions import Self
from typing_extensions import override

from bedspec._bedspec import MISSING_FIELD
from bedspec._bedspec import BedColor
from bedspec._bedspec import BedStrand
from bedspec._bedspec import BedType


class BedReader(TsvStructReader[BedType]):
    """A reader of BED records."""

    @override
    def __init__(
        self,
        handle: TextIOWrapper,
        record_type: type[BedType],
        /,
        has_header: bool = False,
    ):
        """Initialize the BED reader."""
        super().__init__(handle, record_type, has_header=has_header)

    @property
    @override
    def comment_prefixes(self) -> set[str]:
        return {"#", "browser", "track"}

    @staticmethod
    def _build_union(*types: type) -> type | UnionType:
        if len(types) == 1:
            return types[0]
        union: UnionType | type = types[0]
        for t in types[1:]:
            union |= t
        return cast(UnionType, union)

    @override
    def _decode(self, field_type: type[Any] | str | Any, item: Any) -> Any:
        """A callback for overriding the decoding of builtin types and custom types."""
        type_args: tuple[type, ...] = get_args(field_type)
        type_origin: type | None = get_origin(field_type)
        is_union: bool = isinstance(field_type, UnionType)

        if item == MISSING_FIELD and NoneType in type_args:
            return None
        elif field_type is BedColor or BedColor in type_args:
            if item == "0":
                return None
            return json.dumps(as_dict(BedColor.from_string(cast(str, item))))  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        elif field_type is BedStrand or BedStrand in type_args:
            return f'"{item}"'
        elif type_origin in (frozenset, list, tuple, set):
            stripped: str = item.rstrip(",")
            return f"[{stripped}]"
        elif is_union and len(type_args) >= 2 and NoneType in type_args:
            other_types: set[type] = set(type_args) - {NoneType}
            return self._decode(self._build_union(*other_types), item)
        return super()._decode(field_type, item=item)

    @classmethod
    @override
    def from_path(
        cls,
        path: Path | str,
        record_type: type[BedType],
        /,
        has_header: bool = False,
    ) -> Self:
        """Construct a BED reader from a file path."""
        reader = cls(Path(path).open("r"), record_type, has_header=has_header)
        return reader
