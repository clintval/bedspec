from collections.abc import Mapping
from types import MappingProxyType
from typing import Any

from typeline import FieldCodec
from typeline.codecs import delimited

from bedspec._bedspec import BedColor


def _color_from_text(text: str) -> BedColor | None:
    """Read a BED color from its text, where `0` means the record has no color."""
    return None if text == "0" else BedColor.from_string(text)


BED_CODECS: Mapping[Any, FieldCodec[Any]] = MappingProxyType({
    BedColor: FieldCodec(from_text=_color_from_text, into_text=str),
    list[int]: delimited(int),
})
"""How BED fields with their own text formats are read and written, by field type."""
