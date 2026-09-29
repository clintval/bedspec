from typeline import Codecs
from typeline import FieldCodec
from typeline.codecs import delimited
from typeline.codecs import nullable

from bedspec._bedspec import BedColor

BED_CODECS: Codecs = {
    BedColor: nullable(FieldCodec(from_text=BedColor.from_string, into_text=str), missing="0"),
    tuple[int, ...]: delimited(int, container=tuple),
}
"""How BED fields with their own text formats are read and written, by field type."""
