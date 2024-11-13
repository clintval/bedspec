from typing import Any

from typeline import TsvRecordWriter
from typing_extensions import override

from bedspec._bedspec import COMMENT_PREFIXES
from bedspec._bedspec import BedColor
from bedspec._bedspec import BedType


class BedWriter(TsvRecordWriter[BedType]):
    """A writer for writing dataclasses into BED text data."""

    @override
    def _encode(self, item: Any) -> Any:
        """A callback for overriding the encoding of builtin types and custom types."""
        if item is None:
            return "."
        elif isinstance(item, (frozenset, list, set, tuple)):
            return ",".join(map(str, item))  # pyright: ignore[reportUnknownArgumentType]
        elif isinstance(item, BedColor):
            return str(item)
        return super()._encode(item=item)

    def write_comment(self, comment: str) -> None:
        """Write a comment to the BED output."""
        for line in comment.splitlines():
            prefix = "" if any(line.startswith(prefix) for prefix in COMMENT_PREFIXES) else "# "
            _ = self._handle.write(f"{prefix}{line}\n")
