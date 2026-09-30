from typing import TextIO

from typeline import TsvWriter
from typeline import WriterOptions
from typing_extensions import Unpack
from typing_extensions import override

from bedspec._bedspec import COMMENT_PREFIXES
from bedspec._bedspec import MISSING_FIELD
from bedspec._bedspec import BedType
from bedspec._codecs import BED_CODECS

BED_COMMENT_PREFIXES: tuple[str, ...] = ("#", *sorted(COMMENT_PREFIXES - {"#"}))
"""The BED comment prefixes, with the one added to comment lines that have none first."""


class BedWriter(TsvWriter[BedType]):
    """A writer for writing dataclasses into BED text data."""

    @override
    def __init__(self, handle: TextIO, /, **options: Unpack[WriterOptions]) -> None:
        """Instantiate a new BED writer.

        Args:
            handle: a file-like object to write delimited data to.
            options: the options of the writer, with BED defaults for any not given.
        """
        _ = options.setdefault("none_field", MISSING_FIELD)
        _ = options.setdefault("codecs", BED_CODECS)
        _ = options.setdefault("quoting", False)
        _ = options.setdefault("comment_prefixes", BED_COMMENT_PREFIXES)
        super().__init__(handle, **options)
