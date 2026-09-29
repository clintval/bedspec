from io import TextIOWrapper

from typeline import ReaderOptions
from typeline import TsvReader
from typing_extensions import Unpack
from typing_extensions import override

from bedspec._bedspec import COMMENT_PREFIXES
from bedspec._bedspec import MISSING_FIELD
from bedspec._bedspec import BedType
from bedspec._codecs import BED_CODECS


class BedReader(TsvReader[BedType]):
    """A reader of BED records."""

    @override
    def __init__(self, handle: TextIOWrapper, /, **options: Unpack[ReaderOptions]) -> None:
        """Instantiate a new BED reader.

        Args:
            handle: a file-like object to read delimited data from.
            options: the options of the reader, with BED defaults for any not given.
        """
        _ = options.setdefault("header", False)
        _ = options.setdefault("comment_prefixes", COMMENT_PREFIXES)
        _ = options.setdefault("none_field", MISSING_FIELD)
        _ = options.setdefault("codecs", BED_CODECS)
        super().__init__(handle, **options)
