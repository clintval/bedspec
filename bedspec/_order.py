import gzip
import struct
import zlib
from collections.abc import Iterable
from collections.abc import Iterator
from pathlib import Path
from typing import Any
from typing import final

from typing_extensions import Self

from bedspec._bedspec import _REGISTRY  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import _UNREGISTERED  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import Locatable
from bedspec._bedspec import _registered_for  # pyright: ignore[reportPrivateUsage]
from bedspec._bedspec import _spans  # pyright: ignore[reportPrivateUsage]

_FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".fna"})
"""The suffixes of a FASTA file."""

_COMPRESSED_SUFFIXES = frozenset({".gz", ".bgz"})
"""The suffixes of a compressed file, after the suffix of what it holds."""


@final
class ReferenceOrder:
    """An order of reference sequences, to sort records by reference, then start, then end.

    Called with a record, it returns the record's sort key, the place of its reference in this
    order, then its start, then its end, so it can be the key of `sorted()` or of a sorter that
    spills to disk. A record is located as detectors locate a query: by the function registered
    for its class, by its whole extent if it has a `refname`, `start` and `end`, or else by its
    first span. A record with no spans sorts after every other record.

    Args:
        names: the names of the reference sequences, in order.

    Raises:
        ValueError: If a name is given twice.
    """

    __slots__ = ("_places",)

    def __init__(self, names: Iterable[str]) -> None:
        places: dict[str, int] = {}
        for name in names:
            if name in places:
                raise ValueError(f"The reference {name!r} is named twice!")
            places[name] = len(places)
        self._places: dict[str, int] = places

    def __call__(self, record: Locatable) -> tuple[int, int, int]:
        """Return a record's sort key: its reference's place in this order, its start and end.

        Raises:
            ValueError: If the record's reference is not in this order.
            TypeError: If the spans of the record cannot be found.
        """
        located: Any = record
        if not (_REGISTRY and type(record) not in _UNREGISTERED and _registered_for(type(record))):
            try:
                refname, start, end = located.refname, located.start, located.end
            except (AttributeError, KeyError):
                pass
            else:
                return self._place(refname), start, end
        first = next(iter(_spans(record)), None)
        if first is None:
            return len(self._places), 0, 0
        return self._place(first.refname), first.start, first.end

    def _place(self, refname: str) -> int:
        """Return the place of a reference in this order."""
        try:
            return self._places[refname]
        except KeyError:
            raise ValueError(f"The reference {refname!r} is not in this order!") from None

    @classmethod
    def from_path(cls, path: Path | str) -> Self:
        """Read the order of references from a FASTA index, a sequence dictionary, or a SAM or BAM.

        A FASTA file is ordered by the index beside it, its path plus `.fai`, or else by the
        sequence dictionary beside it, its path with `.dict` in place of its suffixes.

        Raises:
            ValueError: If the file is of another kind, a FASTA file has neither beside it, or a
                BAM file is not one.
        """
        path = Path(path)
        suffixes = [suffix.lower() for suffix in path.suffixes]
        last = suffixes[-1] if suffixes else ""
        if last == ".fai":
            return cls(_index_names(path))
        if last in {".dict", ".sam"}:
            with path.open() as handle:
                return cls(_header_names(handle))
        if last == ".bam":
            return cls(_bam_names(path))
        compressed = last in _COMPRESSED_SUFFIXES
        content = suffixes[:-1] if compressed else suffixes
        if content and content[-1] in _FASTA_SUFFIXES:
            index = path.with_name(f"{path.name}.fai")
            if index.exists():
                return cls(_index_names(index))
            stem = path.name
            for _ in range(2 if compressed else 1):
                stem = stem.removesuffix(Path(stem).suffix)
            dictionary = path.with_name(f"{stem}.dict")
            if dictionary.exists():
                with dictionary.open() as handle:
                    return cls(_header_names(handle))
            raise ValueError(
                f"Found neither {index.name} nor {dictionary.name} beside {path}: index it with"
                + " samtools faidx first!"
            )
        raise ValueError(
            f"Cannot read the order of references from {path}: give a .fai, .dict, .sam or .bam"
            + " file, or a FASTA file indexed beside it!"
        )


def _index_names(path: Path) -> Iterator[str]:
    """Yield the names of the references of a FASTA index, in order."""
    with path.open() as handle:
        for line in handle:
            if line.strip():
                yield line.split("\t", 1)[0]


def _header_names(lines: Iterable[str]) -> Iterator[str]:
    """Yield the names of the references of a SAM header or sequence dictionary, in order."""
    for line in lines:
        if not line.startswith("@"):
            return
        if line.startswith("@SQ"):
            fields = line.rstrip("\n").split("\t")
            yield next(field[3:] for field in fields if field.startswith("SN:"))


def _bam_names(path: Path) -> list[str]:
    """Return the names of the references of a BAM file's header, in order.

    Raises:
        ValueError: If the file is not a BAM file, its header is corrupt, or it ends within its
            header.
    """
    try:
        with gzip.open(path, "rb") as handle:
            if handle.read(4) != b"BAM\x01":
                raise ValueError(f"{path} is not a BAM file!")
            _ = _read_exactly(handle, _count(handle, path), path)
            names: list[str] = []
            for _reference in range(_count(handle, path)):
                name = _read_exactly(handle, _count(handle, path), path)
                names.append(name.rstrip(b"\x00").decode())
                _ = _int32(handle, path)
    except (gzip.BadGzipFile, EOFError, zlib.error) as error:
        raise ValueError(f"{path} is not a BAM file!") from error
    return names


def _count(handle: gzip.GzipFile, path: Path) -> int:
    """Read a length or count of a BAM file's header.

    Raises:
        ValueError: If it is negative, which a BAM file's header never holds.
    """
    count = _int32(handle, path)
    if count < 0:
        raise ValueError(f"{path} has a corrupt header!")
    return count


def _read_exactly(handle: gzip.GzipFile, size: int, path: Path) -> bytes:
    """Read exactly this many bytes of a BAM file's header.

    Raises:
        ValueError: If the file ends first.
    """
    data = handle.read(size)
    if len(data) != size:
        raise ValueError(f"{path} ends within its header!")
    return data


def _int32(handle: gzip.GzipFile, path: Path) -> int:
    """Read a little-endian 32-bit integer of a BAM file's header."""
    value: int = struct.unpack("<i", _read_exactly(handle, 4, path))[0]
    return value
