from dataclasses import dataclass
from typing import ClassVar

from shipslog._binding import FromPath

@dataclass(frozen=True, kw_only=True)
class PerTargetCoverage:
    chrom: str = ""
    start: int = 0
    end: int = 0
    length: int = 0
    name: str = ""

@dataclass(frozen=True, kw_only=True)
class PerBaseCoverage:
    chrom: str = ""
    pos: int = 0
    target: str = ""

@dataclass(frozen=True, kw_only=True)
class FingerprintingDetailMetrics:
    snp: str = ""
    chrom: str = ""
    position: int = 0

@dataclass(frozen=True, kw_only=True)
class RrbsCpgDetailMetrics:
    sequence_name: str = ""
    position: int | None = None

class PicardCoverageReader:
    from_path: ClassVar[FromPath]

class PicardReader:
    from_path: ClassVar[FromPath]
