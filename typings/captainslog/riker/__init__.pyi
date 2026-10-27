from dataclasses import dataclass
from typing import ClassVar

from captainslog._binding import FromPath

@dataclass(frozen=True, kw_only=True)
class HybcapMetrics: ...

@dataclass(frozen=True, kw_only=True)
class HybcapPerTarget:
    chrom: str = ""
    start: int = 0
    end: int = 0
    length: int = 0
    name: str = ""

@dataclass(frozen=True, kw_only=True)
class HybcapPerBase:
    chrom: str = ""
    pos: int = 0
    target: str = ""

class RikerReader:
    from_path: ClassVar[FromPath]
