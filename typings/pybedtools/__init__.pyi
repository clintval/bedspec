from collections.abc import Iterator

class Interval:
    def __init__(
        self,
        chrom: str,
        start: int,
        end: int,
        name: str = ".",
        score: str = ".",
        strand: str = ".",
        otherfields: list[str] | None = None,
    ) -> None: ...
    @property
    def chrom(self) -> str: ...
    @property
    def start(self) -> int: ...
    @property
    def end(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def strand(self) -> str: ...

class BedTool:
    def __init__(
        self, fn: str | None = None, from_string: bool = False, remote: bool = False
    ) -> None: ...
    def __iter__(self) -> Iterator[Interval]: ...

def create_interval_from_list(fields: list[str]) -> Interval: ...
