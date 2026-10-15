from Bio.Seq import Seq

class SeqRecord:
    id: str
    def __init__(self, seq: Seq | str, id: str = "<unknown id>") -> None: ...
