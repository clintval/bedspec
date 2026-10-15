from collections.abc import Iterable


def joined_blocks(blocks: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return blocks with those that abut, as across an insertion or a mismatch, joined into one."""
    joined: list[tuple[int, int]] = []
    for start, end in blocks:
        if joined and joined[-1][1] == start:
            joined[-1] = (joined[-1][0], end)
        else:
            joined.append((start, end))
    return joined
