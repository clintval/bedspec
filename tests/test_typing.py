from bedspec._typing import is_union


def test_is_union() -> None:
    """Test that a union type is the conjunction of two types."""
    # TODO: have a positive unit test for is_union
    assert not is_union(type(int))
    assert not is_union(type(None))
