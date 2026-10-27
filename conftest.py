from collections.abc import Iterator
from importlib.util import find_spec

import pytest
from sybil import Sybil
from sybil.parsers.markdown import PythonCodeBlockParser
from sybil.parsers.markdown import SkipParser

pytest_collect_file = Sybil(
    parsers=[PythonCodeBlockParser(), SkipParser()], patterns=["*.md"]
).pytest()

# captainslog is not on PyPI, so its module's doctests run only where it is installed.
collect_ignore: list[str] = (
    []
    if find_spec("captainslog")
    else ["bedspec/ext/_captainslog.py", "bedspec/ext/captainslog.py"]
)


@pytest.fixture(autouse=True)
def forget_registrations() -> Iterator[None]:
    """Forget every registration a test made once it ends, so that tests stay independent."""
    from bedspec import _bedspec  # pyright: ignore[reportPrivateUsage]

    yield
    _bedspec._CHAINS.clear()  # pyright: ignore[reportPrivateUsage]
    _bedspec._REGISTRY.clear()  # pyright: ignore[reportPrivateUsage]
    _bedspec._forget_resolutions()  # pyright: ignore[reportPrivateUsage]
