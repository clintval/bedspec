"""Tests that each module of bedspec.ext names its extra when its library is missing."""

import importlib
import sys

import pytest

BLOCKED: dict[str, tuple[str, ...]] = {
    "biopython": ("Bio", "Bio.Align"),
    "cyvcf2": ("cyvcf2",),
    "fgpyo": ("fgpyo", "fgpyo.sam"),
    "mappy": ("mappy",),
    "pybedtools": ("pybedtools",),
    "pybwa": ("pybwa",),
    "pysam": ("pysam", "pysam.libctabixproxies"),
}
"""The modules each extra installs, blocked so importing them fails as if they were missing."""


@pytest.mark.parametrize(("extra", "modules"), BLOCKED.items())
def test_a_module_of_bedspec_ext_names_its_extra_when_its_library_is_missing(
    extra: str, modules: tuple[str, ...], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Test that importing a module of bedspec.ext without its library names its extra."""
    for module in modules:
        monkeypatch.setitem(sys.modules, module, None)
    monkeypatch.delitem(sys.modules, f"bedspec.ext.{extra}", raising=False)
    with pytest.raises(ImportError, match=rf"pip install 'bedspec\[{extra}\]'"):
        _ = importlib.import_module(f"bedspec.ext.{extra}")
