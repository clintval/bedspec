from typing import Any

from Cython.Build import build_ext
from setuptools import Extension

ext_modules = [
    Extension(
        "cgranges",
        sources=["cgranges/python/cgranges.pyx", "cgranges/cgranges.c"],
        depends=[
            "cgranges/cgranges.h",
            "cgranges/khash.h",
            "cgranges/python/cgranges.pyx"
        ],
        include_dirs=["cgranges"],
    ),
    Extension(
        "bgzip.bgzip_utils",
        sources=["xbgzip/bgzip_utils/bgzip_utils.pyx"],
        depends=["xbgzip/bgzip_utils/bgzip_utils.pyx"],
        include_dirs=["xbgzip"],
    ),
]

def build(setup_kwargs: dict[str, Any]) -> None:
    """This function is mandatory in order to build the extensions."""
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": build_ext},
        }
    )
