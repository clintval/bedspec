from sybil import Sybil
from sybil.parsers.markdown import PythonCodeBlockParser
from sybil.parsers.markdown import SkipParser

pytest_collect_file = Sybil(
    parsers=[PythonCodeBlockParser(), SkipParser()], patterns=["*.md"]
).pytest()
