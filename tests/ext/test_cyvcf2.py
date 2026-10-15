"""Tests of the function that reads cyvcf2's records as spans."""

from pathlib import Path

from cyvcf2 import VCF
from cyvcf2 import Variant

import bedspec.ext.cyvcf2
from bedspec import Bed3
from bedspec import register_locatable
from bedspec.ext.cyvcf2 import variant_span
from bedspec.overlap import TreeDetector

VCF_TEXT = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tACG\tA\t.\t.\t.
chr1\t500\t.\tC\tT\t.\t.\t.
"""


def test_variant_span_covers_the_reference_allele_0_based(tmp_path: Path) -> None:
    """Test that a variant at POS 100 spans its reference allele from 0-based 99."""
    path = tmp_path / "variants.vcf"
    _ = path.write_text(VCF_TEXT)
    deletion, snv = list(VCF(str(path)))
    assert variant_span(deletion) == [Bed3("chr1", start=99, end=102)]
    assert variant_span(snv) == [Bed3("chr1", start=499, end=500)]


def test_registered_variants_are_found_as_themselves(tmp_path: Path) -> None:
    """Test that variants registered by their spans are found and returned as themselves."""
    path = tmp_path / "variants.vcf"
    _ = path.write_text(VCF_TEXT)
    variants = list(VCF(str(path)))
    with register_locatable(Variant, variant_span):
        detector = TreeDetector(variants)
        (found,) = detector.overlapping(Bed3("chr1", start=101, end=102))
        assert found is variants[0]
        assert TreeDetector([Bed3("chr1", start=499, end=500)]).overlaps(variants[1])


def test_bedspec_ext_cyvcf2_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.cyvcf2 holds its functions and none of the names it imports."""
    module = bedspec.ext.cyvcf2
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert public == sorted(module.__all__) == ["variant_span"]
