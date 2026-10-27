"""Tests of the function that reads captainslog's records as spans."""

from pathlib import Path
from typing import Any

import pytest
from captainslog import fgbio
from captainslog import fgumi
from captainslog.fgbio import ErrorRateByReadPositionMetric
from captainslog.fgbio import FgbioReader
from captainslog.fgumi import FgumiReader
from captainslog.picard import FingerprintingDetailMetrics
from captainslog.picard import PerBaseCoverage
from captainslog.picard import PerTargetCoverage
from captainslog.picard import PicardCoverageReader
from captainslog.picard import PicardReader
from captainslog.picard import RrbsCpgDetailMetrics
from captainslog.riker import HybcapMetrics
from captainslog.riker import HybcapPerBase
from captainslog.riker import HybcapPerTarget
from captainslog.riker import RikerReader

import bedspec.ext.captainslog
from bedspec import Bed3
from bedspec import Bed4
from bedspec import register_locatable
from bedspec.ext.captainslog import as_interval
from bedspec.overlap import TreeDetector

PER_TARGET_COVERAGE = """chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage\tnormalized_coverage\tmin_normalized_coverage\tmax_normalized_coverage\tmin_coverage\tmax_coverage\tpct_0x\tread_count
1\t201\t300\t100\trs7555566\t0.41\t77.49\t1.098526\t0.921463\t1.261695\t65\t89\t0\t140
3\t201\t300\t100\trs17272796\t0.38\t76.45\t1.083782\t0.878934\t1.332577\t62\t94\t0\t149
"""  # noqa: E501

PER_BASE_COVERAGE = """chrom\tpos\ttarget\tcoverage
1\t201\trs7555566\t70
1\t300\trs7555566\t82
"""

HYBCAP_PER_TARGET = """sample\tchrom\tstart\tend\tlength\tname\tgc_frac\tmean_coverage\tnormalized_coverage\tmin_coverage\tmax_coverage\tfrac_0x\tread_count
NA12891\t1\t201\t300\t100\trs7555566\t0.41000\t77.49\t1.09853\t65\t89\t0.00000\t140
NA12891\t3\t201\t300\t100\trs17272796\t0.38000\t76.45\t1.08378\t62\t94\t0.00000\t149
"""  # noqa: E501

HYBCAP_PER_BASE = """sample\tchrom\tpos\ttarget\tcoverage
NA12891\t1\t201\trs7555566\t70
NA12891\t1\t300\trs7555566\t82
"""

FINGERPRINTING_DETAIL_METRICS = """## htsjdk.samtools.metrics.StringHeader
# CheckFingerprint INPUT=in/NA12891.over.fingerprints.r1.sam GENOTYPES=in/NA12891.fp.vcf
## htsjdk.samtools.metrics.StringHeader
# Started on: Sat Oct 24 21:43:28 HST 2026

## METRICS CLASS\tpicard.analysis.FingerprintingDetailMetrics
READ_GROUP\tSAMPLE\tSNP\tSNP_ALLELES\tCHROM\tPOSITION\tEXPECTED_GENOTYPE\tOBSERVED_GENOTYPE\tLOD\tOBS_A\tOBS_B
HJYN2CCXX160123.1.TATCCAGG\tNA12891\trs7555566\tAg\t1\t14804874\tAA\tAA\t0.239099\t0\t0
HJYN2CCXX160123.1.TATCCAGG\tNA12891\trs17272796\tCt\t3\t17077268\tCT\tCT\t0.082835\t0\t0


"""  # noqa: E501

RRBS_DETAIL_METRICS = """## htsjdk.samtools.metrics.StringHeader
# CollectRrbsMetrics INPUT=in/input.sam METRICS_FILE_PREFIX=picard/rrbs REFERENCE=in/test.fasta
## htsjdk.samtools.metrics.StringHeader
# Started on: Sat Oct 24 21:43:34 HST 2026

## METRICS CLASS\tpicard.analysis.RrbsCpgDetailMetrics
SEQUENCE_NAME\tPOSITION\tTOTAL_SITES\tCONVERTED_SITES\tPCT_CONVERTED\tSAMPLE\tLIBRARY\tREAD_GROUP
chr1\t20\t4\t0\t0\t\t\t
chr1\t20\t4\t0\t0\tunknown\t\t
chr1\t20\t4\t0\t0\tunknown\tunknown\t
chr1\t20\t4\t0\t0\tunknown\tunknown\tunknown


"""

CHR1 = "TTCATGCTGAAGCCCTCTTACGATCGTACAGATGCAAATATTAACAAACC"
"""The first bases of chr1 of the reference the RRBS rows were collected against."""

REVIEW = """chrom\tpos\tref\tgenotype\tfilters\tA\tC\tG\tT\tN\tconsensus_read\tconsensus_insert\tconsensus_call\tconsensus_qual\ta\tc\tg\tt\tn
chr1\t501\tA\tA/G\tPASS\t2\t0\t1\t0\t0\tS1-lib|S2-lib:4/1\tchr1:452-770 | F1R2\tG\t72\t0\t0\t1\t0\t0
"""  # noqa: E501
"""A consensus read carrying a variant, as fgumi and fgbio both write it."""

CHR1_FROM_490 = "TCATAATCTCAGTGCTGGTG"
"""The bases of chr1 from 0-based 490 of the reference the reviewed reads were aligned to."""


def written(tmp_path: Path, name: str, text: str) -> Path:
    """Return the path of a file in a temporary directory holding text."""
    path = tmp_path / name
    _ = path.write_text(text)
    return path


class Target(PerTargetCoverage):
    """A user's own class of Picard's per-target coverage."""

    @property
    def is_long(self) -> bool:
        """Whether the target is longer than a read."""
        return self.length > 150


def test_as_interval_of_a_target_runs_from_its_1_based_start_to_its_end(tmp_path: Path) -> None:
    """Test that a Picard or Riker target starts a base earlier, 0-based, and keeps its length."""
    with (
        PicardCoverageReader.from_path[PerTargetCoverage](
            written(tmp_path, "s.per_target_coverage.txt", PER_TARGET_COVERAGE)
        ) as picard,
        RikerReader.from_path[HybcapPerTarget](
            written(tmp_path, "s.hybcap-per-target.txt", HYBCAP_PER_TARGET)
        ) as riker,
    ):
        targets: list[PerTargetCoverage | HybcapPerTarget] = [*picard, *riker]
    spans = [span for target in targets for span in as_interval(target)]
    assert spans == 2 * [
        Bed4("1", start=200, end=300, name="rs7555566"),
        Bed4("3", start=200, end=300, name="rs17272796"),
    ]
    assert [span.end - span.start for span in spans] == [target.length for target in targets]


def test_as_interval_of_a_base_covers_its_1_based_position(tmp_path: Path) -> None:
    """Test that a Picard or Riker base covers its position, 0-based, named for its target."""
    with (
        PicardCoverageReader.from_path[PerBaseCoverage](
            written(tmp_path, "s.per_base_coverage.txt", PER_BASE_COVERAGE)
        ) as picard,
        RikerReader.from_path[HybcapPerBase](
            written(tmp_path, "s.hybcap-per-base.txt", HYBCAP_PER_BASE)
        ) as riker,
    ):
        bases: list[PerBaseCoverage | HybcapPerBase] = [*picard, *riker]
    assert [span for base in bases for span in as_interval(base)] == 2 * [
        Bed4("1", start=200, end=201, name="rs7555566"),
        Bed4("1", start=299, end=300, name="rs7555566"),
    ]


def test_as_interval_of_a_fingerprinted_snp_covers_its_1_based_position(tmp_path: Path) -> None:
    """Test that a SNP Picard fingerprinted covers its position, 0-based, named for the SNP."""
    path = written(tmp_path, "s.fingerprinting_detail_metrics", FINGERPRINTING_DETAIL_METRICS)
    with PicardReader.from_path[FingerprintingDetailMetrics](path) as reader:
        assert [span for snp in reader for span in as_interval(snp)] == [
            Bed4("1", start=14804873, end=14804874, name="rs7555566"),
            Bed4("3", start=17077267, end=17077268, name="rs17272796"),
        ]


def test_as_interval_of_a_cpg_covers_its_c_at_the_0_based_position_picard_writes(
    tmp_path: Path,
) -> None:
    """Test that a CpG covers the C of its CpG on the forward strand, at its 0-based position."""
    path = written(tmp_path, "s.rrbs_detail_metrics", RRBS_DETAIL_METRICS)
    with PicardReader.from_path[RrbsCpgDetailMetrics](path) as reader:
        spans = [span for cpg in reader for span in as_interval(cpg)]
    assert spans == 4 * [Bed3("chr1", start=20, end=21)]
    assert CHR1[spans[0].start : spans[0].end + 1] == "CG"


def test_as_interval_of_a_reviewed_variant_covers_its_1_based_position(tmp_path: Path) -> None:
    """Test that a variant fgumi or fgbio reviewed covers its position, 0-based, and is unnamed."""
    path = written(tmp_path, "review.txt", REVIEW)
    with (
        FgbioReader.from_path[fgbio.ConsensusVariantReviewInfo](path) as by_fgbio,
        FgumiReader.from_path[fgumi.ConsensusVariantReviewInfo](path) as by_fgumi,
    ):
        variants: list[fgbio.ConsensusVariantReviewInfo | fgumi.ConsensusVariantReviewInfo] = [
            *by_fgbio,
            *by_fgumi,
        ]
    spans = [span for variant in variants for span in as_interval(variant)]
    assert spans == 2 * [Bed3("chr1", start=500, end=501)]
    assert CHR1_FROM_490[spans[0].start - 490 : spans[0].end - 490] == variants[0].ref


def test_as_interval_of_a_target_or_base_without_a_name_is_unnamed() -> None:
    """Test that a target or a base whose name is empty has an interval without a name."""
    assert as_interval(PerTargetCoverage(chrom="1", start=201, end=300)) == [
        Bed4("1", start=200, end=300, name=None)
    ]
    assert as_interval(HybcapPerBase(chrom="1", pos=201)) == [
        Bed4("1", start=200, end=201, name=None)
    ]


def test_as_interval_of_a_subclass_of_a_target_is_the_interval_of_the_target(
    tmp_path: Path,
) -> None:
    """Test that a row read into a user's subclass of a Picard target keeps its interval."""
    path = written(tmp_path, "s.per_target_coverage.txt", PER_TARGET_COVERAGE)
    with (
        PicardCoverageReader.from_path[PerTargetCoverage](path) as targets,
        PicardCoverageReader.from_path[Target](path) as subclassed,
    ):
        own = [as_interval(target) for target in subclassed]
        assert own == [as_interval(target) for target in targets]
    assert own[0] == [Bed4("1", start=200, end=300, name="rs7555566")]


def test_registered_targets_are_found_by_their_first_and_last_bases_only(tmp_path: Path) -> None:
    """Test that registered targets are found by their first and last bases, not those beyond."""
    with (
        PicardCoverageReader.from_path[PerTargetCoverage](
            written(tmp_path, "s.per_target_coverage.txt", PER_TARGET_COVERAGE)
        ) as reader,
        PicardCoverageReader.from_path[PerBaseCoverage](
            written(tmp_path, "s.per_base_coverage.txt", PER_BASE_COVERAGE)
        ) as bases,
        register_locatable(PerTargetCoverage, as_interval),
        register_locatable(PerBaseCoverage, as_interval),
    ):
        targets = TreeDetector(reader)
        assert [[target.name for target in targets.overlapping(base)] for base in bases] == [
            ["rs7555566"],
            ["rs7555566"],
        ]
        touching = targets.overlapping(Bed3("1", start=150, end=201))
        assert [target.name for target in touching] == ["rs7555566"]
        assert not targets.overlaps(Bed3("1", start=150, end=200))
        assert not targets.overlaps(Bed3("1", start=300, end=350))
        everywhere = targets.overlapping(Bed3("3", start=0, end=1_000))
        assert [target.name for target in everywhere] == ["rs17272796"]


def test_registered_reviewed_variants_are_found_by_their_base_only(tmp_path: Path) -> None:
    """Test that registered reviewed variants are found by their base, not the bases beside it."""
    path = written(tmp_path, "review.txt", REVIEW)
    with (
        FgumiReader.from_path[fgumi.ConsensusVariantReviewInfo](path) as reader,
        register_locatable(fgumi.ConsensusVariantReviewInfo, as_interval),
    ):
        variants = TreeDetector(reader)
        found = variants.overlapping(Bed3("chr1", start=500, end=501))
        assert [variant.consensus_read for variant in found] == ["S1-lib|S2-lib:4/1"]
        assert not variants.overlaps(Bed3("chr1", start=499, end=500))
        assert not variants.overlaps(Bed3("chr1", start=501, end=502))


def test_as_interval_of_a_record_of_another_class_raises_naming_it() -> None:
    """Test that a record of a class without an interval, like Riker's hybcap summary, raises."""
    record: Any = HybcapMetrics()
    with pytest.raises(TypeError, match=r"of type captainslog\.riker\..*HybcapMetrics!"):
        as_interval(record)


def test_as_interval_of_an_error_rate_by_read_position_raises() -> None:
    """Test that fgbio's error rate at a position within a read, not on a reference, raises."""
    record: Any = ErrorRateByReadPositionMetric(position=501)
    with pytest.raises(TypeError, match=r"captainslog\.fgbio\..*ErrorRateByReadPositionMetric!"):
        as_interval(record)


def test_as_interval_of_a_cpg_without_a_position_raises() -> None:
    """Test that a CpG made without a position has no interval, and raises saying so."""
    with pytest.raises(ValueError, match="needs a position"):
        _ = as_interval(RrbsCpgDetailMetrics(sequence_name="chr1"))


def test_bedspec_ext_captainslog_holds_only_the_functions_it_exports() -> None:
    """Test that bedspec.ext.captainslog holds its function and none of the names it imports."""
    module = bedspec.ext.captainslog
    public = sorted(name for name in vars(module) if not name.startswith("_"))
    assert public == sorted(module.__all__) == ["as_interval"]
