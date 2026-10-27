from typing import overload

try:
    from captainslog import fgbio
    from captainslog import fgumi
    from captainslog.picard import FingerprintingDetailMetrics
    from captainslog.picard import PerBaseCoverage
    from captainslog.picard import PerTargetCoverage
    from captainslog.picard import RrbsCpgDetailMetrics
    from captainslog.riker import HybcapPerBase
    from captainslog.riker import HybcapPerTarget
except ImportError as error:
    raise ImportError(
        "bedspec.ext.captainslog needs captainslog: pip install 'bedspec[captainslog]'"
    ) from error

from bedspec import Bed3
from bedspec import Bed4


# fmt: off
@overload
def as_interval(record: PerTargetCoverage) -> list[Bed4]: ...
@overload
def as_interval(record: HybcapPerTarget) -> list[Bed4]: ...
@overload
def as_interval(record: PerBaseCoverage) -> list[Bed4]: ...
@overload
def as_interval(record: HybcapPerBase) -> list[Bed4]: ...
@overload
def as_interval(record: FingerprintingDetailMetrics) -> list[Bed4]: ...
@overload
def as_interval(record: fgbio.ConsensusVariantReviewInfo) -> list[Bed3]: ...
@overload
def as_interval(record: fgumi.ConsensusVariantReviewInfo) -> list[Bed3]: ...
@overload
def as_interval(record: RrbsCpgDetailMetrics) -> list[Bed3]: ...
# fmt: on
def as_interval(record: object) -> list[Bed4] | list[Bed3]:
    """Return the interval of a target, a base, a SNP, a variant, or a CpG, from its coordinates.

    A target runs from its 1-based start to its end and a base or a SNP covers its 1-based position,
    each named for its target or SNP, while a CpG covers its C on the forward strand, at the 0-based
    position Picard writes. A variant fgumi or fgbio reviewed covers its 1-based position, unnamed,
    as its row names a consensus read and not the variant.

    Raises:
        TypeError: If the record is none of captainslog's targets, bases, SNPs, variants, or CpGs.
        ValueError: If the record is a CpG without a position.
    """
    if isinstance(record, PerTargetCoverage | HybcapPerTarget):
        name = record.name or None
        return [Bed4(record.chrom, start=record.start - 1, end=record.end, name=name)]
    if isinstance(record, PerBaseCoverage | HybcapPerBase):
        name = record.target or None
        return [Bed4(record.chrom, start=record.pos - 1, end=record.pos, name=name)]
    if isinstance(record, FingerprintingDetailMetrics):
        name = record.snp or None
        return [Bed4(record.chrom, start=record.position - 1, end=record.position, name=name)]
    if isinstance(record, fgbio.ConsensusVariantReviewInfo | fgumi.ConsensusVariantReviewInfo):
        return [Bed3(record.chrom, start=record.pos - 1, end=record.pos)]
    if isinstance(record, RrbsCpgDetailMetrics):
        if record.position is None:
            raise ValueError("A CpG needs a position to find its interval, not None!")
        return [Bed3(record.sequence_name, start=record.position, end=record.position + 1)]
    kind = type(record)
    raise TypeError(
        f"Cannot find the interval of a record of type {kind.__module__}.{kind.__qualname__}!"
    )
