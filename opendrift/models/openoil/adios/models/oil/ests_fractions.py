"""
dataclass to hold the ESTS Hydrocarbon Fractions

This is kind of a special case, so this nails it down
"""
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json

from .compound import CompoundList


@dataclass_to_json
@dataclass
class ESTSFractions:
    """
    hold the ESTS Hydrocarbon Fractions from env CA

    (https://www.ccme.ca/en/resources/canadian_environmental_quality_guidelines/index.html)
    """
    method: str = None

    saturates: CompoundList = field(default_factory=CompoundList)
    aromatics: CompoundList = field(default_factory=CompoundList)
    GC_TPH: CompoundList = field(default_factory=CompoundList)
