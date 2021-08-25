"""
dataclass to hold the CCME data

CCME is kind of a special case, so this nails it down
"""
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json

from ..common.measurement import MassFraction


@dataclass_to_json
@dataclass
class CCME:
    """
    hold the CCME data from env CA

    (https://www.ccme.ca/en/resources/canadian_environmental_quality_guidelines/index.html)
    """
    F1: MassFraction = None
    F2: MassFraction = None
    F3: MassFraction = None
    F4: MassFraction = None
    method: str = None


