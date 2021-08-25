#
# Model class definitions for embedded content in our oil records
#
# Note: These are deprecated a bit.  An updated SARAFraction is sitting in
#       models.oil
from dataclasses import dataclass

from ..common.utilities import dataclass_to_json

from ..common.measurement import MassFraction


@dataclass_to_json
@dataclass
class Sara:
    method: str = None

    saturates: MassFraction = None
    aromatics: MassFraction = None
    resins: MassFraction = None
    asphaltenes: MassFraction = None
