
"""

Classes for individual classes -not much here!

"""

from ..common.utilities import dataclass_to_json
from ..common.validators import YearValidator

from .validation.warnings import WARNINGS
from .validation.errors import ERRORS

from dataclasses import dataclass


@dataclass_to_json
@dataclass
class ReferenceBase:
    year: int = None
    reference: str = ""


class Reference(ReferenceBase):
    _validator = YearValidator(1600,
                               2050,
                               ERRORS["E012"],
                               )

    @classmethod
    def validate(cls, value):
        year = value.year
        if not year:
            return [WARNINGS["W008"]]
        return cls._validator(year)

