"""
simple class to hold the review status of a record
"""

from datetime import datetime
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json

from ..common.validators import EnumValidator

from .validation.errors import ERRORS
from .validation.warnings import WARNINGS


@dataclass_to_json
@dataclass
class ReviewStatus:
    status: str = "Not Reviewed"
    reviewers: str = ""
    review_date: str = ""
    notes: str = ""

    _status_validator = EnumValidator(["Not Reviewed",
                                       "Under Review",
                                       "Review Complete"],
                                      ERRORS['E013'],
                                      case_insensitive=True)

    def validate(self):
        msgs = []
        if self.review_date:
            try:
                datetime.fromisoformat(self.review_date)
            except ValueError as err:
                msgs.append(WARNINGS["W011"].format("review date", self.review_date, str(err)))

        msgs.extend(self._status_validator(self.status))

        return msgs

