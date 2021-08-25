"""
Validation of a single oil record

The actual validation is done in the Oil object (and sub-objects)

This just provides some  convenient wrappers around the validate call
"""

from ..oil import Oil

import logging


# Putting these all here so we can keep track more easily
# from .warnings import WARNINGS
from .errors import ERRORS

logger = logging.getLogger(__name__)


def validate_json(oil_json):
    """
    validate a json-compatible-python record

    An Oil object is returned, if it's possible to do so.

    The "status" field is updated in place, with no other alterations of the record
    """

    if "oil_id" not in oil_json:
        raise ValueError(ERRORS["E010"])
    try:
        oil = Oil.from_py_json(oil_json)
    except TypeError as err:
        if "argument: 'oil_id'" in err.args[0]:
            raise TypeError(ERRORS["E011"].format(oil_json["oil_id"]))
        else:
            raise

    oil.reset_validation()

    return oil


def validate(oil):
    """
    validate an Oil object

    oil.status is updated in place -- this is simply a wrapper around
    Oil.reset_validation() -- probably no longer needed
    """
    oil.reset_validation()
