"""
dataclass to store the industry properties

These are tricky, as they can be any weird unit
"""
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json, JSON_List

from ..common.measurement import AnyUnit


@dataclass_to_json
@dataclass
class IndustryProperty:
    '''
    This handles the various odd industry properties.

    For example:

    * Total Acid Number
    * Reid Vapor Pressure
    * Aniline Point
    * Cetane Index
    * Cloud Point
    * Smoke Point
    * Freeze Point
    '''
    name: str = ""
    method: str = ""
    measurement: AnyUnit = None
    comment: str=""

class IndustryPropertyList(JSON_List):
    item_type = IndustryProperty
