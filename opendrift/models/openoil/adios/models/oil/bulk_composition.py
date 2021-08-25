"""
dataclass to store the compounds
"""
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json, JSON_List

from ..common.measurement import MassOrVolumeFraction


@dataclass_to_json
@dataclass
class BulkComposition:
    '''
    object to hold bulk composition -- could be mass or volume fraction
    '''
    name: str = ""
    groups: list = field(default_factory=list)
    method: str = ""
    measurement: MassOrVolumeFraction = None
    comment: str = ""


class BulkCompositionList(JSON_List):
    item_type = BulkComposition
