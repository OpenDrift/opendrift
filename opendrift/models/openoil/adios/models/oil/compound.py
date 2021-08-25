"""
dataclass to store the compounds
"""
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json, JSON_List

from ..common.measurement import MassFraction


@dataclass_to_json
@dataclass
class Compound:
    '''
        Some compounds that will be handled by this dataclass:
        - sulfur_mass_fraction: MassFraction = None
        - carbon_mass_fraction: MassFraction = None
        - hydrogen_mass_fraction: MassFraction = None
        - mercaptan_sulfur_mass_fraction: MassFraction = None
        - nitrogen_mass_fraction: MassFraction = None
        - ccr_percent: MassFraction = None  # conradson carbon residue
        - calcium_mass_fraction: MassFraction = None
        - hydrogen_sulfide_concentration: MassFraction = None
        - salt_content: MassFraction = None
        - paraffin_volume_fraction: MassFraction = None
        - naphthene_volume_fraction: MassFraction = None
        - aromatic_volume_fraction: MassFraction = None
    '''
    name: str = ""
    groups: list = field(default_factory=list)
    method: str = ""
    measurement: MassFraction = None
    comment: str = ""


class CompoundList(JSON_List):
    item_type = Compound
