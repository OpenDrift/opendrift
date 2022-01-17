#
# Model class definitions for embedded content in our oil records
#
# Note: These are deprecated a bit.  An updated SARAFraction is sitting in
#       models.oil

from enum import Enum


class ProductTypeEnum(Enum):
    crude = 'crude'
    refined = 'refined'


class SaraTypeEnum(Enum):
    saturates = 'Saturates'
    aromatics = 'Aromatics'
    resins = 'Resins'
    asphaltenes = 'Asphaltenes'


class ToxicityTypeEnum(Enum):
    ec = 'EC'  # effective concentration
    lc = 'LC'  # lethal concentration


# for use in the Emulsion related models
class VisualStabilityEnum(Enum):
    entrained = 'Entrained'
    did_not_form = 'Did not form'
    unstable = 'Unstable'
    stable = 'Stable'
    meso_stable = 'Meso-stable'


# for use in the interfacial tension related models
class InterfaceTypeEnum(Enum):
    air = 'air'
    water = 'water'
    seawater = 'seawater'
