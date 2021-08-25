#
# Model class definitions for embedded content in our oil records
#
# Note: These are deprecated a bit.  An updated SARAFraction is sitting in
#       models.oil
from dataclasses import dataclass
from ..common.utilities import dataclass_to_json

from .enum_types import SaraTypeEnum


@dataclass_to_json
@dataclass
class SARAFraction:
    sara_type: SaraTypeEnum

    fraction: float
    ref_temp_k: float = 273.15
    weathering: float = 0.0

    standard_deviation: float = None
    replicates: float = None
    method: str = None

    def __post_init__(self):
        self.sara_type = SaraTypeEnum(self.sara_type)
        self.fraction = float(self.fraction)
        self.ref_temp_k = float(self.ref_temp_k)
        self.weathering = float(self.weathering)

        if self.standard_deviation is not None:
            self.standard_deviation = float(self.standard_deviation)

        if self.replicates is not None:
            self.replicates = float(self.replicates)

        if self.method is not None:
            self.method = str(self.method)


@dataclass_to_json
@dataclass
class SARADensity:
    sara_type: SaraTypeEnum

    kg_m_3: float
    ref_temp_k: float = 273.15
    weathering: float = 0.0

    def __post_init__(self):
        self.sara_type = SaraTypeEnum(self.sara_type)
        self.kg_m_3 = float(self.kg_m_3)
        self.ref_temp_k = float(self.ref_temp_k)
        self.weathering = float(self.weathering)


@dataclass_to_json
@dataclass
class MolecularWeight:
    sara_type: SaraTypeEnum

    g_mol: float
    ref_temp_k: float = 273.15
    weathering: float = 0.0

    def __post_init__(self):
        self.sara_type = SaraTypeEnum(self.sara_type)
        self.g_mol = float(self.g_mol)
        self.ref_temp_k = float(self.ref_temp_k)
        self.weathering = float(self.weathering)

    @property
    def kg_mol(self):
        return self.g_mol / 1000.0

    @kg_mol.setter
    def length(self, value):
        self.g_mol = value * 1000.0
