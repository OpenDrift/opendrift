'''
    class that represents the demographic data (metadata) of an oil record.
'''
from datetime import datetime
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json
from ..common.measurement import MassFraction, Temperature

from .values import Reference
from .product_type import ProductType, DOESNT_NEED_API
from .location_coordinates import LocationCoordinates

from .validation.warnings import WARNINGS
from .validation.errors import ERRORS


@dataclass_to_json
@dataclass
class MetaData:
    name: str = ''
    source_id: str = ''
    alternate_names: list = field(default_factory=list)
    location: str = ''
    reference: Reference = field(default_factory=Reference)
    sample_date: str = ''
    product_type: ProductType = ''
    API: float = None
    comments: str = ''
    labels: list = field(default_factory=list)
    model_completeness: float = None
    location_coordinates: LocationCoordinates = None
    gnome_suitable: bool = None

    def validate(self):
        msgs = []
        # check for API
        api = self.API
        if api is None:
            if self.product_type in DOESNT_NEED_API:
                msgs.append(WARNINGS["W004"])
            else:
                msgs.append(ERRORS["E030"])
        else:
            if not (-60.0 < api < 100):  # somewhat arbitrary limits
                msgs.append(WARNINGS["W005"].format(api=api))

        # Check for a reasonable name
        # right now, reasonable is more than 5 characters -- we may want to add more later
        if len(self.name.strip()) < 2:
            msgs.append(WARNINGS["W001"].format(self.name))

        # check sample date is valid
        if self.sample_date:
            try:
                datetime.fromisoformat(self.sample_date)
            except ValueError as err:
                msgs.append(WARNINGS["W011"].format("sample date", self.sample_date, str(err)))
        return msgs


@dataclass_to_json
@dataclass
class SampleMetaData:
    name: str = "Fresh Oil Sample"
    short_name: str = None
    sample_id: str = None
    description: str = None
    fraction_evaporated: MassFraction = None
    boiling_point_range: Temperature = None
