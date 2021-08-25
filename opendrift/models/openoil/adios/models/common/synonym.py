from dataclasses import dataclass
from ..common.utilities import dataclass_to_json


@dataclass_to_json
@dataclass
class Synonym:
    name: str
