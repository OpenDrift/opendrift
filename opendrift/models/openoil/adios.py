# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2021, Gaute Hope, MET Norway
"""
Interface to the ADIOS oil database.
"""

import logging

logger = logging.getLogger(__name__)
import requests
from typing import List

ADIOS = "https://adios.orr.noaa.gov/api/oils/"

# The SSL configuration of the ADIOS database does not work outside the browser it seems. Please see
# this issue: https://github.com/NOAA-ORR-ERD/adios_oil_database/issues/2 .
VERIFY = False


class ThinOil:
    """
    Basic Oil object for listing. Upgrade to `class:Oil` object for useful methods.
    """
    id: str
    type: str
    name: str
    API: float
    gnome_suitable: bool
    labels: List[str]
    location: str
    model_completeness: float
    product_type: str
    sample_date: str

    def __init__(self, _id, _type, name, API, gnome_suitable, labels, location,
                 model_completeness, product_type, sample_date):
        self.id = _id
        self.type = _type
        self.name = name
        self.API = API
        self.gnome_suitable = gnome_suitable
        self.labels = labels
        self.location = location
        self.model_completeness = model_completeness
        self.product_type = product_type
        self.sample_date = sample_date

    @staticmethod
    def from_json(d) -> 'ThinOil':
        return ThinOil(d['_id'], d['type'], **d['attributes']['metadata'])

    def __repr__(self):
        return f"[<adios.ThinOil> {self.id}] {self.name}"

    def is_generic(self):
        return 'GENERIC' in self.name

    def make_full(self) -> 'Oil':
        """
        Fetch the full oil from ADIOS.
        """
        return Oil.from_id(self.id)


class Oil(ThinOil):
    data: dict

    def __init__(self, o):
        self.data = o

        data = o['data']
        meta = data['attributes']['metadata']
        self.id = data['_id']
        self.name = meta['name']

        from pprint import pp
        pp(o)

    @staticmethod
    def from_id(_id) -> 'Oil':
        logger.debug(f"Fetching full oil: {_id}")
        o = requests.get(f"{ADIOS}/{_id}", verify=VERIFY).json()
        return Oil(o)

    def __repr__(self):
        return f"[<adios.Oil> {self.id}] {self.name}"

    def density_at_temp(self, t) -> float:
        raise NotImplementedError

    def kvis_at_temp(self, t) -> float:
        raise NotImplementedError

    @property
    def mass_fraction(self) -> float:
        raise NotImplementedError

    def oil_water_surface_tension(self) -> List:
        raise NotImplementedError

    @property
    def bulltime(self) -> float:
        raise NotImplementedError

    @property
    def bullwinkle(self) -> float:
        raise NotImplementedError

    @property
    def emulsion_water_fraction_max(self) -> float:
        raise NotImplementedError

    def vapor_pressure(self, t) -> float:
        raise NotImplementedError

    @property
    def molecular_weight(self) -> float:
        raise NotImplementedError

    @property
    def k0y(self) -> float:
        # XXX: No idea what this is.
        raise NotImplementedError


def oils(limit=50, query='') -> List[ThinOil]:
    """
    Get all oils.

    Args:

        limit: number of oils to retrieve, <= 0 means all available.
        query: search string (name, id, location).


    Returns:

        List of `class:ThinOil`s.
    """
    # The batch size seems to be maximum 205 at the moment.
    MAX_BATCH_SZ = 205

    batch = min(MAX_BATCH_SZ, limit)

    oils = []

    # XXX: This fails when batch size is less or unequal to `batch`.
    while len(oils) < limit or limit <= 0:
        p = int(len(oils) / batch)  # next page, XXX: check for off-by-one?
        logging.debug(
            f"Requesting list of oils from ADIOS, oils: {len(oils)} of {limit}, page: {p}"
        )

        o = requests.get(ADIOS, {
            'dir': 'asc',
            'limit': batch,
            'page': p,
            'sort': 'metadata.name',
            'q': query
        },
                         verify=VERIFY).json()

        oils.extend(o['data'])

        if float(o['meta']['totalPages']) <= p:
            break

    limit = len(oils) if limit <= 0 else limit
    oils = oils[:min(limit, len(oils))]
    oils = [ThinOil.from_json(o) for o in oils]

    return oils
