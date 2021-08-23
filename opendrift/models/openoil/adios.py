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

import requests
from typing import List

ADIOS = "https://adios.orr.noaa.gov/api/oils/"

# The SSL configuration of the ADIOS database does not work outside the browser it seems. Please see
# this issue: https://github.com/NOAA-ORR-ERD/adios_oil_database/issues/2.
VERIFY = False


class Oil:
    _id: str
    _type: str
    name: str
    API: float
    gnome_suitable: bool
    labels: List[str]
    location: str
    model_completeness: float
    product_type: str
    sample_date: str

    @staticmethod
    def from_json(d) -> 'Oil':
        o = Oil()
        o._id = d['_id']
        o._type = d['type']

        meta = d['attributes']['metadata']
        o.name = meta['name']
        o.API = meta['API']
        o.gnome_suitable = meta['gnome_suitable']
        o.labels = meta['labels']
        o.location = meta['location']
        o.model_completeness = meta['model_completeness']
        o.product_type = meta['product_type']
        o.sample_date = meta['sample_date']

        return o

    def __repr__(self):
        return f"[<adios.Oil> {self._id}] {self.name}"


def __oils_per_page__() -> float:
    """
    Number of oils per page.
    """
    m = requests.get(ADIOS, {'limit': 1}, verify=VERIFY).json()
    total = int(m['meta']['total'])
    pages = float(m['meta']['totalPages'])

    return total / pages


def oils(limit=50) -> List[Oil]:
    """
    Get all oils.
    """
    LIMIT = 50
    n = __oils_per_page__()
    assert n > 0

    oils = []

    while len(oils) < limit or limit <= 0:
        p = int(len(oils) / n)  # current page
        o = requests.get(ADIOS, {
            'dir': 'asc',
            'limit': LIMIT,
            'page': p,
            'sort': 'metadata.name'
        },
                         verify=VERIFY).json()
        oils.extend(o['data'])

    oils = oils[:min(limit, len(oils))]
    oils = [Oil.from_json(o) for o in oils]

    return oils
