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

import logging
logger = logging.getLogger(__name__)
import requests
from typing import List

from .oil import ThinOil, OpendriftOil

ADIOS = "https://adios.orr.noaa.gov/api/oils/"

# The SSL configuration of the ADIOS database does not work outside the browser it seems. Please see
# this issue: https://github.com/NOAA-ORR-ERD/adios_oil_database/issues/2 .
VERIFY = False

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

def get_full_oil_from_id(_id) -> 'OpendriftOil':
    logger.debug(f"Fetching full oil: {_id}")
    o = requests.get(f"{ADIOS}/{_id}", verify=VERIFY).json()
    return OpendriftOil(o)

