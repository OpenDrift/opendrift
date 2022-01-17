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

from importlib import resources
from pathlib import Path
import logging
import json
import itertools
from functools import lru_cache

from .oil import OpendriftOil

logger = logging.getLogger(__name__)


@lru_cache(maxsize=None)
def __get_archive__():
    import lzma

    oils = []

    with resources.open_binary('opendrift.models.openoil.adios',
                               'oils.xz') as archive:
        with lzma.open(archive, 'rt') as c:
            oils = json.load(c)

    # Add additional oils
    for f in resources.contents('opendrift.models.openoil.adios.extra_oils'):
        if Path(f).suffix == '.json':
            o = json.loads(resources.read_text('opendrift.models.openoil.adios.extra_oils', f))
            o = { 'data': { 'attributes' : o } }
            o['data']['_id'] = o['data']['attributes']['oil_id']
            o['data']['attributes']['metadata']['location'] = 'NORWAY'
            logger.debug(f"Adding additional oil: {f}..: {o['data']['_id']}, {o['data']['attributes']['metadata']['name']}")
            oils.append(o)

    return oils

def get_oil_names(location=None):
    a = __get_archive__()

    if location is not None:
        a = filter(
            lambda o: o['data']['attributes']['metadata'].get(
                'location', None) == location, a)

    return [o['data']['attributes']['metadata']['name'] for o in a]


def oils(limit=50, query=''):
    oils = filter(
        lambda o: query in o['data']['attributes']['metadata']['name'],
        __get_archive__())
    return list(OpendriftOil(o) for o in itertools.islice(oils, limit))


def find_full_oil_from_name(name) -> 'OpendriftOil':
    logger.info(f'Querying ADIOS database for oil: {name}')
    o = oils(query=name)

    if len(o) > 1:
        logger.warning(
            f"Several oils found with name: {name}: {[oo.id for oo in o]}, using first."
        )
    elif len(o) < 1:
        raise ValueError(f"No oil found with name: {name}")

    return o[0]


def get_full_oil_from_id(_id) -> 'OpendriftOil':
    logger.debug(f"Fetching full oil: {_id}")
    oils = filter(lambda o: _id == o['data']['_id'], __get_archive__())
    oil = next(OpendriftOil(o) for o in oils)
    return oil
