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
import logging
import json
import itertools
from functools import cache

from .oil import OpendriftOil

logger = logging.getLogger(__name__)


# @cache
# def __get_archive__():
#     archive = resources.files('opendrift.models.openoil.adios').joinpath(
#         'oils.tar.xz').open('rb')

#     with tarfile.open(fileobj=archive, mode='r:xz') as c:
#         return [json.load(c.extractfile(f)) for f in c]


@cache
def __get_archive__():
    import lzma
    with resources.files('opendrift.models.openoil.adios').joinpath(
            'oils.xz').open('rb') as archive:
        with lzma.open(archive, 'rt') as c:
            oils = json.load(c)
            return oils
            # return [OpendriftOil(o) for o in oils]

# @cache
# def __get_archive__():
#     with resources.files('opendrift.models.openoil.adios').joinpath(
#             'oils.json').open('rb') as c:
#         oils = json.load(c)
#         return [OpendriftOil(o) for o in oils]

def get_oil_names():
    return [o['data']['attributes']['metadata']['name'] for o in __get_archive__()]


def oils(limit=50, query=''):
    oils = filter(lambda o: query in o['data']['attributes']['metadata']['name'], __get_archive__())
    return list(OpendriftOil(o) for o in itertools.islice(oils, limit))

def find_full_oil_from_name(name) -> 'OpendriftOil':
    logger.info(f'Querying ADIOS database for oil: {name}')
    o = oils(query = name)

    if len(o) > 1:
        logger.warning(f"Several oils found with name: {name}: {[oo.id for oo in o]}, using first.")
    elif len(o) < 1:
        raise ValueError(f"No oil found with name: {name}")

    return o[0]

def get_full_oil_from_id(_id) -> 'OpendriftOil':
    logger.debug(f"Fetching full oil: {_id}")
    oils = filter(lambda o: _id == o['data']['_id'], __get_archive__())
    oil = next(OpendriftOil(o) for o in oils)
    return oil

