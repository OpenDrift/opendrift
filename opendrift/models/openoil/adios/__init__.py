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

from . import api
from . import oil
from .api import oils, find_full_oil_from_name, get_full_oil_from_id

oil_name_alias = {
        'GENERIC BUNKER C': 'Bunker C [1987]',
        'GENERIC DIESEL': 'DIESEL',
        }

