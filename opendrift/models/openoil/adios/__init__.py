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

from . import dirjs
from . import oil
from .dirjs import get_oil_names, oils, find_full_oil_from_name, get_full_oil_from_id

oil_name_alias = {
        'GENERIC BUNKER C': 'Bunker C [1987]',
        'GENERIC DIESEL': 'DIESEL',
        'GENERIC LIGHT CRUDE': 'DIESEL',
        'GENERIC HEAVY CRUDE': 'DIESEL',
        'GENERIC JET FUEL': 'DIESEL',
        'OSEBERG SOR 2000': 'DIESEL',
        'KVITEBJORN 2019': 'DIESEL',
        'VISUND CRUDE OIL 2020': 'DIESEL',
        'IRIS CONDENSATE 2020': 'DIESEL',
        'FOGELBERG CONDENSATE 2021': 'DIESEL',
        'VISUND SOR CONDENSATE 2020': 'DIESEL',
        'GENERIC HOME HEATING OIL': 'DIESEL',
        'VEGA CONDENSATE 2015': 'DIESEL',
        'AASTA HANSTEEN BLEND 2020': 'DIESEL',
        'GENERIC INTERMEDIATE FUEL OIL 300': 'DIESEL',
        'GINA KROG CRUDE (2018)': 'DIESEL',
        'GENERIC MEDIUM CRUDE': 'DIESEL',
        'SKOGUL 2020': 'DIESEL',
        'GENERIC GASOLINE': 'DIESEL',
        'DVALIN 2020': 'DIESEL',
        'ATLA KONDENSAT 2013': 'DIESEL',
        'ODA 2019': 'DIESEL',
        'GUDRUN 2019': 'DIESEL',
        'GENERIC INTERMEDIATE FUEL OIL 180': 'DIESEL',
        'GENERIC FUEL OIL No. 6': 'DIESEL',
        'GENERIC HEAVY FUEL OIL': 'DIESEL',
        'FROSK 2020': 'DIESEL',
        'GENERIC FUEL OIL No.2': 'DIESEL',
        'FENJA (PIL) 2015': 'DIESEL',
        'MARULK (2014)': 'DIESEL',
        'OSEBERG OST 1998': 'DIESEL',
 }

