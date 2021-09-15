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
        'GENERIC LIGHT CRUDE': None,
        'GENERIC HEAVY CRUDE': None,
        'GENERIC JET FUEL': 'FUEL OIL NO.1 (JET FUEL A)',
        'OSEBERG SOR 2000': 'OSEBERG SOR 2013', # TODO: added by us?
        'KVITEBJORN 2019': 'KVITEBJORN 2009',   # TODO: added by us?
        'VISUND CRUDE OIL 2020': 'VISUND 2009', # TODO: added by us?
        'IRIS CONDENSATE 2020': None, # TODO: added by us?
        'FOGELBERG CONDENSATE 2021': None, # TODO: added by us?
        'VISUND SOR CONDENSATE 2020': None, # TODO: added by us?
        'GENERIC HOME HEATING OIL': 'HOME HEATING OIL',
        'VEGA CONDENSATE 2015': None, # TODO: added by us?
        'AASTA HANSTEEN BLEND 2020': None, # TODO: added by us?
        'GENERIC INTERMEDIATE FUEL OIL 300': 'INTERMEDIATE FUEL OIL 300',
        'GINA KROG CRUDE (2018)': None, # TODO: added by us?
        'GENERIC MEDIUM CRUDE': None,
        'SKOGUL 2020': None, # TODO: added by us?
        'GENERIC GASOLINE': 'GASOLINE (LEADED)',
        'DVALIN 2020': None, # TODO: added by us?
        'ATLA KONDENSAT 2013': 'ALTA 2016',  # TODO: added by us?
        'ODA 2019': None, # TODO: added by us?
        'GUDRUN 2019': 'GUDRUN 2012', # TODO: added by us?
        'GENERIC INTERMEDIATE FUEL OIL 180': 'INTERMEDIATE FUEL OIL 180 (SOCSEX)',
        'GENERIC FUEL OIL No. 6': 'FUEL OIL NO.6',
        'GENERIC HEAVY FUEL OIL': None,
        'FROSK 2020': None,
        'GENERIC FUEL OIL No.2': 'FUEL OIL NO.2',
        'FENJA (PIL) 2015': None,
        'MARULK (2014)': None,
        'OSEBERG OST 1998': None,
 }

