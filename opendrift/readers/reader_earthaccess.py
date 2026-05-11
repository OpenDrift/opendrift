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
# Copyright 2026, Knut-Frode Dagestad, MET Norway

import os
import logging
from opendrift.readers.reader_netCDF_CF_generic import Reader as Reader_CF_generic
logger = logging.getLogger(__name__)

class Reader(Reader_CF_generic):
    '''A wrapper around reader_netCDF_CF_generic for NASA eartchaccess datasets'''

    def __init__(self, shortname, start_time=None, end_time=None):

        try:
            import earthaccess
        except:
            raise ValueError('Earthaccess Client is not installed')

        logger.debug('Authenticating...')
        auth = earthaccess.login(strategy="netrc")

        logger.debug('Searching for granules...')
        granules = earthaccess.search_data(short_name=shortname,
                                           temporal=(start_time, end_time))

        open_options = {
            "access": "indirect",  # access to cloud data (faster in AWS with "direct")
            "load": True,  # Load metadata immediately (required for indexing)
            "concat_dim": "time",  # Concatenate files along the time dimension
            "data_vars": "minimal",  # Only load data variables that include the concat_dim
            "coords": "minimal",  # Only load coordinate variables that include the concat_dim
            "compat": "override",  # Avoid coordinate conflicts by picking the first
            "combine_attrs": "override",  # Avoid attribute conflicts by picking the first
        }

        logger.debug('Creating virtual dataset...')
        vds = earthaccess.virtualize(granules, **open_options)

        vds['name'] = str(shortname)

        # Run constructor of parent Reader class
        super(Reader, self).__init__(vds, name=shortname)
