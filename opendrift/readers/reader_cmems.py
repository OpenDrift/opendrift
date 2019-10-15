#!/usr/bin/env python
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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

# This reader pre-downloads to local disk netCDF-files with data 
# from CMEMS, # http://marine.copernicus.eu 
# This first version fetches current data from Mercator global ocean model

import logging

import os
from datetime import datetime, timedelta
import numpy as np

from .reader_netCDF_CF_generic import Reader as NCReader


class Reader(NCReader):

    def __init__(self, username, password,
                 motu='motu-client.py',
                 productID='global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh',
                 motu_URL='http://nrtcmems.mercator-ocean.fr/motu-web/Motu', serviceID='GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS',
                 variables=['uo', 'vo'],
                 lon_min=-10, lon_max=10, lat_min=55, lat_max=62,
                 depth_min=0, depth_max=3,
                 time_start=datetime.now(),
                 time_end=datetime.now() + timedelta(days=1)):

        nc_file = 'opendrift_cmems_download.nc'
        cmd = 'python ' + motu + ' --auth-mode=cas -m %s -s %s -d %s -x %s -X %s -y %s -Y %s -z %s -Z %s -t %s -T %s -v uo -v vo -f %s -u %s -p %s' % (
            motu_URL, serviceID, productID,
            lon_min, lon_max, lat_min, lat_max, depth_min, depth_max,
            time_start.strftime('"%Y-%m-%d %H:%M:%S"'),
            time_end.strftime('"%Y-%m-%d %H:%M:%S"'),
            nc_file, username, password)
        logging.info('Downloading file from CMEMS server, using motu-client:')
        logging.info(cmd)
        os.system(cmd)

        super(Reader, self).__init__(filename=nc_file)
