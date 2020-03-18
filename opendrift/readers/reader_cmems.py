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
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

# This reader pre-downloads to local disk netCDF-files with data 
# from CMEMS, # https://marine.copernicus.eu 
# This first version fetches current data from Mercator global ocean model


import os
import logging
import subprocess
from datetime import datetime, timedelta
from xml.etree import ElementTree
import numpy as np
import pyproj
import isodate


logger = logging.getLogger('opendrift')  # using common logger

try:
    import motuclient
except:
    raise ValueError('Motu client must be installed to use reader_cmems.py: python -m pip install motuclient')

from opendrift.readers.reader_netCDF_CF_generic import Reader as NCReader


motu_URL='https://nrt.cmems-du.eu/motu-web/Motu'
products = {
    'NORTHWESTSHELF_ANALYSIS_FORECAST_PHY_004_013-TDS': 
        {'productID': 'MetO-NWS-PHY-hi-CUR',
         'variables': {'uo': 'x_sea_water_velocity',
                       'vo': 'y_sea_water_velocity'}},
    'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS':
        {'productID': 'global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh',
         'variables': {'uo': 'x_sea_water_velocity',
                       'vo': 'y_sea_water_velocity'}}
        }

class Reader(NCReader):

    def __init__(self, cmems_user=None, cmems_password=None,
                 serviceID='GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS',
                 lon_min=None, lon_max=None, lat_min=None, lat_max=None,
                 depth_min=0, depth_max=3,
                 time_start=datetime.now(), ID='',
                 time_end=datetime.now() + timedelta(days=1)):

        if cmems_user is None:
            if 'CMEMS_USER' in os.environ and 'CMEMS_PASSWORD' in os.environ:
                self.cmems_user = os.environ['CMEMS_USER']
                self.cmems_password = os.environ['CMEMS_PASSWORD']
            else:
                raise ValueError('CMEMS username and password must be provided, '
                                 'or stored as environment variables CMEMS_USER and CMEMS_PASSWORD')
        else:
            self.cmems_user = cmems_user
            self.cmems_password = cmems_password

        content = products[serviceID]
        if serviceID not in products:
            raise ValueError('serviceID must be one of: ' + str(products.keys()))
        else:
            self.productID = content['productID']

        self.serviceID = serviceID
        self.variables = list(content['variables'].values())

        # Downloaded data will be stored in this file, to be overwritten by subsequent downloads
        self.nc_file = self.productID + ID + '.nc'

        # Download xml file specifying content
        content_file = self.productID + '.xml'
        #cmd = 'motuclient -D --auth-mode=cas -m %s -s %s -d %s -u %s -p %s -f %s' % (
        cmd = 'motuclient -D -m %s -s %s -d %s -u %s -p %s -f %s' % (
            motu_URL, serviceID, self.productID, self.cmems_user, self.cmems_password, content_file)
        print(cmd.replace(self.cmems_password, '****'))

        if os.path.exists(content_file) and True is False:  # reloading each time, TBD
            print('Reusing contents file')
        else:
            print('Downloading contents file')
            p = subprocess.Popen(cmd.split(' '))
            try:
                p.wait(60)  # timeout in seconds
            except subprocess.TimeoutExpired:
                print('TIEMOUT!')
                p.kill()

            if os.path.exists(content_file):
                print('Success')
            else:
                print('Failure')

        
        #with open(content_file, 'rt') as e:
        #    tree = ElementTree.parse(e)
        #tc = tree.findall('timeCoverage')[0]

        xml = open(content_file, 'rt').read()
        root = ElementTree.fromstring(xml)
        # Time
        #times = root.find('timeCoverage')
        times = root.find('availableTimes')
        start_time, end_time, time_step = (times.text.split('/'))
        start_time = start_time.split('Z')[0]
        end_time = end_time.split('Z')[0]
        self.time_step = isodate.parse_duration(time_step)

        try:
            self.start_time = datetime.strptime(
                start_time, '%Y-%m-%dT%H:%M:%S')
            self.end_time = datetime.strptime(
                end_time, '%Y-%m-%dT%H:%M:%S')
        except:
            self.start_time = datetime.strptime(
                start_time, '%Y-%m-%d')
            self.end_time = datetime.strptime(
                end_time, '%Y-%m-%d')

        # Depths
        try:
            depths = root.find('availableDepths')
            depths = np.array(depths.text.split(';'))
        except:
            depths = np.array(0)
        self.z = np.array([np.float(d.replace('Surface', '0')) for d in depths])
        self.z = -np.abs(self.z)

        # Axes
        axes = root.find('dataGeospatialCoverage')
        for axis in axes:
            desc = axis.attrib['name']
            axisType = axis.attrib['axisType']
            if axisType == 'Lat':
                latmin = np.float(axis.attrib['lower'])
                latmax = np.float(axis.attrib['upper'])
            if axisType == 'Lon':
                lonmin = np.float(axis.attrib['lower'])
                lonmax = np.float(axis.attrib['upper'])
        # Presently only supporting lonlat-projection, which is most common from CMEMS
        self.xmin = lonmin
        self.xmax = lonmax
        self.ymin = latmin
        self.ymax = latmax

        variables = root.find('variables')
        self.variables = []
        for variable in variables:
            standard_name = variable.attrib['standardName']
            self.variables.append(standard_name)

        self.name = 'CMEMS'
        self.proj4 = '+proj=latlong'
        self.time_step = timedelta(hours=1)

        self.name = self.productID
        self.proj4 = '+proj=latlong'

        super(NCReader, self).__init__()

    def prepare(self, extent, start_time, end_time):
        print('Preparing reader ' + self.name)

        #if os.path.exists(nc_file):
        #    try:
        #        r = NCReader(nc_file)
        #        if (r.xmin <= lon_min and r.xmax >= lon_max and
        #            r.ymin <= lat_min and r.ymin <= lat_max and
        #            r.start_time <= time_start and r.end_time >= time_end):
        #            print('Reusing downloaded file: ' + nc_file)
        #            super(Reader, self).__init__(filename=nc_file)
        #            return
        #        else:
        #            print(r.xmin <= lon_min and r.xmax >= lon_max)
        #            print(r.end_time >= time_end)
        #            print(r.start_time <= time_start)
        #            print('NOCOVER...')
        #    except Exception as e:
        #        print(e)
        #        pass
        #else:
        #    print('Does not exist: ' + nc_file)
        ##stop

        lon_min, lat_min, lon_max, lat_max = extent
        time_start = start_time - timedelta(hours=1)  # Some extra coverage
        time_end = end_time + timedelta(hours=1)  # Some extra coverage
        z_epsilon = 1
        depth_min = np.abs(self.z.max()) - 1
        depth_max = np.abs(self.z.min()) + 1
        #  Using only surface current:
        depth_max = np.abs(self.z.max()) + 1

        # TODO
        # Hardcoded for current, presently
        cmd = 'motuclient --auth-mode=cas -m %s -s %s -d %s -x %s -X %s -y %s -Y %s -z %s -Z %s -t %s -T %s -v uo -v vo -f %s -u %s -p %s' % (
            motu_URL, self.serviceID, self.productID,
            lon_min, lon_max, lat_min, lat_max, depth_min, depth_max,
            time_start.strftime('"%Y-%m-%d %H:%M:%S"'),
            time_end.strftime('"%Y-%m-%d %H:%M:%S"'),
            self.nc_file, self.cmems_user, self.cmems_password)
        print(cmd.replace(self.cmems_password, '*****'), 'CMD')
        self.logger.info('Downloading file from CMEMS server, using motu-client:')
        self.logger.info(cmd)
        os.system(cmd)

        super(Reader, self).__init__(filename=self.nc_file)
