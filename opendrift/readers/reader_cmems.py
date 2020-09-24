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
from netCDF4 import Dataset


logger = logging.getLogger('opendrift')  # using common logger

try:
    import motuclient
except:
    raise ValueError('Motu client must be installed to use reader_cmems.py: python -m pip install motuclient')

from opendrift.readers.reader_netCDF_CF_generic import Reader as NCReader


motu_URL='http://nrt.cmems-du.eu/motu-web/Motu'
products = {
    #  product: { dataset1: {variable_mapping: {var: standard_name}},
    #             dataset2: {variable_mapping: {var: standard_name}} }
    'NORTHWESTSHELF_ANALYSIS_FORECAST_PHY_004_013-TDS': {
        'MetO-NWS-PHY-hi-CUR': {
            'variable_mapping': {'uo': 'x_sea_water_velocity',
                                 'vo': 'y_sea_water_velocity'}}
         },
    'GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS': {
        'global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh': {
            'variable_mapping': {'uo': 'x_sea_water_velocity',
                                 'vo': 'y_sea_water_velocity'}},
        'global-analysis-forecast-phy-001-024-hourly-merged-uv': {
            'variable_mapping': {'utotal': 'x_sea_water_velocity',
                                 'vtotal': 'y_sea_water_velocity',
                                 'vsdx': 'sea_surface_wave_stokes_drift_x_velocity',
                                 'vsdy': 'sea_surface_wave_stokes_drift_y_velocity',
                          }}
        }
    }


class Reader(NCReader):

    def __init__(self, dataset, product=None, variable_mapping=None,
                 cmems_user=None, cmems_password=None,
                 lon_min=None, lon_max=None, lat_min=None, lat_max=None,
                 depth_min=0, depth_max=3,
                 time_start=datetime.now(), ID='',
                 time_end=datetime.now() + timedelta(days=1)):

        if cmems_user is None:
            try:
                import netrc
                n = netrc.netrc()
                self.cmems_user, dummy, self.cmems_password = n.authenticators('cmems')
            except:
                raise ValueError('CMEMS username and password must be provided, '
                                 'or stored in a .netrc file under machine name "cmems"')
        else:
            self.cmems_user = cmems_user
            self.cmems_password = cmems_password

        for productname, product in products.items():
            if productname == dataset:
                raise ValueError(
                    'Product name is provided, please use dataset name.')
            if dataset in product:
                self.product = productname
                self.dataset = dataset
                if variable_mapping is None:
                    self.variable_mapping = product[dataset]['variable_mapping']
                else:  # Using custom variable mapping
                    self.variable_mapping = variable_mapping
                break

        if not hasattr(self, 'dataset'):
            raise ValueError('Dataset not available: ' + dataset)

        self.variables = list(self.variable_mapping.values())

        # Downloaded data will be stored in this file, to be overwritten by subsequent downloads
        self.nc_file = self.dataset + ID + '.nc'

        # Download xml file specifying content
        content_file = self.dataset + ID + '.xml'
        #cmd = 'motuclient -D --auth-mode=cas -m %s -s %s -d %s -u %s -p %s -f %s' % (
        cmd = 'motuclient -D -m %s -s %s -d %s -u %s -p %s -f %s' % (
            motu_URL, self.product, self.dataset,
            self.cmems_user, self.cmems_password, content_file)
        print(cmd.replace(self.cmems_password, '****'))

        print('Downloading contents file')
        p = subprocess.Popen(cmd.split(' '))
        try:
            p.wait(120)  # timeout in seconds
        except subprocess.TimeoutExpired:
            print('TIEMOUT!')
            p.kill()

        xml = open(content_file, 'rt').read()
        root = ElementTree.fromstring(xml)
        # Time
        #times = root.find('timeCoverage')
        times = root.find('availableTimes')
        start_time, end_time, time_step = (times.text.split('/'))[0:3]
        time_step = time_step.split(',')[0]
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
        self.proj4 = '+proj=latlong'
        self.xmin = lonmin
        self.xmax = lonmax
        self.ymin = latmin
        self.ymax = latmax

        #variables = root.find('variables')
        #self.variables = []
        #for variable in variables:
        #    standard_name = variable.attrib['standardName']
        #    self.variables.append(standard_name)

        self.name = self.dataset

        super(NCReader, self).__init__()

    def prepare(self, extent, start_time, end_time):
        print('Preparing reader ' + self.name)

        lon_min, lat_min, lon_max, lat_max = extent
        time_start = start_time - timedelta(hours=1)  # Some extra coverage
        time_end = end_time + timedelta(hours=1)  # Some extra coverage
        z_epsilon = 1
        depth_min = np.abs(self.z.max()) - 1
        depth_max = np.abs(self.z.min()) + 1
        #  Using only surface current:
        depth_max = np.abs(self.z.max()) + 1

        varstring = ''.join(['-v %s ' % v for v in self.variable_mapping])
        cmd = 'motuclient --auth-mode=cas -m %s -s %s -d %s -x %s -X %s -y %s -Y %s -z %s -Z %s -t %s -T %s %s -f %s -u %s -p %s' % (
            motu_URL, self.product, self.dataset,
            lon_min, lon_max, lat_min, lat_max, depth_min, depth_max,
            time_start.strftime('"%Y-%m-%d %H:%M:%S"'),
            time_end.strftime('"%Y-%m-%d %H:%M:%S"'),
            varstring,
            self.nc_file, self.cmems_user, self.cmems_password)
        print(cmd.replace(self.cmems_password, '*****'), 'CMD')
        self.logger.info('Downloading file from CMEMS server, using motu-client:')
        self.logger.info(cmd)
        os.system(cmd)

        # Update standard_name attribute with provided variable mapping
        d = Dataset(self.nc_file, 'a')
        for var, val in self.variable_mapping.items():
            print('Setting standard_name of %s to %s' % (var, val))
            d.variables[var].standard_name = val

        super(Reader, self).__init__(filename=self.nc_file)
