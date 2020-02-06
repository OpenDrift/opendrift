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

#####################################
# NOTE:
# This reader is under development,
# and presently not fully functional
#####################################


import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator

from basereader import BaseReader, pyproj


class Reader(BaseReader):

    def __init__(self, filename=None, name=None, gridfile=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')
        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        # Due to misspelled standard_name in
        # some (Akvaplan-NIVA) FVCOM files
        variable_aliases = {
            'eastward_sea_water_velocity': 'x_sea_water_velocity',
            'Northward_sea_water_velocity': 'y_sea_water_velocity',
            'eastward wind': 'x_wind',
            'northward wind': 'y_wind'
            }

        # Mapping FVCOM variable names to CF standard_name
        fvcom_mapping = {
            'um': 'x_sea_water_velocity',
            'vm': 'y_sea_water_velocity'}

        self.return_block = True

        try:
            # Open file, check that everything is ok
            self.logger.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                self.logger.info('Opening files with MFDataset')
                self.Dataset = MFDataset(filename)
            else:
                self.logger.info('Opening file with Dataset')
                self.Dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

        # We are reading and using lon/lat arrays,
        # and not any projected coordinates
        self.proj4 =  '+proj=latlong'

        self.logger.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        # first check if we have specified a separate grid file
        if gridfile is None: 
            self.gridfile = self.Dataset
        else:
            self.gridfile = Dataset(gridfile)
            self.logger.info('Opening Grid file')
        # now check content of grid- or datafile
        for var_name in self.gridfile.variables:
            var = self.gridfile.variables[var_name]
            if var.ndim > 1:
                continue  # Coordinates must be 1D-array
            attributes = var.ncattrs()
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            CoordinateAxisType = ''
            if 'standard_name' in attributes:
                standard_name = var.__dict__['standard_name']
            if 'long_name' in attributes:
                long_name = var.__dict__['long_name']
            if 'axis' in attributes:
                axis = var.__dict__['axis']
            if 'grid' in attributes:
                grid = var.__dict__['grid']
            if 'units' in attributes:
                units = var.__dict__['units']
            if '_CoordinateAxisType' in attributes:
                CoordinateAxisType = var.__dict__['_CoordinateAxisType']
            # read FVCOM Elements/Center grid ( for u and v):
            if standard_name == 'longitude' and grid == 'Elems' or \
                    var_name == 'lonc':
                self.xname = var_name
                self.numx = var.shape[0]
                x = var[:]
            if standard_name == 'latitude' and grid == 'Elems' or \
                    var_name == 'latc':
                self.yname = var_name
                self.numy = var.shape[0]
                y = var[:]
            if var_name == 'siglayz_center' and grid == 'Elems':
                if 'positive' not in var.ncattrs() or \
                        var.__dict__['positive'] == 'up':
                    self.z = var[:]
                else:
                    self.z = -var[:]

            # todo: read FVCOM Vertices grid ( for tracers)
            #
            
        self.lon = x
        self.lat = y

        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname, 'depth']:
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]
            attributes = var.ncattrs()
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            CoordinateAxisType = ''
            if 'standard_name' in attributes:
                standard_name = var.__dict__['standard_name']
            if 'long_name' in attributes:
                long_name = var.__dict__['long_name']
            if 'axis' in attributes:
                axis = var.__dict__['axis']
            if 'grid' in attributes:
                grid = var.__dict__['grid']
            if 'units' in attributes:
                units = var.__dict__['units']

            if standard_name == 'time' or axis == 'T' or var_name == 'time':
                # Read and store time coverage (of this particular file)
                time = var[:]
                time_units = units
                self.times = num2date(time, time_units)
                self.start_time = self.times[0]
                self.end_time = self.times[-1]
                if len(self.times) > 1:
                    self.time_step = self.times[1] - self.times[0]
                else:
                    self.time_step = None

            if 'standard_name' in attributes:
                standard_name = str(var.__dict__['standard_name'])
                if standard_name in variable_aliases:  # Mapping if needed
                    standard_name = variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)
            elif var_name in fvcom_mapping:
                self.variable_mapping[fvcom_mapping[var_name]] = \
                    str(var_name)

        self.variables = list(self.variable_mapping.keys())

        self.xmin = self.lon.min()
        self.xmax = self.lon.max()
        self.ymin = self.lat.min()
        self.ymax = self.lat.max()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):

        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)


        # Finding a subset around the particles, so that
        # we do not interpolate more points than is needed.
        # Performance is quite dependent on the given buffer,
        # but it should not be made too small to make sure
        # particles are inside box
        buffer = .1  # degrees around given positions

        lonmin = x.min() - buffer
        lonmax = x.max() + buffer
        latmin = y.min() - buffer
        latmax = y.max() + buffer
        c = np.where((self.lon > lonmin) &
                     (self.lon < lonmax) &
                     (self.lat > latmin) &
                     (self.lat < latmax))[0]
        # Making a lon-lat grid onto which data is interpolated
        lonstep = .0004   # hardcoded for now
        latstep = .0002   # hardcoded for now
        lons = np.arange(lonmin, lonmax, lonstep)
        lats = np.arange(latmin, latmax, latstep)
        lonsm, latsm = np.meshgrid(lons, lats)

        # Initialising dictionary to contain data
        variables = {'x': lons, 'y': lats, 'z': z,
                     'time': nearestTime}

        # Reader coordinates of subset
        for par in requested_variables:
            var = self.Dataset.variables[self.variable_mapping[par]]
            print (var)
            if var.ndim == 1:
                data = var[c]
            elif var.ndim == 2:
                data = var[indxTime,c]
            elif var.ndim == 3:
                data = var[indxTime,-1,c]
            else:
                raise ValueError('Wrong dimension of %s: %i' %
                                 (var_name, var.ndim))

            print (data)

            if 'interpolator' not in locals():
                self.logger.debug('Making interpolator...')
                interpolator = LinearNDInterpolator((self.lat[c],
                                                     self.lon[c]),
                                                    data)
            else:
                # Re-use interpolator for other variables
                interpolator.values[:,0] = data
            interpolator((0,0))

            variables[par] = interpolator(latsm, lonsm)

        return variables
