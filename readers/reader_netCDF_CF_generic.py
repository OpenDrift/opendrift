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

import os
import logging
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date

from reader import Reader


class Reader(Reader):

    def __init__(self, filename=None, name=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        if name is None:
            self.name = filename
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            logging.info('Opening dataset: ' + filename)
            self.Dataset = MFDataset(filename, 'r')
        except:
            try:
                logging.info('Could not open with MFDataset, '
                             'trying with Dataset:')
                self.Dataset = Dataset(filename, 'r')
            except:
                raise ValueError('Could not open ' + filename +
                                 ' with netCDF4 library')

        logging.debug('Finding map projection.')
        # Find projection (variable which as proj4 string)
        for var_name in self.Dataset.variables:
            var = self.Dataset.variables[var_name]
            for att in var.ncattrs():
                if 'proj4' in att:
                    self.proj4 = str(var.__getattr__(att))
                    grid_mapping = var_name
        if not hasattr(self, 'proj4'):
            self.proj4 = '+proj=longlat'  # Assuming lonlat
            #raise ValueError('Did not find any proj4 string in dataset')

        logging.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        for var_name in self.Dataset.variables:
            var = self.Dataset.variables[var_name]
            if var.ndim > 1:
                continue  # Coordinates must be 1D-array
            attributes = var.ncattrs()
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            if 'standard_name' in attributes:
                standard_name = var.__dict__['standard_name']
            if 'long_name' in attributes:
                long_name = var.__dict__['long_name']
            if 'axis' in attributes:
                axis = var.__dict__['axis']
            if 'units' in attributes:
                units = var.__dict__['units']
            if standard_name == 'longitude' or \
                    long_name == 'longitude' or \
                    axis == 'X' or \
                    standard_name == 'projection_x_coordinate':
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                x = var[:]*unitfactor
                self.unitfactor = unitfactor
                self.numx = var.shape[0]
            if standard_name == 'latitude' or \
                    long_name == 'latitude' or \
                    axis == 'Y' or \
                    standard_name == 'projection_y_coordinate':
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                y = var[:]*unitfactor
                self.numy = var.shape[0]
            if standard_name == 'depth' or axis == 'Z':
                if not 'positive' in var.ncattrs() or var.__dict__['positive'] == 'up':
                    self.z = var[:]
                else:
                    self.z = -var[:]
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

        if not 'x' in locals():
            raise ValueError('Did not find x-coordinate variable')
        if not 'y' in locals():
            raise ValueError('Did not find y-coordinate variable')
        self.xmin, self.xmax = x.min(), x.max()
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_x = np.abs(x[1] - x[0])
        self.delta_y = np.abs(y[1] - y[0])
        self.x = x  # Store coordinate vectors
        self.y = y

        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname, 'depth']:
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]
            attributes = var.ncattrs()
            if 'standard_name' in attributes:
                standard_name = str(var.__dict__['standard_name'])
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)

        self.variables = self.variable_mapping.keys()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        if hasattr(self, 'z'):
            # Find z-index range
            # NB: may need to flip if self.z is ascending
            indices = np.searchsorted(-self.z, [-z.min(), -z.max()])
            indz = np.arange(np.maximum(0, indices.min() - 1
                                        - self.verticalbuffer),
                             np.minimum(len(self.z), indices.max() + 1
                                        + self.verticalbuffer))
            if len(indz) == 1:
                indz = indz[0]  # Extract integer to read only one layer
        else:
            indz = 0

        # Find indices corresponding to requested x and y
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int)
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int)
        # If x or y coordinates are decreasing, we need to flip
        if self.x[0] > self.x[-1]:
            indx = len(self.x) - indx
        if self.y[0] > self.y[-1]:
            indy = len(self.y) - indy
        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            indx = np.arange(np.max([0, indx.min()-buffer]),
                             np.min([indx.max()+buffer, self.numx]))
            indy = np.arange(np.max([0, indy.min()-buffer]),
                             np.min([indy.max()+buffer, self.numy]))
        else:
            indx[outside[0]] = 0  # To be masked later
            indy[outside[0]] = 0

        variables = {}

        for par in requested_variables:
            var = self.Dataset.variables[self.variable_mapping[par]]

            if var.ndim == 2:
                variables[par] = var[indy, indx]
            elif var.ndim == 3:
                variables[par] = var[indxTime, indy, indx]
            elif var.ndim == 4:
                variables[par] = var[indxTime, indz, indy, indx]
            else:
                raise Exception('Wrong dimension of variable: '
                                + self.variable_mapping[par])

            # If 2D array is returned due to the fancy slicing methods
            # of netcdf-python, we need to take the diagonal
            if variables[par].ndim > 1 and block is False:
                variables[par] = variables[par].diagonal()

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par], ndmin=2, mask=False)
            if block is False:
                variables[par].mask[outside[0]] = True

        # Store coordinates of returned points
        try:
            variables['z'] = self.z[indz]
        except:
            variables['z'] = None
        if block is True:
            variables['x'] = \
                self.Dataset.variables[self.xname][indx]*self.unitfactor
            # Subtracting 1 from indy (not indx) makes Norkyst800
            # fit better with GSHHS coastline - but unclear why
            variables['y'] = \
                self.Dataset.variables[self.yname][indy]*self.unitfactor
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['time'] = nearestTime
        return variables
