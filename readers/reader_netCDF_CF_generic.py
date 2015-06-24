#!/usr/bin/env python

import os
import logging
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset, num2date

from readers import Reader


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
        # Find x, y and depth coordinates
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
                standard_name = var.getncattr('standard_name')
            if 'long_name' in attributes:
                long_name = var.getncattr('long_name')
            if 'axis' in attributes:
                axis = var.getncattr('axis')
            if 'units' in attributes:
                units = var.getncattr('units')
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
                self.depths = var[:]

        if not 'x' in locals():
            raise ValueError('Did not find x-coordinate variable')
        if not 'y' in locals():
            raise ValueError('Did not find y-coordinate variable')
        self.xmin, self.xmax = x.min(), x.max()
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_x = np.abs(x[1] - x[0])
        self.delta_y = np.abs(y[1] - y[0])

        # Read and store time coverage (of this particular file)
        time = self.Dataset.variables['time'][:]  # Assuming name 'time'
        time_units = self.Dataset.variables['time'].getncattr('units')
        self.times = num2date(time, time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
        else:
            self.time_step = None

        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname, 'depth']:
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]
            attributes = var.ncattrs()
            if 'standard_name' in attributes:
                standard_name = str(var.getncattr('standard_name'))
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)

        self.variables = self.variable_mapping.keys()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, depth=None, block=False):

        requested_variables, time, x, y, depth, outside = self.check_arguments(
            requested_variables, time, x, y, depth)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        try:
            ind_depth = self.index_of_closest_depths(depth)[0]
        except:
            ind_depth = 0  # For datasets with no vertical dimension

        # Find indices corresponding to requested x and y
        indx = np.round((x-self.xmin)/self.delta_x).astype(int)
        indy = np.round((y-self.ymin)/self.delta_y).astype(int)
        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = 8
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
                # Temporarily neglecting depth
                variables[par] = var[indxTime, 0, indy, indx]  # NB 0 was 1
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
            variables['depth'] = self.depths[ind_depth]
        except:
            variables['depth'] = None
        if block is True:
            variables['x'] = \
                self.Dataset.variables[self.xname][indx]*self.unitfactor
            variables['y'] = \
                self.Dataset.variables[self.yname][indy]*self.unitfactor
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['time'] = nearestTime
        return variables
