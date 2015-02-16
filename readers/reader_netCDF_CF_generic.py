#!/usr/bin/env python

import os
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

from readers import Reader


class Reader(Reader):

    def __init__(self, filename=None, name=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')
        #if not os.path.exists(filename):
        #    raise IOError('File does not exist: ' + filename)

        if name is None:
            self.name = filename
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            self.Dataset = Dataset(filename)
        except:
            raise ValueError('Could not open ' + filename +
                             ' with netCDF4 library')

        # Find projection (variable which as proj4 string)
        for varName in self.Dataset.variables:
            var = self.Dataset.variables[varName]
            for att in var.ncattrs():
                if 'proj4' in att:
                    self.proj4 = str(var.__getattr__(att))
                    grid_mapping = varName
        if not hasattr(self, 'proj4'):
            raise ValueError('Did not find any proj4 string in dataset')

        # Find x, y and depth coordinates
        for varName in self.Dataset.variables:
            var = self.Dataset.variables[varName]
            try:
                att = var.getncattr('standard_name')
                if att == 'projection_x_coordinate':
                    self.xname = varName
                    # Fix for units; should ideally use udunits package
                    if var.getncattr('units') == 'km':
                        unitfactor = 1000
                    else:
                        unitfactor = 1
                    x = var[:]*unitfactor
                if att == 'projection_y_coordinate':
                    self.yname = varName
                    y = var[:]*unitfactor
                if att == 'depth':
                    self.depths = var[:]
            except:
                pass

        if not 'x' in locals():
            raise ValueError('Did not find x-coordinate variable')
        self.xmin, self.xmax = x.min(), x.max()
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_x = np.abs(x[1] - x[0])
        self.delta_y = np.abs(y[1] - y[0])

        # Read and store time coverage (of this particular file)
        time = self.Dataset.variables['time'][:]  # Assuming name 'time'
        self.times = [datetime.utcfromtimestamp(t) for t in time]
        self.startTime = datetime.utcfromtimestamp(time[0])
        self.endTime = datetime.utcfromtimestamp(time[-1])
        self.timeStep = timedelta(seconds=(time[1] - time[0]))

        # Find all variables (with given grid mapping)
        self.variableMapping = {}
        for varName in self.Dataset.variables:
            var = self.Dataset.variables[varName]
            try:
                if var.getncattr('grid_mapping') == grid_mapping:
                    stdName = str(var.getncattr('standard_name'))
                    self.variableMapping[stdName] = str(varName)
            except:
                pass

        self.variables = self.variableMapping.keys()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, depth=None, block=False):

        requestedVariables, time, x, y, depth, outside = self.check_arguments(
            requestedVariables, time, x, y, depth)

        indxTime, nearestTime = self.index_of_closest_time(time)

        ind_depth = self.index_of_closest_depths(depth)[0]

        # Find indices corresponding to requested x and y
        indx = np.round((x-self.xmin)/self.delta_x).astype(int)
        indy = np.round((y-self.ymin)/self.delta_y).astype(int)
        if block is True:
            # Adding buffer of 1, to be checked
            indx = np.arange(indx.min()-1, indx.max()+1)
            indy = np.arange(indy.min()-1, indy.max()+1)
        else:
            indx[outside[0]] = 0  # To be masked later
            indy[outside[0]] = 0

        variables = {}

        for par in requestedVariables:
            var = self.Dataset.variables[self.variableMapping[par]]
            # Hardcoded checks for variable dimensions below,
            # this should be generalised later
            if par == 'x_wind' or par == 'y_wind':
                variables[par] = var[indxTime, indy, indx]
            elif par == 'sea_floor_depth_below_sea_level':
                variables[par] = var[indy, indx]
            else:
                # Temporarily neglecting depth
               variables[par] = var[indxTime, 1, indy, indx]

            # If 2D array is returned due to the fancy slicing methods
            # of netcdf-python, we need to take the diagonal
            if variables[par].ndim > 1 and block is False:
                variables[par] = variables[par].diagonal()

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par], mask=False)
            if block == False:
                variables[par].mask[outside[0]] = True

        # Store coordinates of returned points
        variables['depth'] = self.depths[ind_depth]
        if block is True:
            variables['x'] = self.Dataset.variables[self.xname][indx]
            variables['y'] = self.Dataset.variables[self.yname][indy]
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        return variables
