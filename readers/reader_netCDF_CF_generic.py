#!/usr/bin/env python

import os
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

from readers import PointReader

class Reader(PointReader):


    def __init__(self, filename=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')
        if not os.path.exists(filename):
            raise IOError('File does not exist: ' + filename)

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

        # Find x and y coordinates
        for varName in self.Dataset.variables:
            var = self.Dataset.variables[varName]
            try:
                att = var.getncattr('standard_name')
                if att == 'projection_x_coordinate':
                    x = var[:]
                if att == 'projection_y_coordinate':
                    y = var[:]
            except:
                pass

        if not 'x' in locals():
            raise ValueError('Did not find x-coordinate variable')
        self.xmin, self.xmax = x.min(), x.max()
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_x = np.abs(x[1] - x[0])
        self.delta_y = np.abs(y[1] - y[0])

        ## Depth
        #self.depths = self.Dataset.variables['depth'][:]

        # Read and store time coverage (of this particular file)
        time = self.Dataset.variables['time'][:] # Assuming it is called 'time'
        self.times = [datetime.utcfromtimestamp(t) for t in time] 
        self.startTime = datetime.utcfromtimestamp(time[0])
        self.endTime = datetime.utcfromtimestamp(time[-1])
        self.timeStep = timedelta(seconds = (time[1] - time[0]))

        # Find all parameters (with given grid mapping)
        self.parameterMapping = {}
        for varName in self.Dataset.variables:
            var = self.Dataset.variables[varName]
            try:
                if var.getncattr('grid_mapping') == grid_mapping:
                    stdName = str(var.getncattr('standard_name'))
                    self.parameterMapping[stdName] = str(varName)
            except:
                pass

        self.parameters = self.parameterMapping.keys()

        # Run constructor of parent PointReader class
        super(Reader, self).__init__()


    def get_parameters(self, requestedParameters, time=None, 
                        x=None, y=None, depth=None):

        if isinstance(requestedParameters, str):
            requestedParameters = [requestedParameters]
        if time == None:
            time = self.startTime # Get data from first timestep, if not given
            indxTime = 0
        else:
            indxTime, nearestTime = self.index_of_closest_time(time)

        self.check_arguments(requestedParameters, time, x, y, depth)

        # Find indices corresponding to requested x and y 
        indx = np.round((x-self.xmin)/self.delta_x).astype(int)
        indy = np.round((y-self.ymin)/self.delta_y).astype(int)

        par = requestedParameters[0] # Temporarily only one parameter
        var = self.Dataset.variables[self.parameterMapping[par]]
        return var[indxTime, 1, indy, indx] # Temporarily neglecting depth
