#!/usr/bin/env python

from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

from readers import Reader

class Reader(Reader):

    # Variables (CF standard names) which
    # can be provided by this model/reader
    variables = ['x_sea_water_velocity',
                 'y_sea_water_velocity',
                 'sea_water_salinity',
                 'sea_water_temperature',
                 'sea_ice_area_fraction']

    def __init__(self):

        f = 'http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be'
        self.name = 'roms_arctic_20km'

        # Open file, check that everything is ok
        self.Dataset = Dataset(f)

        # Read and store projection as proj4 string
        self.proj4 = str(self.Dataset.variables['polar_stereographic'].\
                                    __getattr__('proj4_string'))

        # Read and store min, max and step of x and y 
        x = self.Dataset.variables['X'][:]
        x = x*1000 # NB - fix to convert from km to m!
        self.xmin, self.xmax = x.min(), x.max()
        self.delta_x = np.abs(x[1] - x[0])
        y = self.Dataset.variables['Y'][:]
        y = y*1000 # NB - fix to convert from km to m!
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_y = np.abs(y[1] - y[0])
        # Depth
        self.depths = self.Dataset.variables['depth'][:]

        # Read and store time coverage (of this particular file)
        time = self.Dataset.variables['time'][:]
        self.times = [datetime.utcfromtimestamp(t) for t in time] 
        self.startTime = datetime.utcfromtimestamp(time[0])
        self.endTime = datetime.utcfromtimestamp(time[-1])
        self.timeStep = timedelta(seconds = (time[1] - time[0]))

        # Map variable names to their CF standard names
        self.variableMapping = {
            'sea_water_temperature': 'temperature',
            'sea_water_salinity': 'salinity',
            'x_sea_water_velocity': 'u',
            'y_sea_water_velocity': 'v',
            'sea_ice_area_fraction': 'aice'}

        # Run constructor of parent Reader class
        super(Reader, self).__init__()


    def get_variables(self, requestedVariables, time=None, 
                        x=None, y=None, depth=None):

        if isinstance(requestedVariables, str):
            requestedVariables = [requestedVariables]
        if time == None:
            time = self.startTime # Get data from first timestep, if not given
            indxTime = 0
        else:
            indxTime, nearestTime = self.index_of_closest_time(time)

        self.check_arguments(requestedVariables, time, x, y, depth)

        # Find indices corresponding to requested x and y 
        indx = np.round((x-self.xmin)/self.delta_x).astype(int)
        indy = np.round((y-self.ymin)/self.delta_y).astype(int)

        par = requestedVariables[0] # Temporarily only one variable
        var = self.Dataset.variables[self.variableMapping[par]]
        if par == 'sea_ice_area_fraction':
            return var[indxTime, indy, indx] # Sea ice variable has no depth dimension
        else:
            return var[indxTime, 1, indy, indx]
