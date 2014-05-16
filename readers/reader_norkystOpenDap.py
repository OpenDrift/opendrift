#!/usr/bin/env python

from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

from readers import PointReader

class Reader(PointReader):

    # Parameters (CF standard names) which
    # can be provided by this model/reader
    parameters = ['x_sea_water_velocity',
                  'y_sea_water_velocity',
                  'sea_water_salinity',
                  'sea_water_temperature']

    def __init__(self):

        print 'Initialising...'
        
        # Construct filename (updated daily, hence using present time)
        fileTime = datetime(datetime.now().year,
                            datetime.now().month,
                            datetime.now().day)
        f = 'http://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.fc.' + fileTime.strftime('%Y%m%d') + '00.nc'

        # Open file, check that everything is ok
        self.Dataset = Dataset(f)

        # Read and store projection as proj4 string
        self.proj4 = str(self.Dataset.variables['projection_stere'].\
                                    __getattr__('proj4'))

        # Read and store min, max and step of x and y 
        x = self.Dataset.variables['X'][:]
        self.xmin, self.xmax = x.min(), x.max()
        self.delta_x = np.abs(x[1] - x[0])
        y = self.Dataset.variables['Y'][:]
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_y = np.abs(y[1] - y[0])
        # Depth
        self.depths = self.Dataset.variables['depth'][:]

        # Read and store time coverage (of this particular file)
        # - to be implemented
        time = self.Dataset.variables['time'][:]
        self.startTime = datetime.utcfromtimestamp(time[0])
        self.endTime = datetime.utcfromtimestamp(time[-1])
        self.timeStep = timedelta(seconds = (time[1] - time[0]))

        # Run constructor of parent PointReader class
        super(Reader, self).__init__()


    def get_parameters(self, requestedParameters, time=None, 
                        x=None, y=None, depth=None):

        if isinstance(requestedParameters, str):
            requestedParameters = [requestedParameters]
        if time == None:
            time = self.startTime # Get data from first timestep, if not given
        self.check_arguments(requestedParameters, time, x, y, depth)

        # Find indices corresponding to requested x and y 
        indx = np.round((x-self.xmin)/self.delta_x).astype(int)
        indy = np.round((y-self.ymin)/self.delta_y).astype(int)

        par = requestedParameters[0]
        if par == 'sea_water_temperature':
            par = 'temperature'

        var = self.Dataset.variables[par]
        return var[1, 1, indy, indx] # Temporarily neglecting time and depth
