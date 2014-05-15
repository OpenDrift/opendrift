#!/usr/bin/env python

from datetime import datetime

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

    def __init__(self, startTime=None):

        print 'Initialising...'
        self.startTime = None
        self.endTime = None
        
        # Construct filename
        if startTime is None: # Find file of today, if time is not given
            startTime = datetime(datetime.now().year,
                                 datetime.now().month,
                                 datetime.now().day)
        f = 'http://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.fc.' + startTime.strftime('%Y%m%d') + '00.nc'

        # Open file, check that everything is ok
        self.Dataset = Dataset(f)

        # Read and store projection as proj4 string
        self.proj4 = str(self.Dataset.variables['projection_stere'].\
                                    __getattr__('proj4'))

        # Read and store min, max and step of x and y 
        x = self.Dataset.variables['X'][:]
        self.xmin = x.min()
        self.xmax = x.max()
        self.delta_x = np.abs(x[1] - x[0])
        y = self.Dataset.variables['Y'][:]
        self.ymin = y.min()
        self.ymax = y.max()
        self.delta_y = np.abs(x[1] - x[0])
        # Depth
        self.depths = self.Dataset.variables['depth'][:]

        # Read and store time coverage (of this particular file)
        # - to be implemented

        # Run constructor of parent PointReader class
        super(Reader, self).__init__()


    def get_parameters(self, parameters, time=None, 
                        x=None, y=None, depth=None):

        self.check_arguments(parameters, time, x, y, depth)

        # Fake northward current with velocity 1 m/s
        if 'eastward_sea_water_velocity' in parameters:
            self.easthward_sea_water_velocity = 0
        if 'northward_sea_water_velocity' in parameters:
            self.northward_sea_water_velocity = 1
