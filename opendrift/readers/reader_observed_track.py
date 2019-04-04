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

from datetime import datetime, timedelta
import logging
import time
import numpy as np
import pandas as pd
import os
from readers.idtt import ReadTrack as rt

from readers.basereader import BaseReader

class Reader(BaseReader):

    name = 'observed_track'
    return_block = False  # Vector based, so checks only individual points

    # Variables (CF standard names) which
    # can be provided by this model/reader
    variables = ['x_observed','y_observed']
    
    def __init__(self, filename=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')
            
        self.name = 'observed_track'
        

#        if name is None:
#            self.name = filename
#        else:
#            self.name = name
        
#        try:
            # Open file, check that everything is ok
        logging.info('Opening dataset: ' + filename)
        data = rt(filename)

        data = data.set_index('datetime')
        #print ' hello data \n'
        #print data
        self.Trajectory = {}

        logging.debug('Assuming geodetic lon/lat projection.')
        self.projected = False
        
     
        # Depth
        self.z = None

        if 'PAL_2012_data' or 'PAL_2017_data' in filename:  #  Need to deal with bad observations
            # Implement a 3 point median filter prior to gap filling to remove outliers.
            # This perfectly preserves monotonic sequences. Need to check if this handles both lat and lon 
            #data.set_index('datetime')
            #data = pd.rolling_median(data,3,1,center=True)   This has been deprecated
            print("Doing median interpolation")
            data = pd.rolling_median(data,3,1,center=True)
         #   print 'Hey ',type(data)
         #   print data
            # Resample to regular one hour intervals.  This leaves NaNs for new time points.
#            idata = data.resample('H',loffset= '30Min')*1;  # *1 needed to realize resampler 
            idata = data.resample('H')*1;  # *1 needed to realize resampler 
            #print type(idata); print idata
            # Replace the NaN values with interpolated.  Alternative to cubic is 'time' for linear interp
            idata.interpolate('cubic',inplace=True);  

            thing = np.reshape(idata['longitude'].values,(-1,1))
            self.Trajectory['longitude'] = thing
            thing = np.reshape(idata['latitude'].values,(-1,1))
            self.Trajectory['latitude'] = thing
            self.times = idata.index; 
            self.start_time = idata.index[0]; print('self.start_time '); print(self.start_time)
            self.end_time = idata.index[-1]; print('self.end_time '); print(self.end_time)
            self.time_step = None

#
        # Read and store min, max and step of x and y
        # Got to fix this!!
#        self.xmin, self.ymin = (min(self.Trajectory.get('longitude')), min(self.Trajectory.get('latitude')))
        self.xmin, self.ymin = (min(self.Trajectory.get('longitude')), min(self.Trajectory.get('latitude')))
        self.xmax, self.ymax = (max(self.Trajectory.get('longitude')), max(self.Trajectory.get('latitude')))
        self.delta_x = None
        self.delta_y = None

        
        # Run constructor of parent Reader class

        super(Reader, self).__init__()


    def get_variables(self, requestedVariables, time_=None, block=False):

        if isinstance(requestedVariables, str):
            requestedVariables = [requestedVariables]

        #self.check_arguments(requestedVariables, time)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time_)

        variables = {}
        variables['x_observed'] = self.Trajectory['longitude'][indxTime]
        variables['y_observed'] = self.Trajectory['latitude'][indxTime]

        # print 'variables[x_observed]'
        # print variables['x_observed']

        return variables


