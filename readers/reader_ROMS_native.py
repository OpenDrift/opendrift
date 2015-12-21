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
from netCDF4 import Dataset, num2date
from scipy.ndimage import map_coordinates
from scipy.spatial import KDTree

from readers import Reader


class fakeproj():
    # Since ROMS native grid does not correspond to a SRS and proj4-string,
    # we emulate a pyproj class with needed functions
    def is_latlong(self):
        return False

    def __call__(self, x, y, inverse=False):
        return x, y


class Reader(Reader):

    ROMS_variable_mapping = {
        'mask_psi': 'land_binary_mask',
        'zeta': 'sea_surface_height',
        'u': 'x_sea_water_velocity',
        'v': 'y_sea_water_velocity',
        'temp': 'sea_water_temperature',
        'salt': 'sea_water_salinity',
        'uice': 'sea_ice_x_velocity',
        'vice': 'sea_ice_y_velocity',
        'aice': 'sea_ice_area_fraction',
        'hice': 'sea_ice_thickness'}

    def xy2lonlat(self, x, y):
        y = np.atleast_1d(np.array(y))
        x = np.atleast_1d(np.array(x))
        lon = map_coordinates(self.lon, [x, y], order=1,
                              cval=np.nan, mode='nearest')
        lat = map_coordinates(self.lat, [x, y], order=1,
                              cval=np.nan, mode='nearest')
        return (lon, lat)

    def lonlat2xy(self, lon, lat):
        # nearest point, gives distance and point index
        d,p = self.kdtree.query(zip(np.atleast_1d(lon),
                                    np.atleast_1d(lat)), k=1) 
        i, j = np.unravel_index(p, self.lon.shape)
        return (i, j)

    def __init__(self, filename=None, name=None):

        print '\n\n\n'

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
        
        # Read sigma-coordinate parameters
        Vtransform = self.Dataset.variables['Vtransform'][:]
        Vstretching = self.Dataset.variables['Vstretching'][:]
        theta_s = self.Dataset.variables['theta_s'][:]
        theta_b = self.Dataset.variables['theta_b'][:]
        Tcline = self.Dataset.variables['Tcline'][:]

        # Horizontal oordinates and directions
        self.lat = self.Dataset.variables['lat_rho'][:]
        self.lon = self.Dataset.variables['lon_rho'][:]
        self.angle_between_x_and_east = self.Dataset.variables['angle'][:]

        # Get time coverage
        ocean_time = self.Dataset.variables['ocean_time']
        time_units = ocean_time.getncattr('units')
        self.times = num2date(ocean_time, time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
        else:
            self.time_step = None

        # x and y are rows and columns for unprojected datasets
        self.xmin = 0.
        self.xmax = np.float(len(self.Dataset.dimensions['eta_rho']))
        self.delta_x = 1.
        self.ymin = 0.
        self.ymax = np.float(len(self.Dataset.dimensions['xi_rho']))
        self.delta_y = 1.

        self.proj4 = 'none'
        self.proj = fakeproj()
        self.name = 'roms native'

        # Find all variables having standard_name
        self.variables = []
        for var_name in self.Dataset.variables:
            if var_name in self.ROMS_variable_mapping.keys():
                var = self.Dataset.variables[var_name]
                self.variables.append(self.ROMS_variable_mapping[var_name])

        #import matplotlib.pyplot as plt
        #plt.imshow(h.T)
        #plt.colorbar()
        #plt.show()

        logging.debug('Making KDTree for lon,lat to x,y convertion...')
        self.kdtree = KDTree(zip(self.lon.ravel(), self.lat.ravel()))

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
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int)
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int)

        ## If x or y coordinates are decreasing, we need to flip
        #if self.x[0] > self.x[-1]:
        #    indx = len(self.x) - indx
        #if self.y[0] > self.y[-1]:
        #    indy = len(self.y) - indy
        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = 8
            indx = np.arange(np.max([0, indx.min()-buffer]),
                             np.min([indx.max()+buffer, self.lon.shape[0]]))
            indy = np.arange(np.max([0, indy.min()-buffer]),
                             np.min([indy.max()+buffer, self.lon.shape[1]]))
        else:
            indx[outside[0]] = 0  # To be masked later
            indy[outside[0]] = 0

        variables = {}

        for par in requested_variables:
            varname = [name for name, cf in
                        self.ROMS_variable_mapping.items() if cf == par]
            var = self.Dataset.variables[varname[0]]

            # NB: below x and y are swapped, relative to netCDF_generic-reader
            if var.ndim == 2:
                variables[par] = var[indx, indy]
            elif var.ndim == 3:
                variables[par] = var[indxTime, indx, indy]
            elif var.ndim == 4:
                # Temporarily neglecting depth
                variables[par] = var[indxTime, 0, indx, indy]  # NB 0 was 1
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

        # Return coordinate system orientation, for vector rotation
        variables['angle_between_x_and_east'] = \
            np.degrees(self.angle_between_x_and_east[np.meshgrid(indx, indy)])

        # Store coordinates of returned points
        try:
            variables['depth'] = self.depths[ind_depth]
        except:
            variables['depth'] = None
        if block is True:
            variables['x'] = indx
            variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['time'] = nearestTime
        return variables
