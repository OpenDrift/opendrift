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

#############################################################################################
# Attempt to reproduce the reader using XARRAY rather than MFDataset to read the netcdf files
#############################################################################################

import logging

import numpy as np
# from netCDF4 import Dataset, MFDataset, num2date
import xarray
from opendrift.readers.basereader import BaseReader


class Reader(BaseReader):

    def __init__(self, filename=None, name=None, proj4=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            logging.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                logging.info('Opening files with xarray.open_mfdataset')
                self.Dataset = xarray.open_mfdataset(filename)
            else:
                logging.info('Opening file with xarray.open_dataset')
                self.Dataset = xarray.open_dataset(filename)
        except Exception as e:
            raise ValueError(e)
        logging.debug('Finding coordinate variables.')
        if proj4 is not None:  # If user has provided a projection apriori
            self.proj4 = proj4
        #
        # Find x, y and z coordinates - use xarray's coords
        #
        # using xarray API :
        # coords  Dictionary of xarray.DataArray objects corresponding to coordinate variables
        # data_vars   Dictionary of xarray.DataArray objects corresponding to data variables
        # dims    Mapping from dimension names to lengths.

        for var_name in self.Dataset.coords:
            logging.debug('Parsing coordinate variable: ' +  var_name)
            var = self.Dataset.variables[var_name]
            #if var.ndim > 1:
            #    continue  # Coordinates must be 1D-array
            attributes = var.attrs #attributes = var.ncattrs()
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            CoordinateAxisType = ''
            
            if not hasattr(self, 'proj4'):
                for att in attributes:
                    if 'proj4' in att:
                        self.proj4 = str(var.__getattr__(att))
            if 'standard_name' in attributes:
                standard_name = var.attrs['standard_name']
            if 'long_name' in attributes:
                long_name = var.attrs['long_name']
            if 'axis' in attributes:
                axis = var.attrs['axis']
            if 'units' in attributes:
                units = var.attrs['units']
            if '_CoordinateAxisType' in attributes:
                CoordinateAxisType = var.attrs['_CoordinateAxisType']
            if standard_name == 'longitude' or \
                    CoordinateAxisType == 'Lon' or \
                    long_name == 'longitude':
                self.lon = var.data # var - not sure what is expected when using netCDF4 ??
                lon_var_name = var_name
            if standard_name == 'latitude' or \
                    CoordinateAxisType == 'Lat' or \
                    long_name == 'latitude':
                self.lat = var.data # var
                lat_var_name = var_name
            if axis == 'X' or \
                    standard_name == 'projection_x_coordinate':
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                elif units == '100  km':
                    unitfactor = 100000
                else:
                    unitfactor = 1
                x = var[:]*unitfactor
                self.numx = var.shape[0] 
            if axis == 'Y' or \
                    standard_name == 'projection_y_coordinate':
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                elif units == '100  km':
                    unitfactor = 100000
                else:
                    unitfactor = 1
                self.unitfactor = unitfactor
                y = var[:]*unitfactor
                self.numy = var.shape[0] 
            if standard_name == 'depth' or axis == 'Z':
                if var[:].ndim == 1:
                    if 'positive' not in attributes or \
                            var.attrs['positive'] == 'up':
                        self.z = var[:]
                    else:
                        self.z = -var[:]
            if standard_name == 'time' or axis == 'T' or var_name in ['time', 'vtime']:
                # Read and store time coverage (of this particular file)
                #
                # the reader normally expects num2date(times,units,calendar='standard')= Return datetime objects given numeric time values
                # convert xarray time to datetime objects, or....
                # store as native xarray time and use nice slicing options ??
                
                time_xarray = var #time = var[:] # store xarray time variables - will be used to data selection later on
                # time_units = units
                self.times = time_xarray# self.times = num2date(time, time_units)
                self.start_time = self.times.data[0]# self.start_time = self.times[0]
                self.end_time = self.times.data[-1]# self.end_time = self.times[-1]
                if len(self.times) > 1:
                    self.time_step = self.times[1] - self.times[0]
                else:
                    self.time_step = None
            if standard_name == 'realization':
                self.realizations = var[:]
                logging.debug('%i ensemble members available'
                              % len(self.realizations))
        if 'x' not in locals():
            if self.lon.ndim == 1:
                x = self.lon[:]
                self.xname = lon_var_name
                self.numx = len(x)
            else:
                raise ValueError('Did not find x-coordinate variable')
        if 'y' not in locals():
            if self.lat.ndim == 1:
                y = self.lat[:]
                self.yname = lat_var_name
                self.numy = len(y)
            else:
                raise ValueError('Did not find y-coordinate variable')
        
        if not hasattr(self, 'unitfactor'):
            self.unitfactor = 1
        if 'x' in locals() and 'y' in locals():
            # at this stage (x,y) are xarray variables, could be casted to something consistent with normal reader
            # 
            # self.xmin, self.xmax = x.min(), x.max()
            # self.ymin, self.ymax = y.min(), y.max()
            # self.delta_x = np.abs(x[1] - x[0])
            # self.delta_y = np.abs(y[1] - y[0])
            # rel_delta_x = (x[1::] - x[0:-1])
            # rel_delta_x = np.abs((rel_delta_x.max() -
            #                       rel_delta_x.min())/self.delta_x)
            # rel_delta_y = (y[1::] - y[0:-1])
            # rel_delta_y = np.abs((rel_delta_y.max() -
            #                       rel_delta_y.min())/self.delta_y)

            # using content of x, i.e. x.data 
            self.xmin, self.xmax = x.data.min(), x.data.max()
            self.ymin, self.ymax = y.data.min(), y.data.max()
            self.delta_x = np.abs(x.data[1] - x.data[0])
            self.delta_y = np.abs(y.data[1] - y.data[0])
            rel_delta_x = (x.data[1::] - x.data[0:-1])
            rel_delta_x = np.abs((rel_delta_x.max() -
                                  rel_delta_x.min())/self.delta_x)
            rel_delta_y = (y.data[1::] - y.data[0:-1])
            rel_delta_y = np.abs((rel_delta_y.max() -
                                  rel_delta_y.min())/self.delta_y)

            if rel_delta_x > 0.05:  # Allow 5 % deviation
                print(rel_delta_x)
                print(x[1::] - x[0:-1])
                raise ValueError('delta_x is not constant!')
            if rel_delta_y > 0.05:
                print(rel_delta_y)
                print(y[1::] - y[0:-1])
                raise ValueError('delta_y is not constant!')
            self.x = x  # Store coordinate vectors
            self.y = y

        else:
            if hasattr(self, 'lon') and hasattr(self, 'lat'):
                logging.info('No projection found, using lon/lat arrays')
                self.xname = lon_var_name
                self.yname = lat_var_name
            else:
                raise ValueError('Neither x/y-coordinates or lon/lat arrays found')

        if not hasattr(self, 'proj4'):
            if self.lon.ndim == 1:
                logging.debug('Lon and lat are 1D arrays, assuming latong projection')
                self.proj4 = '+proj=latlong'
            elif self.lon.ndim == 2:
                logging.debug('Reading lon lat 2D arrays, since projection is not given')
                self.lon = self.lon[:]
                self.lat = self.lat[:]
                self.projected = False

        if hasattr(self, 'proj4') and 'latlong' in self.proj4 and hasattr(self, 'xmax') and self.xmax > 360:
            logging.info('Longitudes > 360 degrees, subtracting 360')
            self.xmin -= 360
            self.xmax -= 360
            self.x -= 360
            self.x -= 360

        # change the `variable_aliases` variable to fit with MetOcean convention (note the mapping was initially defined in basereader.py)
        self.variable_aliases = {
             # standard mapping from CF_generic########################################
            'sea_water_potential_temperature': 'sea_water_temperature',
            'eastward_wind': 'x_wind',
            'northward_wind': 'y_wind',
            'x_wind_10m': 'x_wind',
            'y_wind_10m': 'y_wind',
            'eastward_surface_stokes_drift': 'sea_surface_wave_stokes_drift_x_velocity',
            'northward_surface_stokes_drift': 'sea_surface_wave_stokes_drift_y_velocity',
            'sea_water_x_velocity': 'x_sea_water_velocity',
            'sea_water_y_velocity': 'y_sea_water_velocity',
            'eastward_sea_water_velocity': 'x_sea_water_velocity',
            'northward_sea_water_velocity': 'y_sea_water_velocity',
            'surface_eastward_sea_water_velocity': 'x_sea_water_velocity',
            'surface_northward_sea_water_velocity': 'y_sea_water_velocity',
            'eastward_current_velocity': 'x_sea_water_velocity',
            'northward_current_velocity': 'y_sea_water_velocity',
            'eastward_tidal_current': 'x_sea_water_velocity',
            'northward_tidal_current': 'y_sea_water_velocity',
            'eastward_ekman_current_velocity': 'x_sea_water_velocity',
            'northward_ekman_current_velocity': 'y_sea_water_velocity',
            'eastward_geostrophic_current_velocity': 'x_sea_water_velocity',
            'northward_geostrophic_current_velocity': 'y_sea_water_velocity',
            'eastward_eulerian_current_velocity': 'x_sea_water_velocity',
            'northward_eulerian_current_velocity': 'y_sea_water_velocity',
            'surface_geostrophic_eastward_sea_water_velocity': 'x_sea_water_velocity',
            'surface_geostrophic_northward_sea_water_velocity': 'y_sea_water_velocity',
            'surface_geostrophic_eastward_sea_water_velocity_assuming_sea_level_for_geoid': 'x_sea_water_velocity',
            'surface_geostrophic_northward_sea_water_velocity_assuming_sea_level_for_geoid': 'y_sea_water_velocity',
            'surface_eastward_geostrophic_sea_water_velocity_assuming_sea_level_for_geoid': 'x_sea_water_velocity',
            'surface_northward_geostrophic_sea_water_velocity_assuming_sea_level_for_geoid': 'y_sea_water_velocity',
            # 'mask_rho': 'land_binary_mask',
            # 'mask_psi': 'land_binary_mask',
            ###########################################################################
            # standard mapping MetOcean Reader ########################################
            'dep': 'sea_floor_depth_below_sea_level',
            'el': 'sea_surface_height',   # el = et + ssh
            'et': 'sea_surface_height',   # tidal elevation
            'ssh': 'sea_surface_height',  # non-tidal
            'um': 'x_sea_water_velocity', # depth-averaged total / eastward_sea_water_velocity
            'vm': 'y_sea_water_velocity', # depth-averaged total / northward_sea_water_velocity
            'us': 'x_sea_water_velocity', # surface current total / surface_eastward_sea_water_velocity
            'vs': 'y_sea_water_velocity', # surface current total / surface_northward_sea_water_velocity
            'umo': 'x_sea_water_velocity',# depth-averaged non-tidal / eastward_sea_water_velocity_assuming_no_tide
            'vmo': 'y_sea_water_velocity',# depth-averaged non-tidal / northward_sea_water_velocity_assuming_no_tide
            'ut': 'x_sea_water_velocity', # tide eastward_tidal_current
            'vt': 'y_sea_water_velocity', # tide northward_tidal_current
            'uo': 'x_sea_water_velocity', # non-tidal current 3D / eastward_geostrophic_current_velocity geostrophic_eastward_sea_water_velocity
            'vo': 'y_sea_water_velocity', # non-tidal current 3D / northward_geostrophic_current_velocity
            'u': 'x_sea_water_velocity',  # total current 3D
            'v': 'y_sea_water_velocity',  # total current 3D
            'w': 'upward_sea_water_velocity',
            'uso': 'x_sea_water_velocity', #surface layer non-tidal current (5 m below sea level) / surface_geostrophic_eastward_sea_water_velocity
            'vso': 'y_sea_water_velocity', #surface layer non-tidal current (5 m below sea level) / surface_geostrophic_northward_sea_water_velocity
            'sst': 'sea_water_temperature', # sea surface temperature - cannot be used at the moment
            'temp': 'sea_water_temperature', # 3D water temperature
            'salt': 'sea_water_salinity',    # 3D salinity
            'ugrd10m': 'x_wind', # x_wind_10m
            'vgrd10m': 'y_wind', # y_wind_10m        
            'eastward_surface_stokes_drift': 'sea_surface_wave_stokes_drift_x_velocity',
            'northward_surface_stokes_drift': 'sea_surface_wave_stokes_drift_y_velocity',
            'uuss': 'sea_surface_wave_stokes_drift_x_velocity',
            'vuss': 'sea_surface_wave_stokes_drift_y_velocity',
            'hs' : 'sea_surface_wave_significant_height',
            'tp' : 'sea_surface_wave_period_at_variance_spectral_density_maximum',
            'dpm': 'sea_surface_wave_from_direction',
            'tm01' : 'sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment',
            'tm02' : 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
            'lev' : 'z',                    # This should allow correct loading of z levels - if present. TO CHECK 
                                           # lev :   vertical levels for ocean models or data
            ###########################################################################
            # Additional mapping for top-copy files on servers ########################
            }



        
        # Find all variables having standard_name - use xarray's data_vars
        self.variable_mapping = {}
        for var_name in self.Dataset.data_vars: #self.Dataset.variables:
            if var_name in [self.xname, self.yname, 'depth']:
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]
            attributes = var.attrs#attributes = var.ncattrs()
            if 'standard_name' in attributes:
                # standard_name = str(var.__dict__['standard_name'])
                standard_name = str(var.attrs['standard_name'])
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)

        self.variables = self.variable_mapping.keys()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False,
                      indrealization=None):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        if hasattr(self, 'z') and (z is not None):
            # Find z-index range
            # NB: may need to flip if self.z is ascending
            indices = np.searchsorted(-self.z, [-z.min(), -z.max()])
            indz = np.arange(np.maximum(0, indices.min() - 1 -
                                        self.verticalbuffer),
                             np.minimum(len(self.z), indices.max() + 1 +
                                        self.verticalbuffer))
            if len(indz) == 1:
                indz = indz[0]  # Extract integer to read only one layer
        else:
            indz = 0

        if indrealization == None:
            if hasattr(self, 'realizations'):
                indrealization = range(len(self.realizations))
            else:
                indrealization = None

        # Find indices corresponding to requested x and y
        if hasattr(self, 'clipped'):
            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped
        # If x or y coordinates are decreasing, we need to flip
        if hasattr(self, 'x'):
            if self.x[0] > self.x[-1]:
                indx = len(self.x) - indx
            if self.y[0] > self.y[-1]:
                indy = len(self.y) - indy
        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            if self.global_coverage():
                indx = np.arange(indx.min()-buffer, indx.max()+buffer)
            else:
                indx = np.arange(np.max([0, indx.min()-buffer]),
                                 np.min([indx.max()+buffer, self.numx]))
            indy = np.arange(np.max([0, indy.min()-buffer]),
                             np.min([indy.max()+buffer, self.numy]))
        else:
            indx[outside] = 0  # To be masked later
            indy[outside] = 0
        if indx.min() <= 0 and indx.max() >= self.numx:
            indx = np.arange(0, self.numx)

        variables = {}

        if indx.min() < 0 and indx.max() > 0:
            logging.debug('Requested data block is not continous in file'+
                          ', must read two blocks and concatenate.')
            indx_left = indx[indx<0] + self.numx  # Shift to positive indices
            indx_right = indx[indx>=0]
            continous = False
        else:
            continous = True
        for par in requested_variables:
            var = self.Dataset.variables[self.variable_mapping[par]]

            ensemble_dim = None
            if continous is True:
                if var.ndim == 2:
                    variables[par] = var[indy, indx]
                elif var.ndim == 3:
                    variables[par] = var[indxTime, indy, indx]
                elif var.ndim == 4:
                    variables[par] = var[indxTime, indz, indy, indx]
                elif var.ndim == 5:  # Ensemble data
                    variables[par] = var[indxTime, indz, indrealization, indy, indx]
                    ensemble_dim = 0  # Hardcoded ensemble dimension for now
                else:
                    raise Exception('Wrong dimension of variable: ' +
                                    self.variable_mapping[par])
            else:  # We need to read left and right parts separately
                if var.ndim == 2:
                    left = var[indy, indx_left]
                    right = var[indy, indx_right]
                    variables[par] = np.ma.concatenate((left, right), 1)
                elif var.ndim == 3:
                    left = var[indxTime, indy, indx_left]
                    right = var[indxTime, indy, indx_right]
                    variables[par] = np.ma.concatenate((left, right), 1)
                elif var.ndim == 4:
                    left = var[indxTime, indz, indy, indx_left]
                    right = var[indxTime, indz, indy, indx_right]
                    variables[par] = np.ma.concatenate((left, right), 2)
                elif var.ndim == 5:  # Ensemble data
                    left = var[indxTime, indz, indrealization,
                               indy, indx_left]
                    right = var[indxTime, indz, indrealization,
                                indy, indx_right]
                    variables[par] = np.ma.concatenate((left, right), 3)

            # If 2D array is returned due to the fancy slicing
            # methods of netcdf-python, we need to take the diagonal
            if variables[par].ndim > 1 and block is False:
                variables[par] = variables[par].diagonal()

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par],
                                         ndmin=2, mask=False)
            if block is False:
                variables[par].mask[outside] = True

            # Mask extreme values which might have slipped through
            variables[par] = np.ma.masked_outside(
                variables[par], -30000, 30000)

            # Ensemble blocks are split into lists
            if ensemble_dim is not None:
                num_ensembles = variables[par].shape[ensemble_dim]
                logging.debug('Num ensembles: %i ' % num_ensembles)
                newvar = [0]*num_ensembles
                for ensemble_num in range(num_ensembles):
                    newvar[ensemble_num] = \
                        np.take(variables[par],
                                ensemble_num, ensemble_dim)
                variables[par] = newvar

        # Store coordinates of returned points
        try:
            variables['z'] = self.z[indz]
        except:
            variables['z'] = None
        if block is True:
            if self.projected is True:
                variables['x'] = \
                    self.Dataset.variables[self.xname][indx]*self.unitfactor
                variables['y'] = \
                    self.Dataset.variables[self.yname][indy]*self.unitfactor
            else:
                variables['x'] = indx
                variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y
        if self.global_coverage():
            if self.xmax + self.delta_x >= 360 and self.xmin > 180:
                variables['x'][variables['x']>180] -= 360
        
        variables['time'] = nearestTime

        return variables
