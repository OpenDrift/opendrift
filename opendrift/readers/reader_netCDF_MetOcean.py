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

import logging

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from bisect import bisect_left

from basereader import BaseReader

# This is a reader class for the netcdf file used internally at MetOcean
# It is adapted from reader_netCDF_CF_Generic.py (found in same folder)
# using a new variable name "mapping" matrix, and an option to specify which variables 
# from the file should be used. This is expected to be used mostly when more than one (u,v) pairs are available,
# to avoid (re)mapping differents u or v components to the same 'x_sea_water_velocity' or 'y_sea_water_velocity' 
# in OpenDrift, or when we want to read only a single variable from a file (e.g. depth).
# S.Weppe  

class Reader(BaseReader):

    def __init__(self, filename=None, name=None, variables_to_use = None, **kwargs):
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
                logging.info('Opening files with MFDataset')
                self.Dataset = MFDataset(filename)
            else:
                logging.info('Opening file with Dataset')
                self.Dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

        logging.debug('Finding coordinate variables.')
        # Find x, y and z coordinates

        for var_name in self.Dataset.variables:
            logging.debug('Parsing variable: ' +  var_name)
            var = self.Dataset.variables[var_name]
            #if var.ndim > 1:
            #    continue  # Coordinates must be 1D-array
            attributes = var.ncattrs()
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
                standard_name = var.__dict__['standard_name']
            if 'long_name' in attributes:
                long_name = var.__dict__['long_name']
            if 'axis' in attributes:
                axis = var.__dict__['axis']
            if 'units' in attributes:
                units = var.__dict__['units']
            if '_CoordinateAxisType' in attributes:
                CoordinateAxisType = var.__dict__['_CoordinateAxisType']
            if standard_name == 'longitude' or \
                    long_name == 'longitude' or standard_name == 'lon' or var_name == 'lon':
                self.lon = var
                lon_var_name = var_name
            if standard_name == 'latitude' or \
                    long_name == 'latitude' or standard_name == 'lat' or var_name == 'lat':
                self.lat = var
                lat_var_name = var_name
            if axis == 'X' or \
                    CoordinateAxisType == 'Lon' or \
                    standard_name == 'projection_x_coordinate':
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                x = var[:]*unitfactor
                self.numx = var.shape[0] 
            if axis == 'Y' or \
                    CoordinateAxisType == 'Lat' or \
                    standard_name == 'projection_y_coordinate':
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                self.unitfactor = unitfactor
                y = var[:]*unitfactor
                self.numy = var.shape[0] 
            if standard_name == 'depth' or axis == 'Z' or var_name == 'lev': # probably refers to vertical levels ?
                if var[:].ndim == 1:
                    if 'positive' not in attributes or \
                            var.__dict__['positive'] == 'up':
                        self.z = var[:]
                    else:
                        self.z = -var[:]
            if standard_name == 'time' or axis == 'T' or var_name in ['time', 'vtime']:
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
        self.xmin, self.xmax = x.min(), x.max()
        self.ymin, self.ymax = y.min(), y.max()
        self.delta_x = np.abs(x[1] - x[0])
        self.delta_y = np.abs(y[1] - y[0])
        rel_delta_x = (x[1::] - x[0:-1])
        rel_delta_x = np.abs((rel_delta_x.max() -
                              rel_delta_x.min())/self.delta_x)
        rel_delta_y = (y[1::] - y[0:-1])
        rel_delta_y = np.abs((rel_delta_y.max() -
                              rel_delta_y.min())/self.delta_y)
        if rel_delta_x > 0.05:  # Allow 5 % deviation
            print rel_delta_x
            print x[1::] - x[0:-1]
            raise ValueError('delta_x is not constant!')
        if rel_delta_y > 0.05:
            print rel_delta_y
            print y[1::] - y[0:-1]
            raise ValueError('delta_y is not constant!')
        self.x = x  # Store coordinate vectors
        self.y = y

        if not hasattr(self, 'proj4'):
            if self.lon.ndim == 1:
                logging.debug('Lon and lat are 1D arrays, assuming latlong projection')
                self.proj4 = '+proj=latlong'
            elif self.lon.ndim == 2:
                logging.debug('Reading lon lat 2D arrays, since projection is not given')
                self.lon = self.lon[:]
                self.lat = self.lat[:]
                self.projected = False

        
        # Map MetOcean variable names to CF standard_name used in OpenDrift
        # 
        # see infos on variables names here : https://wiki.metocean.co.nz/display/OPS/Parameter+definitions
        
        # change the `variable_aliases` variable to fit with MetOcean convention (note the mapping was initially defined in basereader.py)
        self.variable_aliases = {
            # 'mask_rho': 'land_binary_mask',
            # 'mask_psi': 'land_binary_mask',
            'dep': 'sea_floor_depth_below_sea_level',
            'el': 'sea_surface_height',   # el = et + ssh
            'et': 'sea_surface_height',   # tidal elevation
            'ssh': 'sea_surface_height',  # non-tidal
            'um': 'x_sea_water_velocity', # depth-averaged total
            'vm': 'y_sea_water_velocity', # depth-averaged total
            'us': 'x_sea_water_velocity',  # surface current total 
            'vs': 'y_sea_water_velocity',  # surface current total
            'umo': 'x_sea_water_velocity',# depth-averaged non-tidal
            'vmo': 'y_sea_water_velocity',# depth-averaged non-tidal
            'ut': 'x_sea_water_velocity', # tide eastward_tidal_current
            'vt': 'y_sea_water_velocity', # tide northward_tidal_current
            'uo': 'x_sea_water_velocity', # non-tidal current 3D / eastward_geostrophic_current_velocity
            'vo': 'y_sea_water_velocity', # non-tidal current 3D / northward_geostrophic_current_velocity
            'u': 'x_sea_water_velocity',  # total current 3D
            'v': 'y_sea_water_velocity',  # total current 3D
            'w': 'upward_sea_water_velocity',
            'uso': 'x_sea_water_velocity', #surface layer non-tidal current (5 m below sea level) / surface_geostrophic_eastward_sea_water_velocity
            'vso': 'y_sea_water_velocity', #surface layer non-tidal current (5 m below sea level) / surface_geostrophic_northward_sea_water_velocity
            'sst': 'sea_water_temperature',
            'salt': 'c',
            'ugrd10m': 'x_wind', # x_wind_10m
            'vgrd10m': 'y_wind', # y_wind_10m
            'lev' : 'z'                    # This should allow correct loading of z levels - if present. TO CHECK 
                                           # lev :   vertical levels for ocean models or data
            }

        self.variable_mapping = {} 
        # variable_mapping provides the names of the variable to be used for each of the opendrift-convention variable names
        # this is the "opposite" of the self.variable_aliases above
        # i.e. if you need 'sea_surface_height' within opendrift, use variable 'el' in the netcdf file

        for var_name in self.Dataset.variables:
            
            if var_name in [self.xname, self.yname,'lev'] :
                 continue  # Skip coordinate variables 
            # if var_name in [self.xname, self.yname] 'depth']: //  not sure why 'dep' is skipped in original reader_netCDF_CF_generic.py - or does this refer to levels, rather than depth?
            # this is needed for sea_floor_depth_below_sea_level
                
            if variables_to_use is not None :
            #filter which variables to use    
                if var_name not in variables_to_use and var_name not in ['time','z']:
                    # skip variables not in variables_to_use, unless they are 'time' or 'z' which are needed anyways
                    logging.debug('Skipping variable : %s ; not used' % (var_name))
                    continue

            var = self.Dataset.variables[var_name]
            attributes = var.ncattrs()
            
            if 'standard_name' in attributes:
                standard_name = str(var.__dict__['standard_name'])
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)
            else: # MetOcean dataset may not have a standard_name so use var_name instead to "map" variables names to opendrift convention
                if var_name in self.variable_aliases:  # Mapping if needed, not this will be true only if var_name is found in the self.variable_aliases first column
                    standard_name = self.variable_aliases[var_name]
                self.variable_mapping[standard_name] = str(var_name) 
        
        # For now we can't differenciate the different kinds of currents (i.e. tide, resiudal, total)...as it seems opendrift only expect (u,v) currents
        # all different currents are therefore mapped to x_sea_water_velocity, y_water_velocity 
        # which means that some overwriting can occur if there are more than 2 (u,v) pairs in the file - an option may be to allow specifying which variables to use 
        # e.g. see here : https://github.com/OpenDrift/opendrift/blob/master/opendrift/readers/basereader.py line 89
        
        self.variables = self.variable_mapping.keys() # check that it does the right thing here

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
            indx[outside] = 0  # To be masked later
            indy[outside] = 0

        variables = {}

        for par in requested_variables:
            var = self.Dataset.variables[self.variable_mapping[par]]

            ensemble_dim = None
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

    def nearest_time(self, time):
        """Return nearest times before and after the requested time.

        Returns:
            nearest_time: datetime
            time_before: datetime
            time_after: datetime
            indx_nearest: int
            indx_before: int
            indx_after: int

        The function is copied from basereader.py, and is only modified to
        allow for field that are set as "always_valid" such as bathys

        """
        if self.start_time == self.end_time:
            return self.start_time, self.start_time, None, 0, 0, 0
        if self.start_time is None:
            return None, None, None, None, None, None
        # field was set as always valid
        # return the first timestep regardless of actual time queried 
        # as the actual time in original file may not fit with the simulation time
        if self.always_valid:
            # return self.start_time, self.start_time, None, 0, 0, 0
            return None,None,None, 0, 0, 0 
            # trick it by settting start_time to "time", this prevents any time interpolation in basereader.py line 470

        if hasattr(self, 'times'):  # Time as array, possibly with holes
            indx_before = np.max((0, bisect_left(self.times, time) - 1))
            if self.times[indx_before + 1] == time:
                # Correction needed when requested time exists in times
                indx_before = indx_before + 1
            time_before = self.times[indx_before]
            indx_after = np.minimum(indx_before + 1,
                                    len(self.times) - 1)  # At the end
            time_after = self.times[indx_after]
            if (time - time_before) < (time_after - time):
                indx_nearest = indx_before
            else:
                indx_nearest = indx_after
            nearest_time = self.times[indx_nearest]
        else:  # Time step is constant (no holes)
            indx = float((time - self.start_time).total_seconds()) / \
                float(self.time_step.total_seconds())
            indx_nearest = int(round(indx))
            nearest_time = self.start_time + indx_nearest*self.time_step
            indx_before = int(np.floor(indx))
            time_before = self.start_time + indx_before*self.time_step
            indx_after = int(np.ceil(indx))
            time_after = self.start_time + indx_after*self.time_step
        return nearest_time, time_before, time_after,\
            indx_nearest, indx_before, indx_after