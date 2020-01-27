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

from opendrift.readers.basereader import BaseReader
from opendrift.readers.interpolation import ReaderBlock
import sys

# This is a reader class for the netcdf file used internally at MetOcean
# It is adapted from reader_netCDF_CF_Generic.py (found in same folder)
# using a new variable name "mapping" matrix, and an option to specify which variables 
# from the file should be used. This is expected to be used mostly when more than one (u,v) pairs are available,
# to avoid (re)mapping differents u or v components to the same 'x_sea_water_velocity' or 'y_sea_water_velocity' 
# in OpenDrift, or when we want to read only a single variable from a file (e.g. depth).
#
# For reference, there is actually some option to select variables in model.add_reader(variables=['var1']), but this is not exactly doing what 
# we need hereto_longitude_0_360
#
# S.Weppe  

# needed to be copied from basereader.py
# Some valid (but extreme) ranges for checking that values are reasonable
standard_names = {
    'x_wind': {'valid_min': -50, 'valid_max': 50},
    'y_wind': {'valid_min': -50, 'valid_max': 50},
    'x_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'y_sea_water_velocity': {'valid_min': -10, 'valid_max': 10}}

# Identify x-y vector components/pairs for rotation (NB: not east-west pairs!)
vector_pairs_xy = [
    ['x_wind', 'y_wind'],
    ['x_sea_water_velocity', 'y_sea_water_velocity'],
    ['sea_surface_wave_stokes_drift_x_velocity',
     'sea_surface_wave_stokes_drift_y_velocity']
    ]


class Reader(BaseReader):

    def __init__(self, filename=None, name=None, variables_to_use = None, **kwargs):

        self.use_log_profile = False # No extrapolation using logarithmic profile by default
        self.reader_water_depth = None # optional "water depth reader" passed to reader

        # allow direct initialization of class attribute using kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)
        #------------------------------------------------------------
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
                self.Dataset = MFDataset(filename, aggdim = 'time') # added aggdim = 'time'
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
                # adding a flag to inform on longitude convention [0,360] or [-180,180]
                # used further below on covers_positions and longitude index selection
                # search "has_lon_0_360" if needed

                if  (self.lon[:]>=0).all() and (self.lon[:]<=360).all() : 
                # reader longitude using convention :  0<lon<360, or has longitude between 0 and +180 only
                    self.has_lon_0_360 = True
                else:
                    self.has_lon_0_360 = False    

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
                    # levels in MetOcean file convention are 0 at the surface, positive down i.e lev = 5 means 5 m below surface
                    # need to be converted to negative down to fit with opendrift convention
                    self.z = -var[:] 
                    # original code-----------------------  
                    # if 'positive' not in attributes or \
                    #         var.__dict__['positive'] == 'up':
                    #     self.z = var[:]
                    # else:
                    #     self.z = -var[:] 
                    # --------------------------------- 
            if standard_name == 'time' or axis == 'T' or var_name in ['time', 'vtime']:
                # Read and store time coverage (of this particular file)
                time = var[:] 
                # time will be incorrect if we read our top copy files that have inconsistent time reference - works good with UDS-queried files though
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
            # elif self.lon.ndim == 2: # added
            #     x = self.lon[:]
            #     self.xname = lon_var_name
            #     self.numx = x.shape[0]#len(x)
            else:
                raise ValueError('Did not find x-coordinate variable')
        if 'y' not in locals():
            if self.lat.ndim == 1:
                y = self.lat[:]
                self.yname = lat_var_name
                self.numy = len(y)
            # elif self.lat.ndim == 2: # added
            #     y = self.lat[:]
            #     self.yname = lat_var_name
            #     self.numy = y.shape[1]#len(x)
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
            print(rel_delta_x)
            print(x[1::] - x[0:-1])
            raise ValueError('delta_x is not constant!')
        if rel_delta_y > 0.05:
            print(rel_delta_y)
            print(y[1::] - y[0:-1])
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
        # http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
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
            'um': 'x_sea_water_velocity', # depth-averaged total / eastward_sea_water_velocity
            'vm': 'y_sea_water_velocity', # depth-averaged total / northward_sea_water_velocity
            'us': 'x_sea_water_velocity', # surface current total / surface_eastward_sea_water_velocity
            'vs': 'y_sea_water_velocity', # surface current total / surface_northward_sea_water_velocity
            'uso': 'x_sea_water_velocity', # surface current non-tidal / surface_eastward_sea_water_velocity_assuming_no_tide
            'vso': 'y_sea_water_velocity', # surface current non-tidal  / surface_northward_sea_water_velocity_assuming_no_tide
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
            'northward_wind': 'x_wind', # x_wind_10m
            'eastward_wind': 'y_wind', # y_wind_10m            
            'eastward_surface_stokes_drift': 'sea_surface_wave_stokes_drift_x_velocity',
            'northward_surface_stokes_drift': 'sea_surface_wave_stokes_drift_y_velocity',
            'uuss': 'sea_surface_wave_stokes_drift_x_velocity',
            'vuss': 'sea_surface_wave_stokes_drift_y_velocity',
            'hs' : 'sea_surface_wave_significant_height',
            'tp' : 'sea_surface_wave_period_at_variance_spectral_density_maximum',
            'dpm': 'sea_surface_wave_from_direction',
            'tm01' : 'sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment',
            'tm02' : 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
            'lev' : 'z'                    # This should allow correct loading of z levels - if present. TO CHECK 
                                           # lev :   vertical levels for ocean models or data
            }

        self.variable_mapping = {} 
        # variable_mapping provides the names of the variable to be used for each of the opendrift-convention variable names
        # this is the "opposite" of the self.variable_aliases above
        # i.e. if you need 'sea_surface_height' within opendrift, use variable 'el' in the netcdf file
        
        for var_name in self.Dataset.variables:

            if var_name in [self.xname, self.yname,'lev','time'] :
                 continue  # Skip coordinate variables 
            # if var_name in [self.xname, self.yname] 'depth']: //  not sure why 'dep' is skipped in original reader_netCDF_CF_generic.py, probably refers to levels, rather than depth?
            # 'depth' or 'dep' should not be skipepd to allow mapping 'sea_floor_depth_below_sea_level'
                
            if variables_to_use is not None :
            #filter which variables to use    
                if var_name not in variables_to_use and var_name not in ['time','z']:
                    # skip variables not in variables_to_use, unless they are 'time' or 'z' which are needed anyways
                    logging.debug('Skipping variable : %s ; not used' % (var_name))
                    continue

            var = self.Dataset.variables[var_name]
            attributes = var.ncattrs()

            if False :#'standard_name' in attributes:
                # for MetOcean datasets force the use of the 'var_name' to map the alias
                # instead of standard_name - see below
                standard_name = str(var.__dict__['standard_name'])
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)
            else: # MetOcean dataset may not have a standard_name so use var_name instead to "map" variables names to opendrift convention
                if var_name in self.variable_aliases:  # Mapping if needed, not this will be true only if var_name is found in the self.variable_aliases first column
                    standard_name = self.variable_aliases[var_name]
                self.variable_mapping[standard_name] = str(var_name) 
            
        # For now we can't differenciate the different kinds of currents (i.e. tide, residual, total)...as it seems opendrift only expect general (u,v) currents
        # All different currents are therefore mapped to x_sea_water_velocity, y_water_velocity for now
        # Here we have added an option variables_to_use to allow specifying which variables to use in a file
        # example : reader_roms_cnz_depth = reader_netCDF_MetOcean.Reader('C:\metocean\cnz19800801_00z_surf.nc',variables_to_use = ['dep','um','vm']) 
        # other variables will be discarded  
        
        self.variables = self.variable_mapping.keys() # check that it does the right thing here
 
        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def __add__(self,reader_to_add):
        # function to overload the standard "add" (or '+') function
        # see https://github.com/OpenDrift/opendrift/issues/73
        # 
        # this will allow doing something this :
        # reader_cur_tot = reader_cur_tide + reader_cur_res
        # where reader_cur_tot will yields the sum of the currents returned by 
        # reader_cur_tide and reader_cur_res 
        # 
        # Only works for current velocities for now :
        # ['x_sea_water_velocity','y_sea_water_velocity']
        # 
        # WARNING : For now the "readers_to_add" fields is added to the first reader
        # input in the addition. If the same reader object is re-used (without deleting it before)
        # then it will continue to add readers (even if called by its initial name e.g. reader_cur_tide)
        # 
        # 
        # can easily be extended for other additions where relevant
        if 'x_sea_water_velocity' in self.variables and \
           'x_sea_water_velocity' in reader_to_add.variables and \
           'y_sea_water_velocity' in self.variables and \
           'y_sea_water_velocity' in reader_to_add.variables :

            logging.debug('Adding currents from readers : %s and %s' % (self.name,reader_to_add.name) )
            if not hasattr(self,'readers_to_add'):
                self.readers_to_add = [] # pre-allocate
                self.use_multireader_interp = True # switch used later in get_variables_interpolated()

            self.readers_to_add.append(reader_to_add) # add field to existing reader (self)
        
        # this is probably not a real "add" overloading in pythonic sense
        # since we simply add a new attribute "readers_to_add" to self containing the other readers to add
        # These can be accessed as self.readers_to_add[0],self.readers_to_add[1] etc...
        # and will be used later in the get_variables_interpolated function
        #
        return self


    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False,
                      indrealization=None):
        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)
            
        if hasattr(self, 'z') and (z is not None): # 3D interpolation

            # "self.z"  : contains the vertical levels of the dataset, surface=0, negative down
            #           : expected to be ordered as surface-most to deeper-most
            #      
            # "z"  : are the particle depths - z=0 at surface, z<0 if below water
            # 
            
            # Find z-index ranges that includes all particles
            # NB: may need to flip if self.z is ascending (original comment)
            # 
            # np.searchsorted(a,v) finds indices where elements V should be inserted in A to maintain order, assuming A is sorted in ascending order.
            # 
            if self.z[0]>self.z[1] : # self.z is "descending" - surface-most layer first then going deeper
                indices = np.searchsorted(-self.z, [-z.min(), -z.max()]) 
            else : # self.z is "ascending" - bottom-most layer first then going towards surface
                indices = np.searchsorted(self.z, [z.min(), z.max()])

            indz = np.arange(np.maximum(0, indices.min() - 1 - self.verticalbuffer),
                             np.minimum(len(self.z), indices.max() + 1 + self.verticalbuffer))
            if len(indz) == 1:
                # in that case indz is single closest layer to particle z
                indz = indz[0]  # Extract integer to read only one layer
                
        else: # 2D interpolation
            indz = 0

        if indrealization == None:
            if hasattr(self, 'realizations'):
                indrealization = range(len(self.realizations))
            else:
                indrealization = None
        ##########################################################################################################        
        # Find indices corresponding to requested x and y
        # indx = np.floor((x-self.xmin)/self.delta_x).astype(int) # original code
        if self.has_lon_0_360 :
            indx = np.floor((self.to_longitude_0_360(x)-self.xmin)/self.delta_x).astype(int) 
        else:
            indx = np.floor((x-self.xmin)/self.delta_x).astype(int) 
       
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int)
        ########################################################################################################## 

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
            # convert self.z to non masked array if needed
            if isinstance(self.z, np.ma.MaskedArray):
                variables['z'] = self.z[indz].data
            else:
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

 
    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader.
           
           Overloads covers_positions() from basereader.py
           to account for netcdf files where 0<longitude<360

        """
        # Calculate x,y coordinates from lon,lat
        x, y = self.lonlat2xy(lon, lat)
        # Only checking vertical coverage if zmin, zmax is defined
        zmin = -np.inf
        zmax = np.inf
        if hasattr(self, 'zmin') and self.zmin is not None:
            zmin = self.zmin
        if hasattr(self, 'zmax') and self.zmax is not None:
            zmax = self.zmax

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        else:
            # at this stage 'lon' follows the convention -180<lon<180
            # need to account for reader(s) that may use 0<lon<360 instead
            if not self.has_lon_0_360 :
            
                indices = np.where((x >= self.xmin) & (x <= self.xmax) &
                               (y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
            elif self.has_lon_0_360: #reader longitude using convention :  0<lon<360, or has longitude between 0 and +180 only
                x360 = self.to_longitude_0_360(x) # convert opendrift x ([-180,180]) to [0,360] convention
                indices = np.where((x360 >= self.xmin) & (x360 <= self.xmax) &
                                   (y >= self.ymin) & (y <= self.ymax) &
                                   (z >= zmin) & (z <= zmax))[0]
        try:
            return indices, x[indices], y[indices]
        except:
            return indices, x, y   


    # this will overload the get_variables_interpolated() from basereader.py
    # to allow applying log profile correction, adding readers etc..
    # eventually this may be a nice feature to add the basereader ?
    def get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):        
        self.timer_start('total')
        self.timer_start('preparing')
        
        # convert lon,lat positions to fit with reader conventions i.e. 0<lon<360 or -180<lon<180
        # so that interpolation is handled correctly
        # opendrit internally uses -180<lon<180 i.e. at this stage -180<lon<180
        if self.has_lon_0_360:
            lon = self.to_longitude_0_360(lon) 

        #----------------------------------------------------------------------------------------------------
        # Bypass the function to get to get_variables_interpolated_multireaders() if there are readers to add
        if hasattr(self,'readers_to_add') and self.use_multireader_interp and 'x_sea_water_velocity' in variables: 
            # "special" case when we need to add currents ['x_sea_water_velocity','y_sea_water_velocity'] from two or more readers
            # by passed otherwise
            # We use the function get_variables_interpolated_multireaders() which calls 
            # get_variables_interpolated() for the two or more readers to add.
            env,env_profiles = self.get_variables_interpolated_multireaders(variables, profiles=profiles, \
                                   profiles_depth=profiles_depth, time=time, \
                                   lon=lon, lat=lat, z=z, \
                                   block=block, rotate_to_proj=rotate_to_proj)
            # import pdb;pdb.set_trace()
            return env, env_profiles
        #----------------------------------------------------------------------------------------------------

        # Raise error if time not not within coverage of reader
        if not self.covers_time(time):
            raise ValueError('%s is outside time coverage (%s - %s) of %s' %
                             (time, self.start_time, self.end_time, self.name))

       # Check which particles are covered (indep of time) - align with new code in basereader.py
        ind_covered, reader_x, reader_y = self.covers_positions(lon, lat, z)

        if len(ind_covered) == 0:
            raise ValueError(('All %s particles (%.2f-%.2fE, %.2f-%.2fN) ' +
                              'are outside domain of %s (%s)') %
                             (len(lon), lon.min(), lon.max(), lat.min(),
                              lat.max(), self.name, self.coverage_string()))

        # Find reader time_before/time_after
        time_nearest, time_before, time_after, i1, i2, i3 = \
            self.nearest_time(time)
        logging.debug('Reader time:\n\t\t%s (before)\n\t\t%s (after)' %
                      (time_before, time_after))
        if time == time_before:
            time_after = None
        
        z = z.copy()[ind_covered]  # Send values and not reference
                                   # to avoid modifications

        if block is False or self.return_block is False:
            # Analytical reader, continous in space and time
            self.timer_end('preparing')
            env_before = self._get_variables(variables, profiles,
                                             profiles_depth,
                                             time,
                                             #time_before,
                                             reader_x, reader_y, z,
                                             block=block)
            logging.debug('Fetched env-before')
            self.timer_start('preparing')

        else:
            # Swap before- and after-blocks if matching times
            if str(variables) in self.var_block_before:
                block_before_time = self.var_block_before[
                    str(variables)].time
                if str(variables) in self.var_block_after:
                    block_after_time = self.var_block_after[
                        str(variables)].time
                    if block_before_time != time_before:
                        if block_after_time == time_before:
                            self.var_block_before[str(variables)] = \
                                self.var_block_after[str(variables)]
                    if block_after_time != time_after:
                        if block_before_time == time_before:
                            self.var_block_after[str(variables)] = \
                                self.var_block_before[str(variables)]
            # Fetch data, if no buffer is available
            if (not str(variables) in self.var_block_before) or \
                    (self.var_block_before[str(variables)].time !=
                     time_before):
                self.timer_end('preparing')
                reader_data_dict = \
                    self._get_variables(variables, profiles,
                                        profiles_depth, time_before,
                                        reader_x, reader_y, z,
                                        block=block)
                self.timer_start('preparing')
                self.var_block_before[str(variables)] = \
                    ReaderBlock(reader_data_dict,
                                interpolation_horizontal=self.interpolation)
                try:
                    len_z = len(self.var_block_before[str(variables)].z)
                except:
                    len_z = 1
                logging.debug(('Fetched env-block (size %ix%ix%i) ' +
                              'for time before (%s)') %
                              (len(self.var_block_before[str(variables)].x),
                               len(self.var_block_before[str(variables)].y),
                               len_z, time_before))
            if not str(variables) in self.var_block_after or \
                    self.var_block_after[str(variables)].time != time_after:
                if time_after is None:
                    self.var_block_after[str(variables)] = \
                        self.var_block_before[str(variables)]
                else:
                    self.timer_end('preparing')
                    reader_data_dict = \
                        self._get_variables(variables, profiles,
                                            profiles_depth, time_after,
                                            reader_x, reader_y, z,
                                            block=block)
                    self.timer_start('preparing')
                    self.var_block_after[str(variables)] = \
                        ReaderBlock(
                            reader_data_dict,
                            interpolation_horizontal=self.interpolation)
                    try:
                        len_z = len(self.var_block_after[str(variables)].z)
                    except:
                        len_z = 1

                    logging.debug(('Fetched env-block (size %ix%ix%i) ' +
                                  'for time after (%s)') %
                                  (len(self.var_block_after[
                                       str(variables)].x),
                                   len(self.var_block_after[
                                       str(variables)].y),
                                   len_z, time_after))

            if self.var_block_before[str(variables)].covers_positions(
                reader_x, reader_y) is False or \
                self.var_block_after[str(variables)].covers_positions(
                    reader_x, reader_y) is False:
                logging.warning('Data block from %s not large enough to '
                                'cover element positions within timestep. '
                                'Buffer size (%s) must be increased.' %
                                (self.name, str(self.buffer)))
                # import pdb;pdb.set_trace()
                # note the covers_positions() of a "reader block" object may not work correctly
                # 0<lon<360 i.e. return false while it is should be true
                # Not critical - simply ensure buffer is larhe enough by setting 
                # o.max_speed = 5 m/s for example

            self.timer_end('preparing')
            ############################################################
            # Interpolate before/after blocks onto particles in space
            ############################################################
            self.timer_start('interpolation')
            logging.debug('Interpolating before (%s) in space  (%s)' %
                          (self.var_block_before[str(variables)].time,
                           self.interpolation))
            env_before, env_profiles_before = self.var_block_before[
                str(variables)].interpolate(
                    reader_x, reader_y, z, variables,
                    profiles, profiles_depth)

            if (time_after is not None) and (time_before != time):
                logging.debug('Interpolating after (%s) in space  (%s)' %
                              (self.var_block_after[str(variables)].time,
                               self.interpolation))
                env_after, env_profiles_after = self.var_block_after[
                    str(variables)].interpolate(
                        reader_x, reader_y, z, variables,
                        profiles, profiles_depth)

            self.timer_end('interpolation')

        #######################
        # Time interpolation
        #######################
        self.timer_start('interpolation_time')
        env_profiles = None
        if (time_after is not None) and (time_before != time):
            weight_after = ((time - time_before).total_seconds() /
                            (time_after - time_before).total_seconds())
            logging.debug(('Interpolating before (%s, weight %.2f) and'
                           '\n\t\t      after (%s, weight %.2f) in time') %
                          (self.var_block_before[str(variables)].time,
                           1 - weight_after,
                           self.var_block_after[str(variables)].time,
                           weight_after))
            env = {}
            for var in variables:
                # Weighting together, and masking invalid entries
                env[var] = np.ma.masked_invalid((env_before[var] *
                                                (1 - weight_after) +
                                                env_after[var] * weight_after))

                if var in standard_names.keys():
                    invalid = np.where((env[var] < standard_names[var]['valid_min']) 
                               | (env[var] > standard_names[var]['valid_max']))[0]
                    if len(invalid) > 0:
                        logging.warning('Invalid values found for ' + var)
                        logging.warning(env[var][invalid])
                        logging.warning('(allowed range: [%s, %s])' %
                                        (standard_names[var]['valid_min'],
                                         standard_names[var]['valid_max']))
                        logging.warning('Replacing with NaN')
                        env[var][invalid] = np.nan
            # Interpolating vertical profiles in time
            if profiles is not None:
                env_profiles = {}
                logging.info('Interpolating profiles in time')
                # Truncating layers not present both before and after
                numlayers = np.minimum(len(env_profiles_before['z']),
                                       len(env_profiles_after['z']))
                env_profiles['z'] = env_profiles_before['z'][0:numlayers+1]
                for var in env_profiles_before.keys():
                    if var == 'z':
                        continue
                    env_profiles[var] = (
                        env_profiles_before[var][0:numlayers, :] *
                        (1 - weight_after) +
                        env_profiles_after[var][0:numlayers, :]*weight_after)
            else:
                env_profiles = None

        else:
            logging.debug('No time interpolation needed - right on time.')
            env = env_before
            if profiles is not None:
                if 'env_profiles_before' in locals():
                    env_profiles = env_profiles_before
                else:
                    # Copying data from environment to vertical profiles
                    env_profiles = {'z': profiles_depth}
                    for var in profiles:
                        env_profiles[var] = np.ma.array([env[var], env[var]])
        self.timer_end('interpolation_time')

        ####################
        # Rotate vectors
        ####################
        if rotate_to_proj is not None:
            if self.simulation_SRS is True:
                logging.debug('Reader SRS is the same as calculation SRS - '
                              'rotation of vectors is not needed.')
            else:
                vector_pairs = []
                for var in variables:
                    for vector_pair in vector_pairs_xy:
                        if var in vector_pair:
                            counterpart = list(set(vector_pair) -
                                               set([var]))[0]
                            if counterpart in variables:
                                vector_pairs.append(vector_pair)
                            else:
                                sys.exit('Missing component of vector pair:' +
                                         counterpart)
                # Extract unique vector pairs
                vector_pairs = [list(x) for x in set(tuple(x)
                                for x in vector_pairs)]

                if len(vector_pairs) > 0:
                    self.timer_start('rotating vectors')
                    for vector_pair in vector_pairs:
                        env[vector_pair[0]], env[vector_pair[1]] = \
                            self.rotate_vectors(reader_x, reader_y,
                                                env[vector_pair[0]],
                                                env[vector_pair[1]],
                                                self.proj, rotate_to_proj)
                        if profiles is not None and vector_pair[0] in profiles:
                            sys.exit('Rotating profiles of vectors '
                                     'is not yet implemented')
                    self.timer_end('rotating vectors')

        # Masking non-covered pixels
        self.timer_start('masking')
        if len(ind_covered) != len(lon):
            logging.debug('Masking %i elements outside coverage' %
                          (len(lon)-len(ind_covered)))
            for var in variables:
                tmp = np.nan*np.ones(lon.shape)
                tmp[ind_covered] = env[var].copy()
                env[var] = np.ma.masked_invalid(tmp)
                # Filling also fin missing columns
                # for env_profiles outside coverage
                if env_profiles is not None and var in env_profiles.keys():
                    tmp = np.nan*np.ones((env_profiles[var].shape[0],
                                          len(lon)))
                    tmp[:, ind_covered] = env_profiles[var].copy()
                    env_profiles[var] = np.ma.masked_invalid(tmp)

        self.timer_end('masking')

        ##################################
        # Apply log profile, if applicable
        ##################################
        self.timer_start('log_profile')
        # Need to build a "readerblock" for 'sea_floor_depth_below_sea_level'
        # to get total water depth a particle locations
        # **this is needed only if we use a log profile, and only once (assuming depth doesnt change over time for now)
        if self.use_log_profile and 'x_sea_water_velocity' in variables:
            # build water depth interpolator - only once
            if not hasattr(self,'water_depth_interp'):
                logging.debug('Building ''sea_floor_depth_below_sea_level'' interpolator for logarithmic profile extrapolation')
                # water depth is available within file
                if 'sea_floor_depth_below_sea_level' in self.variables : 
                    water_depth = self._get_variables(['sea_floor_depth_below_sea_level'], profiles,profiles_depth,time,reader_x, reader_y, z,block=block)
                    self.water_depth_interp = ReaderBlock(water_depth,interpolation_horizontal=self.interpolation)
                    self.water_depth_interp.z = None
                    del water_depth
                # 'water_depth_reader' was passed when defining the reader
                elif hasattr(self,'reader_water_depth'): 
                    water_depth = self.reader_water_depth._get_variables(['sea_floor_depth_below_sea_level'], profiles,profiles_depth,time,reader_x, reader_y, z,block=block)
                    self.water_depth_interp = ReaderBlock(water_depth,interpolation_horizontal=self.interpolation)
                    self.water_depth_interp.z = None
                    del water_depth  
                # no water depth info available - cannot apply log profile             
                else:
                    logging.debug('required variable ''sea_floor_depth_below_sea_level'' not available in reader %s' % (self.name))
                    logging.error('can not apply logarithmic profile')
                    logging.error('set use_log_profile = False, add depth variable to netcdf files, or pass reader_water_depth')
                    import pdb;pdb.set_trace()
                    pass

            # get total water depth at particle positions (lon,lat) covered by reader
            water_depth_at_part,tmp = self.water_depth_interp.interpolate(lon,lat,variables = ['sea_floor_depth_below_sea_level'])
            
            del tmp
            log_fac = np.array([1.]*len(lon))# allocate
            # compute log factor value for particles covered by reader
            log_fac_covered = self.logarithmic_current_profile(z, water_depth_at_part['sea_floor_depth_below_sea_level'][ind_covered])
            # update log_fac array with computed values
            log_fac[ind_covered] = log_fac_covered
            logging.debug('Applying logarithmic profile to depth-averaged currents from reader %s' %(self.name) 
                + ' to %s particles within its coverage' %len(ind_covered) )
            logging.debug('\t\t%s   <- log_ratios[-] ->   %s' % (np.nanmin(log_fac), np.nanmax(log_fac)))
            env['x_sea_water_velocity'] = env['x_sea_water_velocity'] * log_fac
            env['y_sea_water_velocity'] = env['y_sea_water_velocity'] * log_fac
                
        self.timer_end('log_profile')
        #####################################

        self.timer_end('total')

        ###########################
        # Return interpolated data 
        ###########################
        return env, env_profiles
        
    def get_variables_interpolated_multireaders(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):
        # function to loop through the different readers and add interpolated variables
        # used when a reader is made up of an addition of readers i.e. cur_tot = cur_tide + cur_res
        # only functional for currents for now : [x_sea_water_velocity,y_sea_water_velocity]
        
        self.use_multireader_interp = False 
        # switch to False so that the standard get_variables_interpolated() is used 
        # i.e. it doesnt go to get_variables_interpolated_multireaders()

        # first reader
        env,env_profiles = self.get_variables_interpolated(variables, profiles=profiles, \
                                   profiles_depth=profiles_depth, time=time, \
                                   lon=lon, lat=lat, z=z, \
                                   block=block, rotate_to_proj=rotate_to_proj)

        # readers to add
        for reader in self.readers_to_add :
            # get interpolated data from that reader
            reader.use_multireader_interp = False 
            env_add,env_profiles_add = reader.get_variables_interpolated(['x_sea_water_velocity','y_sea_water_velocity'], profiles=profiles, \
                                   profiles_depth=profiles_depth, time=time, \
                                   lon=lon, lat=lat, z=z, \
                                   block=block, rotate_to_proj=rotate_to_proj)
            logging.debug('Adding currents from reader %s' % (reader.name))
            # add to existing env,env_profiles data
            env['x_sea_water_velocity'] += env_add['x_sea_water_velocity']
            env['y_sea_water_velocity'] += env_add['y_sea_water_velocity']

            if env_profiles is not None:
                env_profiles['x_sea_water_velocity'] += env_profiles_add['x_sea_water_velocity']
                env_profiles['y_sea_water_velocity'] += env_profiles_add['y_sea_water_velocity']
        
        self.use_multireader_interp = True 
        # set back to true for next interpolation round, i.e. will go to get_variables_interpolated_multireaders()

        ###########################
        # Return added interpolated data 
        ###########################
        return env, env_profiles


    def logarithmic_current_profile(self, particle_z, total_depth):
        ''' 
        Extrapolation of depth-averaged currents to any vertical 
        level of the water column assuming a logarithmic profile


        Inputs :
            particle_z : vertical position of particle in water column
            total_depth : total water depth at particle position
            z0 : roughness length, in meters (default, z0 = 0.001m )

        Returns : 
            Factors to be apply to interpolated raw depth-averaged currents

        Reference :
            Van Rijn, 1993. Principles of Sediment Transport in Rivers,
            Estuaries and Coastal Seas

        '''

        # Opendrift convention : particle_z is 0 at the surface, negative down
        # we need the height of particle above seabed (positive)
        part_z_above_seabed = np.abs(total_depth) + particle_z 
        if not hasattr(self,'z0'): 
            self.z0 = 0.001 # typical value for sandy seabed
        log_fac = ( np.log(part_z_above_seabed / self.z0) ) / ( np.log(np.abs(total_depth)/self.z0)-1 ) # total_depth must be positive, hence the abs()
        return log_fac



    def to_longitude_0_360(self,lon180):
        '''
        converts  longitude to different convention
        -180<longitude<180 toward 0<longitude<360

        '''

        lon360 = lon180.copy() # convert lon to [0-360] convention
        lon360[np.where(lon360<0)] = 360+ lon360[np.where(lon360<0)]
        return lon360