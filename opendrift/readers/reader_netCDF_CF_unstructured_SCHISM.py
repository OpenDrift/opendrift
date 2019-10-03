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

#####################################
# NOTE:
# This reader is under development,
# and presently not fully functional
#####################################

import logging

import numpy as np
from future.utils import iteritems
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree #cython-based KDtree for quick nearest-neighbor search
from scipy.spatial import ConvexHull # convex hull of grid points
from matplotlib.path import Path # convex hull of grid points

from opendrift.readers.basereader import BaseReader, pyproj

# Some valid (but extreme) ranges for checking that values are reasonable
standard_names = {
    'x_wind': {'valid_min': -50, 'valid_max': 50},
    'y_wind': {'valid_min': -50, 'valid_max': 50},
    'x_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'y_sea_water_velocity': {'valid_min': -10, 'valid_max': 10}}

# Identify x-y vector components/pairs for rotation (NB: not east-west pairs!)
vector_pairs_xy = [
    ['x_wind', 'y_wind'],
    ['sea_ice_x_velocity', 'sea_ice_y_velocity'],
    ['x_sea_water_velocity', 'y_sea_water_velocity'],
    ['sea_surface_wave_stokes_drift_x_velocity',
     'sea_surface_wave_stokes_drift_y_velocity']
    ]


class Reader(BaseReader):

    def __init__(self, filename=None, name=None, proj4=None, use_3d = None):
        """Initialise reader_netCDF_CF_unstructured_SCHISM

        Args:
            filename:  name of SCHISM netcdf file (can have wildcards)
            name: name of reader - optional, taken as filename if not input
                  o.readers['name']
            proj4: proj4 string defining spatial reference system. 
                   find string here : https://spatialreference.org/ref/epsg/
            use_3d: switch to use 3d flows (if available)
        """
        if filename is None:
            raise ValueError('Need filename as argument to constructor')
        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        # Due to misspelled standard_name in
        # some (Akvaplan-NIVA) FVCOM files
        variable_aliases = {
            'eastward_sea_water_velocity': 'x_sea_water_velocity',
            'Northward_sea_water_velocity': 'y_sea_water_velocity',
            'eastward wind': 'x_wind',
            'northward wind': 'y_wind'
            }

        # Mapping FVCOM variable names to CF standard_name
        fvcom_mapping = {
            'um': 'x_sea_water_velocity',
            'vm': 'y_sea_water_velocity'}
        
        # [name_used_in_schism : equivalent_CF_name]
        schism_mapping = {
            'dahv': 'x_sea_water_velocity',
            'dahv': 'y_sea_water_velocity',
            'hvel': 'x_sea_water_velocity',
            'hvel': 'y_sea_water_velocity',
            'depth': 'sea_floor_depth_below_sea_level',
            'elev' : 'sea_surface_height',
            'temp' : 'sea_water_temperature',
            'salt' : 'sea_water_salinity',
            'zcor' : 'zcor', # time-varying vertical coordinates
            'sigma': 'ocean_s_coordinate'}
            # diffusivity
            # viscosity
            # vertical_velocity

        self.return_block = True

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

        # Define projection of input data 
        if proj4 is not None: #  user has provided a projection apriori
            self.proj4 = proj4
        else:                 # no input assumes latlon
            logging.debug('Lon and lat are 1D arrays, assuming latlong projection: proj4 = ''+proj=latlong''')
            self.proj4 = '+proj=latlong'

        # check if 3d data is available and if we should use it
        self.use_3d = use_3d
        if self.use_3d is None : # not specified, use 3d data by default (if available)
            if 'hvel' in self.Dataset.variables:
                self.use_3d = True
            else:
                self.use_3d = False
        if self.use_3d and 'hvel' not in self.Dataset.variables:
            logging.debug('No 3D velocity data in file - cannot find variable ''hvel'' ')


        logging.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        for var_name in self.Dataset.variables:
            
            if var_name in ['SCHISM_hgrid_face_x','SCHISM_hgrid_face_y','SCHISM_hgrid_edge_x','SCHISM_hgrid_edge_y']:
                # all these variables have the same standard name projection_x_coordinate,projection_y_coordinate
                # we need to only use :
                # SCHISM_hgrid_node_x as projection_x_coordinate
                # SCHISM_hgrid_node_y as projection_y_coordinate
                # which are the nodes (also called vertices)
                continue

            var = self.Dataset.variables[var_name]
            if var.ndim > 1:
                continue  # Coordinates must be 1D-array
            attributes = var.ncattrs()
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            CoordinateAxisType = ''
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
                    long_name == 'longitude' or \
                    var_name == 'longitude' or \
                    axis == 'X' or \
                    CoordinateAxisType == 'Lon' or \
                    standard_name == 'projection_x_coordinate':
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                x = var[:]*unitfactor
                self.unitfactor = unitfactor
                self.numx = var.shape[0]
            if standard_name == 'latitude' or \
                    long_name == 'latitude' or \
                    var_name == 'latitude' or \
                    axis == 'Y' or \
                    CoordinateAxisType == 'Lat' or \
                    standard_name == 'projection_y_coordinate':
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                y = var[:]*unitfactor
                self.numy = var.shape[0]
            if standard_name == 'depth' or axis == 'Z':
                if 'positive' not in var.ncattrs() or \
                        var.__dict__['positive'] == 'up':
                    self.z = var[:]
                else:
                    self.z = -var[:]
            if standard_name == 'time' or axis == 'T' or var_name == 'time':
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

        if 'x' not in locals():
            raise ValueError('Did not find x-coordinate variable')
        if 'y' not in locals():
            raise ValueError('Did not find y-coordinate variable')

        self.lon = x
        self.lat = y

        # computing KDtree and convex hull of reader's nodes
        # This is done using cython-based cKDTree from scipy for quick nearest-neighbor search
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
        logging.debug('Building KDtree for nearest-neighbor search')
        self.KDtree = cKDTree(np.vstack((self.lon,self.lat)).T)

        # build convex hull of points for particle-in-mesh checks
        logging.debug('Building convex hull of nodes for particle''s in-mesh checks')
        # *** this a very rough way of doing it...ideally we would be able to 
        # have the actual contour of the grid....
        hull_obj = ConvexHull(np.vstack([self.lon.data,self.lat.data]).T)
        self.hull_path = Path(np.vstack([self.lon[hull_obj.vertices],self.lat[hull_obj.vertices]]).T)  # Matplotlib Path object
        self.hull = np.vstack([self.lon[hull_obj.vertices],self.lat[hull_obj.vertices]]).T


        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname, 'depth']:
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]
            attributes = var.ncattrs()
            # if 'standard_name' in attributes:
            #     standard_name = str(var.__dict__['standard_name'])
            #     if standard_name in variable_aliases:  # Mapping if needed
            #         standard_name = variable_aliases[standard_name]
            #     self.variable_mapping[standard_name] = str(var_name)
            # elif var_name in fvcom_mapping:
            #     self.variable_mapping[fvcom_mapping[var_name]] = \
            #         str(var_name)
            if var_name in schism_mapping:
                # current velocity variables ['dahv','hvel'] need special treatment because
                # [u,v] components are saved in the same variable i.e.
                # u = nc.variables['hvel'][0,:,:]
                # v = nc.variables['hvel'][1,:,:]
                # 
                # this means both x_sea_water_velocity and y_sea_water_velocity
                # will be mapped to same name. Correct data will be then extracted
                # in get_variables()
                # 
                if var_name == 'hvel' and self.use_3d : # then use 3d data
                    self.variable_mapping['x_sea_water_velocity'] = str(var_name)
                    self.variable_mapping['y_sea_water_velocity'] = str(var_name)   
                elif var_name == 'hvel' and not self.use_3d: # skip
                    pass
                elif var_name == 'dahv' and self.use_3d : # skip
                    pass   
                elif var_name == 'dahv' and not self.use_3d: # then use depth-averaged data
                    self.variable_mapping['x_sea_water_velocity'] = str(var_name)
                    self.variable_mapping['y_sea_water_velocity'] = str(var_name)                 
                else: # standard mapping                                    
                    self.variable_mapping[schism_mapping[var_name]] = \
                        str(var_name) 
                   
        self.variables = list(self.variable_mapping.keys())

        self.xmin = self.lon.min()
        self.xmax = self.lon.max()
        self.ymin = self.lat.min()
        self.ymax = self.lat.max()
        
        # Run constructor of parent Reader class
        super(Reader, self).__init__()
  
    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        """ The function extracts 'requested_variables' from the native SCHISM files
            which will then be used in _get_variables() inside get_variables_interpolated() 
            to initialise the ReaderBlockUnstruct objects used to interpolate data in space and time

            For now the function will extract the entire slice of data of 'requested_variables' at given 'time'

            We could consider extracting a subset of data around particles clouds to have less data in memory...
            but that means we'll need to recompute the KDtree of data every time in ReaderBlockUnstruct.
            ...Not sure which option would be fastest...to check

        """

        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

        # Initialising dictionary to contain data
        variables = {'x': self.lon, 'y': self.lat, 'z': 0.*self.lat,
                     'time': nearestTime}

        # extracts the full slices of requested_variables at time indxTime
        for par in requested_variables:
            if par not in ['x_sea_water_velocity','y_sea_water_velocity'] :
                # standard case - for all variables except current velocities
                var = self.Dataset.variables[self.variable_mapping[par]]
                if var.ndim == 1:
                    data = var[:]
                elif var.ndim == 2:
                    data = var[indxTime,:]
                elif var.ndim == 3:
                    data = var[indxTime,0,:]
                else:
                    raise ValueError('Wrong dimension of %s: %i' %
                                     (self.variable_mapping[par], var.ndim))
            else :
                # requested variables are current velocities
                # In SCHISM netcdf filesboth [u,v] components are saved 
                # as two different dimensions of the same variable.
                var = self.Dataset.variables[self.variable_mapping[par]]
                if var.ndim == 3: # 2D current data , dimensions [time,node,2]
                   if par == 'x_sea_water_velocity':
                       data = var[indxTime,:,0]
                   elif par == 'y_sea_water_velocity':
                       data = var[indxTime,:,1]
                elif var.ndim == 4: # 3D current data , dimensions [time,node,zcor,2]
                   logging.debug('3D interp not supported yet by reader_netCDF_CF_unstructured_SCHISM')
                   raise ValueError('Not supported - quitting')
                   import pdb;pdb.set_trace()
                   # Need to extract zcor here - find best way to do it...
                   # we will need to update the dictionary variables['z']
                   # 

                   # if par == 'x_sea_water_velocity':
                   #     data = var[indxTime,c,level,0]
                   # elif par == 'y_sea_water_velocity':
                   #     data = var[indxTime,c,level,1]
                else:
                    raise ValueError('Wrong dimension of %s: %i' %
                                     (self.variable_mapping[par], var.ndim))

            variables[par] = data # save all data slice to dictionary with key 'par'
        
        #print variables        
        return variables 


    def get_variables_depreciated(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        """This version mimics what would be output from a regular 
           grid reader (e.g. reader_netCDF_CF_generic) by returning 
           data from unstructured grid interpolated to a regular grid
           that includes all particles (+some buffer)

           The get_variables() method is called inside get_variables_interpolated()
           which is invoked by model class. The get_variables_interpolated()
           returns variable data interpolated in space anb time at particle positions

           > Depreciated now - see new get_variables() method that is fetching the native SCHISM data
           
        """

        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

        # Finding a subset around the particles, so that
        # we do not interpolate more points than is needed.
        # Performance is quite dependent on the given buffer,
        # but it should not be made too small to make sure
        # particles are inside box
        deg2meters = 1.1e5 # approx
        buffer = .1  # degrees around given positions
        buffer = .05  # degrees around given positions

        if x.max() <= 360 :
            buffer = buffer_deg # latlong reference system
        else:
            buffer = buffer * deg2meters# projected reference system in meters

        lonmin = x.min() - buffer
        lonmax = x.max() + buffer
        latmin = y.min() - buffer
        latmax = y.max() + buffer
        c = np.where((self.lon > lonmin) &
                     (self.lon < lonmax) &
                     (self.lat > latmin) &
                     (self.lat < latmax))[0]
        

        # Making a lon-lat grid onto which data is interpolated
        lonstep = .01 * 0.25 # hardcoded for now
        latstep = .01 * 0.25 # hardcoded for now
        if x.max() >= 360 :
            lonstep = lonstep * deg2meters
            latstep = latstep * deg2meters
        
        lons = np.arange(lonmin, lonmax, lonstep)
        lats = np.arange(latmin, latmax, latstep)
        lonsm, latsm = np.meshgrid(lons, lats)

        # Initialising dictionary to contain data
        variables = {'x': lons, 'y': lats, 'z': z,
                     'time': nearestTime}

        # Reader coordinates of subset
        for par in requested_variables:
            if par not in ['x_sea_water_velocity','y_sea_water_velocity'] :
                # standard case - for all variables except current velocities
                var = self.Dataset.variables[self.variable_mapping[par]]
                if var.ndim == 1:
                    data = var[c]
                elif var.ndim == 2:
                    data = var[indxTime,c]
                elif var.ndim == 3:
                    data = var[indxTime,0,c]
                else:
                    raise ValueError('Wrong dimension of %s: %i' %
                                     (self.variable_mapping[par], var.ndim))
            else :
                # requested variables are current velocities
                # In SCHISM netcdf filesboth [u,v] components are saved 
                # as two different dimensions of the same variable.
                var = self.Dataset.variables[self.variable_mapping[par]]
                if var.ndim == 3: # 2D current data , dimensions [time,node,2]
                   if par == 'x_sea_water_velocity':
                       data = var[indxTime,c,0]
                   elif par == 'y_sea_water_velocity':
                       data = var[indxTime,c,1]
                elif var.ndim == 4: # 3D current data , dimensions [time,node,zcor,2]
                   logging.debug('3D interp not supported yet by reader_netCDF_CF_unstructured_SCHISM')
                   raise ValueError('Not supported - quitting')
                   import pdb;pdb.set_trace()
                   # if par == 'x_sea_water_velocity':
                   #     data = var[indxTime,c,level,0]
                   # elif par == 'y_sea_water_velocity':
                   #     data = var[indxTime,c,level,1]
                else:
                    raise ValueError('Wrong dimension of %s: %i' %
                                     (self.variable_mapping[par], var.ndim))

            # now interpolate extracted data points onto regular grid
            if 'interpolator' not in locals():
                logging.debug('Making interpolator...')
                interpolator = LinearNDInterpolator((self.lat[c],
                                                     self.lon[c]),
                                                    data)
            else:
                # Re-use interpolator for other variables
                interpolator.values[:,0] = data
            interpolator((0,0))

            variables[par] = interpolator(latsm, lonsm)
        
        #print variables        
        return variables

    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader."""
        
        # weird bug in which x,y become 1e30 sometimes...
        # need to print stuff to make it work...

        # Calculate x,y coordinates from lon,lat
        print(lon.max()) # if I comment this...[x,y] may become 1e+30..have no idea why
                         # it runs fine if I have that print statement ....
        # print(lat.max())
        x, y = self.lonlat2xy(lon, lat)      
        # Only checking vertical coverage if zmin, zmax is defined
        zmin = -np.inf
        zmax = np.inf
        if hasattr(self, 'zmin') and self.zmin is not None:
            zmin = self.zmin
        if hasattr(self, 'zmax') and self.zmax is not None:
            zmax = self.zmax

        # sometimes x,y will are = 1e30 at this stage..and have no idea why
        # recomputing them below seems to fix the issue        
        # x, y = self.lonlat2xy(lon, lat)

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        else:
            indices = np.where((x >= self.xmin) & (x <= self.xmax) &
                               (y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        
        try:
            return indices, x[indices], y[indices]
        except:
            return indices, x, y


    def get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):
        """This function will overloads the get_variables_interpolated() methods
           available in basereader.py. 

           The get_variables_interpolated() from basereader.py used regularly gridded 
           data "blocks" extracted from the netcdf files ( which may possibly be "cached" 
           for speed improvements - see code for more detail). A workaround for unstructured 
           grids was to create a regularly gridded block from unstructured data 
           (see get_variables() methods above), and then pass it to the 
           basereader's get_variables_interpolated(). 

           This function instead uses native SCHISM netcdf files to interpolate data 
           in space and time and returns 'variables' at particle positions [lon,lat] 
           This will make better use of the high-resolution provided by unstructured grids.

        """

# >> objective is to have a method here that overload get_variables_interpolated()
# using directly the data from SCHISM, and scipy etc..or any other 3d irreg interp routine

        self.timer_start('total')
        self.timer_start('preparing')
        # Raise error if time not not within coverage of reader
        if not self.covers_time(time):
            raise ValueError('%s is outside time coverage (%s - %s) of %s' %
                             (time, self.start_time, self.end_time, self.name))

        # Check which particles are covered (indep of time)
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
        # For variables which are not time dependent, we do not care about time
        static_variables = ['sea_floor_depth_below_sea_level', 'land_binary_mask']
        if time == time_before or all(v in static_variables for v in variables):
            time_after = None

        z = z.copy()[ind_covered]  # Send values and not reference
        

        # start interpolation procedure

        # 
        # general idea is to create a  new "ReaderBlockUnstruct" class that will be called instead of
        # the regular "ReaderBlock" as in basereader.py


        # for now we are not doing any kind of caching or re-using, only reading "slices" of netcdf files
        # and interpolating from that 
        # This could be potentially improved using a similar approach as what is done
        # in get_variables_interpolated() from basereader.py
     
        block_before = block_after = None
        blockvariables_before = variables
        blockvars_before = str(variables)
        blockvariables_after = variables
        blockvars_after = str(variables)
        for blockvars in self.var_block_before:
            if all(v in blockvars for v in variables):
                block_before = self.var_block_before[blockvars]
                blockvariables_before = block_before.data_dict.keys()
                blockvars_before = blockvars
                break
            blockvariables_before = variables
            blockvars_before = str(variables)
        for blockvars in self.var_block_after:
            if all(v in blockvars for v in variables):
                block_after = self.var_block_after[blockvars]
                blockvariables_after = block_after.data_dict.keys()
                blockvars_after = blockvars
                break

        # Swap before- and after-blocks if matching times
        if block_before is not None and block_after is not None:
            if block_before.time != time_before:
                if block_after.time == time_before:
                    block_before = block_after
                    self.var_block_before[blockvars_before] = block_before
            if block_after.time != time_after:
                if block_before.time == time_before:
                    block_after = block_before
                    self.var_block_after[blockvars_after] = block_after
        # Fetch data, if no buffer is available
        if block_before is None or \
                block_before.time != time_before:
            self.timer_end('preparing')
            reader_data_dict = \
                self._get_variables(blockvariables_before, profiles,
                                    profiles_depth, time_before,
                                    reader_x, reader_y, z,
                                    block=block)
            # now use reader_data_dict to initialize 
            # a ReaderBlockUnstruct
            self.timer_start('preparing')
            self.var_block_before[blockvars_before] = \
                ReaderBlockUnstruct(reader_data_dict,
                    KDtree = self.KDtree,
                    interpolation_horizontal=self.interpolation)
            try:
                len_z = len(self.var_block_before[blockvars_before].z)
            except:
                len_z = 1
            logging.debug(('Fetched env-block (size %ix%ix%i) ' +
                          'for time before (%s)') %
                          (len(self.var_block_before[blockvars_before].x),
                           len(self.var_block_before[blockvars_before].y),
                           len_z, time_before))
            block_before = self.var_block_before[blockvars_before]
        if block_after is None or block_after.time != time_after:
            if time_after is None:
                self.var_block_after[blockvars_after] = \
                    block_before
            else:
                self.timer_end('preparing')
                reader_data_dict = \
                    self._get_variables(blockvariables_after, profiles,
                                        profiles_depth, time_after,
                                        reader_x, reader_y, z,
                                        block=block)
                self.timer_start('preparing')
                self.var_block_after[blockvars_after] = \
                    ReaderBlockUnstruct(
                        reader_data_dict,
                        interpolation_horizontal=self.interpolation)
                try:
                    len_z = len(self.var_block_after[blockvars_after].z)
                except:
                    len_z = 1

                logging.debug(('Fetched env-block (size %ix%ix%i) ' +
                              'for time after (%s)') %
                              (len(self.var_block_after[blockvars_after].x),
                               len(self.var_block_after[blockvars_after].y),
                               len_z, time_after))
                block_after = self.var_block_after[blockvars_after]
        
        if (block_before is not None and block_before.covers_positions(
            reader_x, reader_y) is False) or (\
            block_after is not None and block_after.covers_positions(
                reader_x, reader_y) is False):
            import pdb;pdb.set_trace()
            logging.warning('Data block from %s not large enough to '
                            'cover element positions within timestep. '
                            'Buffer size (%s) must be increased.' %
                            (self.name, str(self.buffer)))

        self.timer_end('preparing')

        # here we need to create a block_before block_after interpolator 

        ############################################################
        # Interpolate before/after blocks onto particles in space
        ############################################################
        self.timer_start('interpolation')
        logging.debug('Interpolating before (%s) in space  (%s)' %
                      (block_before.time, self.interpolation))
        env_before, env_profiles_before = block_before.interpolate(
                reader_x, reader_y, z, variables,
                profiles, profiles_depth)

        if (time_after is not None) and (time_before != time):
            logging.debug('Interpolating after (%s) in space  (%s)' %
                          (block_after.time, self.interpolation))
            env_after, env_profiles_after = block_after.interpolate(
                    reader_x, reader_y, z, variables,
                    profiles, profiles_depth)

        self.timer_end('interpolation')

        #######################
        # Time interpolation
        #######################
        self.timer_start('interpolation_time')
        env_profiles = None
        if (time_after is not None) and (time_before != time) and self.return_block is True:
            weight_after = ((time - time_before).total_seconds() /
                            (time_after - time_before).total_seconds())
            logging.debug(('Interpolating before (%s, weight %.2f) and'
                           '\n\t\t      after (%s, weight %.2f) in time') %
                          (block_before.time, 1 - weight_after,
                           block_after.time, weight_after))
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
                    env_profiles_before[var]=np.atleast_2d(env_profiles_before[var])
                    env_profiles_after[var]=np.atleast_2d(env_profiles_after[var])
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
                            env_profiles[vector_pair[0]], env_profiles[vector_pair[1]] = \
                                    self.rotate_vectors (reader_x, reader_y,
                                            env_profiles[vector_pair[0]],
                                            env_profiles[vector_pair[1]],
                                            self.proj,
                                            rotate_to_proj)


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
        self.timer_end('total')
        return env, env_profiles



###########################
# ReaderBlockUnstruct
###########################

# horizontal_interpolation_methods = {
#     'nearest': Nearest2DInterpolator,
#     'ndimage': NDImage2DInterpolator,
#     'linearND': LinearND2DInterpolator,
#     'linearNDFast': Linear2DInterpolator}


# vertical_interpolation_methods = {
#     'nearest': Nearest1DInterpolator,
#     'linear': Linear1DInterpolator}


class ReaderBlockUnstruct():
    """Class to store and interpolate the data from an *unstructured* reader.
       This is the equivalent of ReaderBlock (regular grid) for *unstructured* grids.

       arguments: (in addition to ReaderBlock ones)

           KDtree : for nearest-neighbor search (initialized using SCHISM nodes in reader's _init_() )
                    This is read from reader object so that it is not recomputed every time

    """

    def __init__(self, data_dict, 
                 KDtree = None,
                 interpolation_horizontal='linearNDFast',
                 interpolation_vertical='linear'):

        # Make pointers to data values, for convenience
        self.x = data_dict['x']
        self.y = data_dict['y']
        self.time = data_dict['time']
        self.data_dict = data_dict
        del self.data_dict['x']
        del self.data_dict['y']
        del self.data_dict['time']
        try:
            self.z = data_dict['z']
            del self.data_dict['z']
        except:
            self.z = None

        # save KDtree
        self.KDtree = KDtree

        # Mask any extremely large values, e.g. if missing netCDF _Fill_value
        filled_variables = set()
        for var in self.data_dict:
            if isinstance(self.data_dict[var], np.ma.core.MaskedArray):
                self.data_dict[var] = np.ma.masked_outside(
                    np.ma.masked_invalid(self.data_dict[var]), -1E+9, 1E+9)
                # Convert masked arrays to numpy arrays
                self.data_dict[var] = np.ma.filled(self.data_dict[var],
                                                   fill_value=np.nan)
            # Fill missing data towards seafloor if 3D
            if isinstance(self.data_dict[var], (list,)):
                logging.warning('Ensemble data currently not extrapolated towards seafloor')
            elif self.data_dict[var].ndim == 3:
                filled = fill_NaN_towards_seafloor(self.data_dict[var])
                if filled is True:
                    filled_variables.add(var)
                
        if len(filled_variables) > 0:
            logging.debug('Filled NaN-values toward seafloor for :'
                          + str(list(filled_variables)))
        
        # below may not be relevant any longer
        if False:
            # 
            # 
            # Set 1D (vertical) and 2D (horizontal) interpolators
            try:
                self.Interpolator2DClass = \
                    horizontal_interpolation_methods[interpolation_horizontal]
            except Exception:
                raise NotImplementedError(
                    'Valid interpolation methods are: ' +
                    str(horizontal_interpolation_methods.keys()))

            try:
                self.Interpolator1DClass = \
                    vertical_interpolation_methods[interpolation_vertical]
            except Exception:
                raise NotImplementedError(
                    'Valid interpolation methods are: ' +
                    str(vertical_interpolation_methods.keys()))

        if 'land_binary_mask' in self.data_dict.keys() and \
                interpolation_horizontal != 'nearest':
            logging.debug('Nearest interpolation will be used '
                          'for landmask, and %s for other variables'
                          % interpolation_horizontal)

    def _initialize_interpolator(self, x, y, z=None):
        logging.debug('Initialising interpolator.')
        self.interpolator2d = self.Interpolator2DClass(self.x, self.y, x, y)
        if self.z is not None and len(np.atleast_1d(self.z)) > 1:
            self.interpolator1d = self.Interpolator1DClass(self.z, z)

    def interpolate(self, x, y, z=None, variables=None,
                    profiles=[], profiles_depth=None):
        # Use the KDtree to interpolate data to [x,y,z] particle positions

        # self._initialize_interpolator(x, y, z)

        env_dict = {}
        if profiles is not []:
            profiles_dict = {'z': self.z}
        for varname, data in iteritems(self.data_dict):
            nearest = False
            # land mask 
            if varname == 'land_binary_mask':
                nearest = True
                self.interpolator2d_nearest = Nearest2DInterpolator(self.x, self.y, x, y)
            # ensemble data
            if type(data) is list:
                num_ensembles = len(data)
                logging.debug('Interpolating %i ensembles for %s' % (num_ensembles, varname))
                if data[0].ndim == 2:
                    horizontal = np.zeros(x.shape)*np.nan
                else:
                    horizontal = np.zeros((len(self.z), len(x)))*np.nan
                ensemble_number = np.remainder(range(len(x)), num_ensembles)
                for en in range(num_ensembles):
                    elnum = ensemble_number == en
                    int_full = self._interpolate_horizontal_layers(data[en], nearest=nearest)
                    if int_full.ndim == 1:
                        horizontal[elnum] = int_full[elnum]
                    else:
                        horizontal[:, elnum] = int_full[:, elnum]
            # horizontal interpolation
            else:
                
                # use KDtree to find nearest neighbours and interpolate based on distance
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html#scipy.spatial.cKDTree.query
                nb_closest_nodes = 3
                DMIN=1.e-10
                try:
                    dist,i=self.KDtree.query(np.vstack((x,y)).T,nb_closest_nodes, n_jobs=-1) #quick nearest-neighbor lookup
                except:
                    import pdb;pdb.set_trace()

                    ..looks like reference to KDtree is lost somehow...
                    self.KDtree becomes None ???

                #dist = distance to nodes
                #i = index of nodes 
                
                # 
                # if i.max()>=np.vstack((x,y)).T.shape[0]: # if max node_id is larger than number of nodes
                #   raise DataException('Finite element interpolation out of range')

                dist[dist<DMIN]=DMIN
                fac=(1./dist)
                
                # spatial interpolation by applying
                # distance-based 'fac' to data at nearest nodes 
                if data.ndim ==1 : # no vertical levels
                    horizontal = (fac*data.take(i)).sum(-1)/fac.sum(-1)
                elif dat.ndim==2: # data has vertical levels
                    # TO DO ! need to use zcor
                    tmp=(fac*data.take(i,1)).sum(-1)/fac.sum(-1) 
                    horizontal = interpz(tmp.astype('f'),p[:,2],self.lev)
                else:
                    import pdb;pdb.set_trace()
                    horizontal = (fac*data.take(i)).sum(-1)/fac.sum(-1)
                
                # TO DO : mask out of grid points - setting velocities to 0.0
                # datout[numpy.where(~self.ingrid(p) )] = 0.0 
 
                # horizontal = self._interpolate_horizontal_layers(data, nearest=nearest)

            if profiles is not None and varname in profiles:
                profiles_dict[varname] = horizontal
            if horizontal.ndim > 1:
                env_dict[varname] = self.interpolator1d(horizontal)
            else:
                env_dict[varname] = horizontal
        if 'z' in profiles_dict:
            profiles_dict['z'] = np.atleast_1d(profiles_dict['z'])
        
        self.KDtree
        import pdb;pdb.set_trace()

        return env_dict, profiles_dict

    def _interpolate_horizontal_layers(self, data, nearest=False):
        '''Interpolate all layers of 3d (or 2d) array.'''
    
        if nearest is True:
            interpolator2d = self.interpolator2d_nearest
        else:
            interpolator2d = self.interpolator2d
        if data.ndim == 2:
            return interpolator2d(data)
        if data.ndim == 3:
            num_layers = data.shape[0]
            # Allocate output array
            result = np.ma.empty((num_layers, len(interpolator2d.x)))
            for layer in range(num_layers):
                result[layer, :] = self.interpolator2d(data[layer, :, :])
            return result

    def covers_positions(self, x, y, z=None):
        '''Check if given positions are covered by this reader block.'''
        
        indices = np.where((x >= self.x.min()) & (x <= self.x.max()) &
                           (y >= self.y.min()) & (y <= self.y.max()))[0]
        if len(indices) == len(x):
            return True
        else:
            return False