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
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator

from opendrift.readers.basereader import BaseReader, pyproj


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
        """For now this mimics what would be output from a regular 
           grid reader (e.g. reader_netCDF_CF_generic) by returning 
           data from unstructured grid interpolated to a regular grid
           that includes all particles (+some buffer)

           The get_variables() method is called inside get_variables_interpolated()
           which is invoked by model class. The get_variables_interpolated()
           returns variable data interpolated in space anb time at particle positions
           
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
        # import pdb;pdb.set_trace()
        # Calculate x,y coordinates from lon,lat
        x, y = self.lonlat2xy(lon, lat)
        print(x.max())
        # print(y.max())
        # Only checking vertical coverage if zmin, zmax is defined
        zmin = -np.inf
        zmax = np.inf
        if hasattr(self, 'zmin') and self.zmin is not None:
            zmin = self.zmin
        if hasattr(self, 'zmax') and self.zmax is not None:
            zmax = self.zmax

        # sometimes x,y will are = 1e30 at this stage..and have no idea why
        # recomputing them below seems to fix the issue        
        x, y = self.lonlat2xy(lon, lat)

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        else:
            indices = np.where((x >= self.xmin) & (x <= self.xmax) &
                               (y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        
        # print(x.max())
        # print(y.max())
        # print(indices)
        if len(indices)==0: import pdb;pdb.set_trace()
        try:
            return indices, x[indices], y[indices]
        except:
            return indices, x, y