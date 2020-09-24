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
#
#####################################
# This reader support native selfie output files and will 
# handle interpolation of current velocities and other variables
# in 2D and 3D.
# The interpolation is done using a KDtree-approach defined
# form the 2D or 3D mesh nodes
# 
# NOTE:
# This reader is under development,
# and presently not fully functional
# In particular, the selfe model time is referenced to seconds since the start of the model run,
# line 267 this needs to be updated between applications to update for the new timestamp
# 
# Author: Simon Weppe. MetOcean Solution, MetService New Zealand
# Author: Ammended by Alice Goward Brown, MetOcean Solutions, MetService New Zealand
#########################################################################

import logging
import numpy as np
from datetime import datetime
from future.utils import iteritems
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree #cython-based KDtree for quick nearest-neighbor search
from scipy.spatial import ConvexHull # convex hull of grid points
from matplotlib.path import Path # convex hull of grid points
from opendrift.readers.basereader import BaseReader, pyproj

try:
    import xarray as xr
    has_xarray = True
except:
    has_xarray = False


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
        """Initialise reader_netCDF_CF_unstructured_selfie

        Args:
            filename:  name of selfie netcdf file (can have wildcards)
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
        
        # [name_used_in_selfie : equivalent_CF_name]
        selfie_mapping = {
            'U_dahv': 'x_sea_water_velocity',
            'V_dahv': 'y_sea_water_velocity',
            'U_hvel': 'x_sea_water_velocity',
            'V_hvel': 'y_sea_water_velocity',
            'h': 'sea_floor_depth_below_sea_level',
            'elev' : 'sea_surface_height',
            'temp' : 'sea_water_temperature',
            'salt' : 'sea_water_salinity',
            'zcor' : 'vertical_levels', # time-varying vertical coordinates
            'siglay': 'ocean_s_coordinate',
            'U_wind':'x_wind',
            'V_wind':'y_wind'}

        self.return_block = True

        try:
            # Open file, check that everything is ok
            logging.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                logging.info('Opening files with MFDataset')
                if has_xarray:
                    self.Dataset = xr.open_mfdataset(filename)
                else:
                    self.Dataset = MFDataset(filename,aggdim='time')

            else:
                logging.info('Opening file with Dataset')
                if has_xarray:
                    self.Dataset = xr.open_dataset(filename)
                else:
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
            if 'U_hvel' or 'V_hvel' in self.Dataset.variables:
                self.use_3d = True
            else:
                self.use_3d = False

        if self.use_3d and 'U_hvel' and 'V_hvel' not in self.Dataset.variables:
            logging.debug('No 3D velocity data in file - cannot find variable ''hvel'' ')
        elif self.use_3d and 'U_hvel' and 'V_hvel' in self.Dataset.variables:
            if 'zcor' in self.Dataset.variables: # both hvel and zcor in files - all good
                self.nb_levels = self.Dataset.variables['U_hvel'].shape[1] #hvel dimensions : [time,lev,node]
            else:
                logging.debug('No vertical level information present in file ''zcor'' ... stopping')
                raise ValueError('variable ''zcor'' must be present in netcdf file to be able to use 3D currents')
                        

        logging.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        for var_name in self.Dataset.variables:
            
            if var_name in ['Easting','Northing']:
                # all these variables have the same standard name projection_x_coordinate,projection_y_coordinate
                # we need to only use :
                # SCHISM_hgrid_node_x as projection_x_coordinate
                # SCHISM_hgrid_node_y as projection_y_coordinate
                # which are the nodes (also called vertices)
                continue

            var = self.Dataset.variables[var_name]
            if var.ndim > 1:
                continue  # Coordinates must be 1D-array
            if has_xarray:
                attributes = var.attrs
                att_dict = var.attrs
            else:
                attributes = var.ncattrs()
                att_dict = var.__dict__
            # attributes = var.ncattrs()
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            CoordinateAxisType = ''
            # add checks on projection here ? 
            # as in reader_netCDF_CF_generic.py
            if 'standard_name' in attributes:
                standard_name = att_dict['standard_name']
            if 'long_name' in attributes:
                long_name = att_dict['long_name']
            if 'axis' in attributes:
                axis = att_dict['axis']
            if 'units' in attributes:
                units = att_dict['units']
            if '_CoordinateAxisType' in attributes:
                CoordinateAxisType = att_dict['_CoordinateAxisType']

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
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                x = var_data*unitfactor
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
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                y = var_data*unitfactor
                self.numy = var.shape[0]
            if standard_name == 'depth' or axis == 'Z':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]

                if 'positive' not in var.ncattrs() or \
                        var.__dict__['positive'] == 'up': 
                    self.z = var_data
                else:
                    self.z = -var_data

            if standard_name == 'time' or axis == 'T' or var_name == 'time':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                # Read and store time coverage (of this particular file)
                time = var_data
                time_units = units
                # self.times = num2date(time, time_units)
                if has_xarray:
                    self.times = [datetime.utcfromtimestamp((OT -
                        np.datetime64('1970-01-01T00:00:00Z')
                            ) / np.timedelta64(1, 's')) for OT in time]
                else:

                    self.times = num2date(time, "seconds since 2002-01-01T00:00:00") # temporary fix to check

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
        self.reader_KDtree = cKDTree(np.vstack((self.lon,self.lat)).T) 
        #do we need copy_data=True ..probably not since self.lon,self.lat wont change.

        # build convex hull of points for particle-in-mesh checks
        logging.debug('Building convex hull of nodes for particle''s in-mesh checks')
        # *** this a very approximate way of doing it...ideally we would be able to 
        # have the actual contour of the grid....
        hull_obj = ConvexHull(np.vstack([self.lon.data,self.lat.data]).T)
        self.hull_path = Path(np.vstack([self.lon[hull_obj.vertices],self.lat[hull_obj.vertices]]).T)  # Matplotlib Path object
        self.hull = np.vstack([self.lon[hull_obj.vertices],self.lat[hull_obj.vertices]]).T


        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname]: #'depth'
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]

            if has_xarray:
                attributes = var.attrs
                att_dict = var.attrs
            else:
                attributes = var.ncattrs()
                att_dict = var.__dict__

            if var_name in selfie_mapping:
                self.variable_mapping[selfie_mapping[var_name]] = str(var_name) 
                    
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
            # standard case - for all variables except current velocities
            var = self.Dataset.variables[self.variable_mapping[par]]
            if var.ndim == 1:
                data = var[:] # e.g. depth
                logging.debug('reading constant data from unstructured reader %s' % (par))
            elif var.ndim == 2: 
                data = var[indxTime,:] # e.g. dahv
                logging.debug('reading 2D data from unstructured reader %s' % (par))
            elif var.ndim == 3:
                data = var[indxTime,:,:] # e.g. 3d salt [time,lev,node] different to SCHISM which is [time,node,lev]
                logging.debug('reading 3D data from unstructured reader %s' % (par))
                # convert 3D data matrix to one column array and define corresponding data coordinates [x,y,z]
                # (+ update variables dictionary with 3d coords if needed)
                data,variables = self.convert_3d_to_array(indxTime,data,variables)
            else:
                raise ValueError('Wrong dimension of %s: %i' %
                                (self.variable_mapping[par], var.ndim))
                    
            variables[par] = data # save all data slice to dictionary with key 'par'

            if has_xarray is True:
                variables[par] = np.asarray(variables[par])
            
            # Store coordinates of returned points
            #  >> done in previous steps here, line 380 and in convert_3d_to_array()

        return variables 


    def convert_3d_to_array(self,id_time,data,variable_dict):
        ''' 
        The function reshapes a data matrix of dimensions = [node,vertical_levels] (i.e. data at vertical levels, at given time step) 
        into a one-column array and works out corresponding 3d coordinates [lon,lat,z] using the time-varying
        'zcor' variable in SCHISM files (i.e. vertical level positions). 

        These 3D coordinates will be used to build the 3D KDtree for data interpolation and will be added to the 'variable_dict' 
        which is eventually passed to get_variables_interpolated().

        args:
            -id_time
            -data
            -variable_dict
        out :
            -flattened 'data' array
            -addition of ['x_3d','y_3d','z_3d'] items to variable_dict if needed.

        '''
        try:
            vertical_levels = self.Dataset.variables['zcor'][id_time,:,:]
            # depth are negative down consistent with convention used in OpenDrift 
            # if using the netCDF4 library, vertical_levels is masked array where "masked" levels are those below seabed  (= 9.9692100e+36)
            # if using the xarray library, vertical_levels is nan for levels are those below seabed
            if has_xarray:
                # convert to masked array to be consistent with what netCDF4 lib returns
                vertical_levels = np.ma.array(vertical_levels, mask = np.isnan(vertical_levels.data)) 
                data = np.asarray(data)
                # vertical_levels.mask = np.isnan(vertical_levels.data) # masked using nan's when using xarray
        except:
            logging.debug('no vertical level information present in file ''zcor'' ... stopping')
            raise ValueError('variable ''zcor'' must be present in netcdf file to be able to use 3D currents')
        # flatten 3D data 
        data = np.ravel(data[~vertical_levels.mask])

        # add corresponding 3D coordinates to the 'variable_dict' which is eventually passed to get_variables_interpolated()
        # Note: these 3d coordinates will be the same for all 3d data (current,salt/temp etc..), and will change at each reader time step
        # They are saved as ['x_3d','y_3d','z_3d'] rather than ['x','y','z'] as it would break things when both 2d and 3d data are requested.
        if 'z_3d' not in variable_dict.keys():
            # now make a long 3-column array [lon,lat,z] [n_nodes x  3] at which 'data' is defined
            # tile lon/lat data so that array shape match vertical_levels shape
            lon_tiled = np.tile(self.lon,(self.nb_levels,1)) # dimensions [node,vertical_levels]
            lat_tiled = np.tile(self.lat,(self.nb_levels,1))
            # arrays are tiled so that lon_tiled[0,:] = return same value i.e. same lon/lat for all z levels 
            if lon_tiled.shape != vertical_levels.shape:
                import pdb;pdb.set_trace()
            # convert to masked array consistent with vertical_levels
            lon_tiled_ma = np.ma.array(lon_tiled, mask = vertical_levels.mask) 
            lat_tiled_ma = np.ma.array(lat_tiled, mask = vertical_levels.mask)
            # flatten arrays and add to dictionary
            variable_dict['x_3d'] = np.ravel(lon_tiled_ma[~vertical_levels.mask]) 
            variable_dict['y_3d'] = np.ravel(lat_tiled_ma[~vertical_levels.mask])
            variable_dict['z_3d'] = np.ravel(vertical_levels[~vertical_levels.mask])
        
        return data,variable_dict
        
    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader.
        
        For an unstructured reader, this is done by checking that if particle positions are within the convex hull
        of the mesh nodes. This means it is NOT using the true polygon bounding the mesh ..but this is generally good enough
        to locate the outer boudary (which is often circular). 
        >> The convex hull will not be good for the shorelines though, but this 
        is not critical since coast interaction will be handled by either by the landmask (read from file or using global_landmask) 
        Data tnterpolation might be an issue for particles reaching the land, and not actually flagged as out-of-bounds...
        To Check...

        Better alternatives could be: 
            - get the true mesh polygon from netCDF files..doesnt seem to be available.
            - compute an alpha-shape of mesh nodes (i.e. enveloppe) and use that as polygon. better than convexhull but not perfect 
                    (shp = alphaShape(double(x),double(y),10000); work well in matlab for example)
            - workout that outer mesh polygon from a cloud of points...maybe using a more evolved trimesh package ?

        """        
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
            pass
            # unlikely to be used for SCHISM domain
        else:
            # indices = np.where((x >= self.xmin) & (x <= self.xmax) &
            #        (y >= self.ymin) & (y <= self.ymax) &
            #        (z >= zmin) & (z <= zmax))[0]

            # use the Path object (matplotlib) defined in __init__() : self.hull_path
            in_hull =self.hull_path.contains_points(np.vstack([x,y]).T)
            indices = np.where(in_hull) 

        try:
            return indices, x[indices], y[indices]
        except:
            return indices, x, y

    def get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):
        """This function will overload the get_variables_interpolated() methods
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
        # This allows re-using almost the same code a basereader.py for the block_before/block_after
     
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
            logging.debug('initialize ReaderBlockUnstruct var_block_before')
            self.var_block_before[blockvars_before] = \
                ReaderBlockUnstruct(reader_data_dict,
                    KDtree = self.reader_KDtree,
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
                logging.debug('initialize ReaderBlockUnstruct var_block_after')
                self.var_block_after[blockvars_after] = \
                    ReaderBlockUnstruct(
                        reader_data_dict,
                        KDtree = self.reader_KDtree,
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

    def plot(self, variable=None, vmin=None, vmax=None,
             filename=None, title=None, buffer=1, lscale='auto'):
        """Plot geographical coverage of reader."""

        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from opendrift_landmask_data import Landmask
        fig = plt.figure()

        #####################################
        # In Dev - plot unstructured mesh
        # 
        # To do 
        #  - plot in correct cartopy projection, consistent with basereader.py plot()
        #  - plot a given timestep for flows
        #  - incorporate that to animation as in basemodel.py etc...
        # 
        #####################################
        import matplotlib.tri as mtri
        from mpl_toolkits.mplot3d import Axes3D
        X=self.Dataset.variables['Easting'][:]
        Y=self.Dataset.variables['Northing'][:]
        Z=self.Dataset.variables['depth'][:]
        face=self.Dataset.variables['nv'][:,0:3]-1
        # build triangulation
        triang =mtri.Triangulation(X,Y, triangles=face, mask=None)
        # quads=Triangulation(gr.x,gr.y,quads)
        #
        # plot the triangular mesh 
        plt.triplot(triang) 
        # plot the variable - here depth for now
        plt.tricontourf(triang, Z)
        # more complex 3D plots 
        # ax = fig.add_subplot(1,1,1, projection='3d')
        # ax.plot_trisurf(triang,Z, cmap='jet')
        # or 
        # from https://github.com/metocean/schism-public/blob/master/LSC/lsc.py
        # tricon=ax.tricontourf(triang,Z,cmap='jet')
        # ax.view_init(azim=-90, elev=90)
        plt.ion()
        plt.show()
        import pdb;pdb.set_trace()
        #####################################


        corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                 [self.ymax, self.ymin, self.ymax, self.ymin])
        lonmin = np.min(corners[0]) - buffer*2
        lonmax = np.max(corners[0]) + buffer*2
        latmin = np.min(corners[1]) - buffer
        latmax = np.max(corners[1]) + buffer
        latspan = latmax - latmin

        # Initialise map
        if latspan < 90:
            # Stereographic projection centred on domain, if small domain
            x0 = (self.xmin + self.xmax) / 2
            y0 = (self.ymin + self.ymax) / 2
            lon0, lat0 = self.xy2lonlat(x0, y0)
            sp = ccrs.Stereographic(central_longitude=lon0, central_latitude=lat0)
            ax = fig.add_subplot(1, 1, 1, projection=sp)
            corners_stere = sp.transform_points(ccrs.PlateCarree(), np.array(corners[0]), np.array(corners[1]))
        else:
            # Global map if reader domain is large
            sp = ccrs.Mercator()
            ax = fig.add_subplot(1, 1, 1, projection=sp)
            #map = Basemap(np.array(corners[0]).min(), -89,
            #              np.array(corners[0]).max(), 89,
            #              resolution='c', projection='cyl')

        # GSHHS coastlines
        f = cfeature.GSHHSFeature(scale=lscale, levels=[1],
                                  facecolor=cfeature.COLORS['land'])
        ax.add_geometries(
            f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
            ccrs.PlateCarree(),
            facecolor=cfeature.COLORS['land'],
            edgecolor='black')

        gl = ax.gridlines(ccrs.PlateCarree())
        gl.xlabels_top = False

        # Get boundary
        npoints = 10  # points per side
        x = np.array([])
        y = np.array([])
        x = np.concatenate((x, np.linspace(self.xmin, self.xmax, npoints)))
        y = np.concatenate((y, [self.ymin]*npoints))
        x = np.concatenate((x, [self.xmax]*npoints))
        y = np.concatenate((y, np.linspace(self.ymin, self.ymax, npoints)))
        x = np.concatenate((x, np.linspace(self.xmax, self.xmin, npoints)))
        y = np.concatenate((y, [self.ymax]*npoints))
        x = np.concatenate((x, [self.xmin]*npoints))
        y = np.concatenate((y, np.linspace(self.ymax, self.ymin, npoints)))
        # from x/y vectors create a Patch to be added to map
        lon, lat = self.xy2lonlat(x, y)
        lat[lat>89] = 89.
        lat[lat<-89] = -89.
        p = sp.transform_points(ccrs.PlateCarree(), lon, lat)
        xsp = p[:, 0]
        ysp = p[:, 1]

        if variable is None:

            boundary = Polygon(list(zip(xsp, ysp)), alpha=0.5, ec='k', fc='b',
                               zorder=100)
            ax.add_patch(boundary)
            buf = (xsp.max()-xsp.min())*.1  # Some whitespace around polygon
            buf = 0
            try:
                ax.set_extent([xsp.min()-buf, xsp.max()+buf, ysp.min()-buf, ysp.max()+buf], crs=sp)
            except:
                pass
        if title is None:
            plt.title(self.name)
        else:
            plt.title(title)
        plt.xlabel('Time coverage: %s to %s' %
                   (self.start_time, self.end_time))

        if variable is not None:
            rx = np.array([self.xmin, self.xmax])
            ry = np.array([self.ymin, self.ymax])
            data = self.get_variables_derived(variable, self.start_time,
                                      rx, ry, block=True)
            rx, ry = np.meshgrid(data['x'], data['y'])
            rlon, rlat = self.xy2lonlat(rx, ry)
            #map_x, map_y = map(rlon, rlat, inverse=False)
            data[variable] = np.ma.masked_invalid(data[variable])
            if hasattr(self, 'convolve'):
                from scipy import ndimage
                N = self.convolve
                if isinstance(N, (int, np.integer)):
                    kernel = np.ones((N, N))
                    kernel = kernel/kernel.sum()
                else:
                    kernel = N
                self.logger.debug('Convolving variables with kernel: %s' % kernel)
                data[variable] = ndimage.convolve(
                            data[variable], kernel, mode='nearest')
            #map.pcolormesh(map_x, map_y, data[variable], vmin=vmin, vmax=vmax)
            ax.pcolormesh(rlon, rlat, data[variable], vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
            #cbar = map.colorbar()
            #cbar.set_label(variable)

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()




###########################
# ReaderBlockUnstruct class
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

       arguments: (in addition to ReaderBlock)

           KDtree : for nearest-neighbor search (initialized using selfie's nodes in reader's _init_() )
                    This is read from reader object, so that it is not recomputed every time

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
            # probably not valid ..since z are different for each point in selfie..rather than fixed
            # This is used for the profile interpolation that is not yet functional 
            #
            del self.data_dict['z']
        except:
            self.z = None    
        # if some 3d data is provided, save additional dict entries
        if 'z_3d' in data_dict.keys():
            self.x_3d = data_dict['x_3d']
            self.y_3d = data_dict['y_3d']
            self.z_3d = data_dict['z_3d'] 
            del self.data_dict['x_3d']
            del self.data_dict['y_3d']
            del self.data_dict['z_3d']      

        # Initialize KDtree(s) 
        # > save the 2D one by default (initizalied during reader __init__()
        # > compute and save the time-varying 3D KDtree if relevant     
        
        logging.debug('saving reader''s 2D (horizontal) KDtree to ReaderBlockUnstruct')
        self.block_KDtree = KDtree # KDtree input to function = one computed during reader's __init__()

        # If we eventually use subset of nodes rather than full mesh, we'll need to re-compute the 2D KDtree
        # as well, instead of re-using the "full" one available from reader's init
        # logging.debug('Compute time-varying KDtree for 2D nearest-neighbor search') 
        # self.block_KDtree = cKDTree(np.vstack((self.x,self.y)).T)  # KDtree input to function = one computed during reader's __init__()

        if hasattr(self,'z_3d'):
            # we need to compute a new KDtree for that time step using vertical coordinates at that time step
            logging.debug('Compute time-varying KDtree for 3D nearest-neighbor search (i.e using ''zcor'') ')
            # clean arrays if needed, especially z_3d (get rid of nan's) - keep only non-nan
            # 
            # check for nan's 
            # if (self.z_3d != self.z_3d).any(): # not required
            #     self.x_3d = self.x_3d[np.where(self.z_3d == self.z_3d)]
            #     self.y_3d = self.y_3d[np.where(self.z_3d == self.z_3d)]
            #     self.z_3d = self.z_3d[np.where(self.z_3d == self.z_3d)]
            # check for infinite values
            if np.isinf(self.z_3d).any() :
                self.z_3d[np.where(np.isinf(self.z_3d))] = 15.0 #limit to +15.0m i.e. above msl

            self.block_KDtree_3d = cKDTree(np.vstack((self.x_3d,self.y_3d,self.z_3d)).T) 
            # do we need copy_data=True ..probably not since "data" [self.x_3d,self.y_3d,self.z_3d] 
            # will not change without the KDtree being recomputed

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
        
        # below probably not be relevant any longer
        if False:
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
            profiles_dict = {'z': self.z} # probably not valid...
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
            # standard data 2D or 3D
            else:
                # use KDtree to find nearest neighbours and interpolate based on distance, on 2D or 3D
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html#scipy.spatial.cKDTree.query
                # print(varname)
                # import pdb;pdb.set_trace()
                nb_closest_nodes = 3
                DMIN=1.e-10
                if data.shape[0] == self.x.shape[0] : # 2D data- full slice
                    #2D KDtree
                    dist,i=self.block_KDtree.query(np.vstack((x,y)).T,nb_closest_nodes, n_jobs=-1) #quick nearest-neighbor lookup
                    # dist = distance to nodes / i = index of nodes
                elif hasattr(self,'z_3d') and (data.shape[0] == self.x_3d.shape[0]) : #3D data
                    #3D KDtree
                    dist,i=self.block_KDtree_3d.query(np.vstack((x,y,z)).T,nb_closest_nodes, n_jobs=-1) #quick nearest-neighbor lookup
                    # dist = distance to nodes / i = index of nodes
                dist[dist<DMIN]=DMIN
                fac=(1./dist)
                data_interpolated = (fac*data.take(i)).sum(-1)/fac.sum(-1) 

     
                # # CHECK#################################################
                # data_KD = np.vstack((self.x_3d,self.y_3d,self.z_3d)).T
                # data_points_to_find = np.vstack((x,y,z)).T
                # data_KD[i[0][0],:]
                # data_points_to_find[0]
                # ########################################################
 
                # horizontal = self._interpolate_horizontal_layers(data, nearest=nearest)
            
            if profiles is not None and varname in profiles:
                # not functional yet...
                # need to lookup what actually is expected here
                profiles_dict[varname] = data_interpolated # horizontal

            # if horizontal.ndim > 1:
            #     env_dict[varname] = self.interpolator1d(data_interpolated) #self.interpolator1d(horizontal)
            # else:

            env_dict[varname] = data_interpolated #horizontal

        if 'z' in profiles_dict:
            profiles_dict['z'] = np.atleast_1d(profiles_dict['z'])

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