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
##########################################################################
# This reader supports native SCHISM output netcdf files and will 
# handle interpolation of current velocities and other variables
# in 2D and 3D.
# The interpolation is done using a cKDtree-approach defined
# form the 2D or 3D mesh nodes
# 
# Author: Simon Weppe. MetOcean Solution, MetService New Zealand
##########################################################################

##########################################################################
# This version aims to be more consistent with the new Opendrift release, 
# and re-uses some of the methods in unstructured.py, and better follow 
# the new logics
# 
# Note the interpolation approach is closer to what is done in 
# structured.py with cached reader "blocks" that are re-used for 
# successive interpolation between two model time steps. 
# unstructured.py instead read data at each interpolation cycle
# 
# 
# Ongoing work on speed improvements....
# 
#   - consider using pykdtree (https://github.com/storpipfugl/pykdtree)
#     for faster KDtree computations - test efficiency
# 
##########################################################################


import logging
logger = logging.getLogger(__name__)

import numpy as np
from datetime import datetime
from future.utils import iteritems
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree #cython-based KDtree for quick nearest-neighbor search
# consider using https://github.com/storpipfugl/pykdtree for KDtree computations - test efficiency
import pyproj
from opendrift.readers.basereader import BaseReader, UnstructuredReader
from opendrift.readers.basereader.consts import *
import xarray as xr
import shapely


class Reader(BaseReader,UnstructuredReader):

    def __init__(self, filename=None, name=None, proj4=None, use_3d = None , use_model_landmask = False,**kwargs):
        """Initialise reader_netCDF_CF_unstructured_SCHISM

        Args:
            filename    :   name of SCHISM netcdf file (can have wildcards)
            name        :   name of reader - optional, taken as filename if not input
                            o.readers['name']
            proj4       :   proj4 string defining spatial reference system. 
                            find string here : https://spatialreference.org/ref/epsg/
            use_3d      :   switch to use 3d flows (if available)
            use_model_landmask  : switch to use time-varying landmask from model wetdry_elem 
                                  (False by default)
            kwargs      : shoreline_file , allows adding a shoreline file that will be used alongside model 
                          landmask to flag particles on land in _get_variables_interpolated_()
                          Required in some cases to avoid particles going 
                          on land
        """
        if filename is None:
            raise ValueError('Need filename as argument to constructor')
        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        # Default interpolation method, see function interpolate_block()
        self.interpolation = 'linearNDFast'
        self.convolve = None  # Convolution kernel or kernel size
        
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
            'zcor' : 'vertical_levels', # time-varying vertical coordinates
            'sigma': 'ocean_s_coordinate',
            'vertical_velocity' : 'upward_sea_water_velocity',
            'wetdry_elem': 'land_binary_mask' }
            # diffusivity
            # viscosity

        self.return_block = True

        try:
            # Open file, check that everything is ok
            logger.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                logger.info('Opening files with open_mfdataset')
                self.dataset = xr.open_mfdataset(filename,chunks={'time': 1})
                # in case of issues with file ordering consider inputting an explicit filelist 
                # to reader, for example:
                # 
                # ordered_filelist_1 = []
                # ordered_filelist_1.extend(glob.glob(data_path + 'schism_marl2008*_00z_3D.nc'))# month block 1
                # ordered_filelist_1.sort()
            else:
                logger.info('Opening file with dataset')
                self.dataset = xr.open_dataset(filename,chunks={'time': 1})

        except Exception as e:
            raise ValueError(e)

        # Define projection of input data
        if proj4 is not None: #  user has provided a projection apriori
            self.proj4 = proj4
        else:                 # no input assumes latlon
            logger.debug('Lon and lat are 1D arrays, assuming latlong projection: proj4 = ''+proj=latlong''')
            self.proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' #'+proj=latlong'

        # check if 3d data is available and if we should use it
        self.use_3d = use_3d
        if self.use_3d is None : # not specified, use 3d data by default (if available)
            if 'hvel' in self.dataset.variables:
                self.use_3d = True
            else:
                self.use_3d = False

        if self.use_3d and 'hvel' not in self.dataset.variables:
            logger.debug('No 3D velocity data in file - cannot find variable ''hvel'' ')
        elif self.use_3d and 'hvel' in self.dataset.variables:
            if 'zcor' in self.dataset.variables: # both hvel and zcor in files - all good
                self.nb_levels = self.dataset.variables['hvel'].shape[2] #hvel dimensions : [time,node,lev,2]
            else:
                logger.debug('No vertical level information present in file ''zcor'' ... stopping')
                raise ValueError('variable ''zcor'' must be present in netcdf file to be able to use 3D currents')
        
        # set options for particles on-land check, using model landmask or not 
        self.use_model_landmask = use_model_landmask
        if self.use_model_landmask : 
            logger.debug('Using time-varying landmask from SCHISM for on-land particles checks (wetdry_elem variable) ')
            if 'shore_file' in kwargs:
                logger.debug('Also adding static shoreline file %s for additionnal on-land particles' %  kwargs['shore_file'])
                self.shore_file = kwargs['shore_file']
                self.shore_landmask = self.load_shoreline_landmask()
                        
        logger.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        for var_name in self.dataset.variables:

            if var_name in ['SCHISM_hgrid_face_x','SCHISM_hgrid_face_y','SCHISM_hgrid_edge_x','SCHISM_hgrid_edge_y']:
                # all these variables have the same standard name projection_x_coordinate,projection_y_coordinate
                # we need to only use :
                # SCHISM_hgrid_node_x as projection_x_coordinate
                # SCHISM_hgrid_node_y as projection_y_coordinate
                # which are the nodes (also called vertices)
                continue

            var = self.dataset.variables[var_name]

            attributes = var.attrs
            att_dict = var.attrs
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
                if var.ndim == 2:
                    # When datasets are concatenated by mfdataset(), coordinates vector (1D) may
                    # be tiled to a 2D array of size (time,node), keep only one vector for x,y  
                    var = var[0,:]
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                var_data = var.values
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
                if var.ndim == 2:
                    # When datasets are concatenated by mfdataset(), coordinates vector (1D) may
                    # be tiled to a 2D array of size (time,node), keep only one vector for x,y  
                    var = var[0,:]
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                else:
                    unitfactor = 1
                var_data = var.values
                y = var_data*unitfactor
                self.numy = var.shape[0]
            if standard_name == 'depth' or axis == 'Z':
                var_data = var.values
                if 'positive' not in var.ncattrs() or \
                        var.__dict__['positive'] == 'up':
                    self.z = var_data
                else:
                    self.z = -var_data
            if standard_name == 'time' or axis == 'T' or var_name == 'time':
                var_data = var.values
                # Read and store time coverage (of this particular file)
                time = var_data
                time_units = units
                # self.times = num2date(time, time_units)
                # convert from numpy.datetime64 to datetime
                self.times = [datetime.utcfromtimestamp((OT -
                    np.datetime64('1970-01-01T00:00:00Z')
                        ) / np.timedelta64(1, 's')) for OT in time]

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

        self.x = x
        self.y = y

        if not (self.x>360.).any() and self.use_3d :
            logger.debug('Native (x,y) coordinates are lon/lat - cannot use 3D data, setting use_3d = False')
            self.use_3d = False
            import pdb;pdb.set_trace()
            # the 3D interpolation currently doesnt work if the x,y coordinates in native netcdf files
            # are not cartesian but geopgraphic. If that is the case, when doing 3D interpolation and tree search, 
            # the vertical distance unit is meter, while the horizontal distance unit is degrees, which will return
            # erroneous "closest" nodes in ReaderBlockUnstruct,interpolate()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

        # compute CKDtree of (static) 2D nodes using _build_ckdtree_() from unstructured.py
        logger.debug('Building CKDtree of static 2D nodes for nearest-neighbor search')
        self.reader_KDtree = self._build_ckdtree_(self.x,self.y) 

        # build convex hull of points for particle-in-mesh checks using _build_boundary_polygon_() from unstructured.py
        logger.debug('Building convex hull of nodes for particle''s in-mesh checks')
        self.boundary = self._build_boundary_polygon_(self.x,self.y)

        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.dataset.variables:
            if var_name in [self.xname, self.yname]: #'depth'
                continue  # Skip coordinate variables
            var = self.dataset.variables[var_name]
            attributes = var.attrs
            att_dict = var.attrs

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

        self.xmin = self.x.min()
        self.xmax = self.x.max()
        self.ymin = self.y.min()
        self.ymax = self.y.max()
        # self.xmin = self.lon.min()
        # self.xmax = self.lon.max()
        # self.ymin = self.lat.min()
        # self.ymax = self.lat.max()

        # Run constructor of parent Reader class
        super(Reader, self).__init__()
        
        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present

    def build_ckdtree(self,x,y):
        # This is done using cython-based cKDTree from scipy for quick nearest-neighbor search
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
        # self.reader_KDtree = cKDTree(np.vstack((self.lon,self.lat)).T) 
        return cKDTree(np.vstack((x,y)).T) 

    def build_boundary_path(self,x,y):
        from scipy.spatial import ConvexHull # convex hull of grid points
        from matplotlib.path import Path # convex hull of grid points
        # 
        # Build a polygon of the boundary of the mesh, for in-grid checks
        # 
        # simple approximation : use convex hull of mesh nodes. This will correctly include the open boundary
        # (semi-circular) but it will miss the details of the shorelines. This is usually not a problem
        # since the shoreline details will be covered by the landmask.
        # 
        # **** Depreciated **** 
        # >> Now using prepared geometry rather than Matplotlib Path (as in unstructured .py)
        # https://sparkgeo.com/blog/using-prepared-geometries-in-shapely/

        hull_obj = ConvexHull(np.vstack([x,y]).T)
        hull_path = Path(np.vstack([x[hull_obj.vertices],y[hull_obj.vertices]]).T)  # Matplotlib Path object
        hull = np.vstack([x[hull_obj.vertices],y[hull_obj.vertices]]).T
        return hull_path # Matplotlib Path object
  
    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):

        """ The function extracts 'requested_variables' from the native SCHISM files
            which will then be used in _get_variables_interpolated_() to initialise the ReaderBlockUnstruct objects 
            used to interpolate data in space and time

            For now the function will extract the entire slice of data of 'requested_variables' at given 'time'

            There is an option to extract only a subset of data around particles clouds to have less data but
            it means we need to recompute the KDtree of the subset nodes every time in ReaderBlockUnstruct.
            
            Speed gain to be tested ...
            
        """
        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        
        variables = {'x': self.x, 'y': self.y, 'z': 0.*self.y,
                     'time': nearestTime}

        # extracts the full slices of requested_variables at time indxTime
        for par in requested_variables:
            if par not in ['x_sea_water_velocity','y_sea_water_velocity','land_binary_mask'] :
                # standard case - for all variables except current velocities
                var = self.dataset.variables[self.variable_mapping[par]]
                if var.ndim == 1:
                    data = var[:] # e.g. depth
                    logger.debug('reading constant data from unstructured reader %s' % (par))
                elif var.ndim == 2: 
                    data = var[indxTime,:] # e.g. 2D temperature
                    logger.debug('reading 2D data from unstructured reader %s' % (par))
                elif var.ndim == 3:
                    data = var[indxTime,:,:] # e.g. 3D salt [time,node,lev]
                    logger.debug('reading 3D data from unstructured reader %s' % (par))
                    # convert 3D data matrix to one column array and define corresponding data coordinates [x,y,z]
                    # (+ update variables dictionary with 3d coords if needed)
                    data,variables = self.convert_3d_to_array(indxTime,data,variables)
                else:
                    raise ValueError('Wrong dimension of %s: %i' %
                                     (self.variable_mapping[par], var.ndim))
            elif par in ['x_sea_water_velocity','y_sea_water_velocity'] :               
                # requested variables are current velocities
                # In SCHISM netcdf filesboth [u,v] components are saved 
                # as two different dimensions of the same variable.
                var = self.dataset.variables[self.variable_mapping[par]]
                if var.ndim == 3: # depth-averaged current data 'dahv' defined at each node and time [time,node,2]
                    if par == 'x_sea_water_velocity':
                       data = var[indxTime,:,0]
                    elif par == 'y_sea_water_velocity':
                       data = var[indxTime,:,1] 
                    logger.debug('reading 2D velocity data from unstructured reader %s' % (par))

                elif var.ndim == 4: # #3D current data 'hvel' defined at each node, level, and time [time,node,zcor,2]
                    if par == 'x_sea_water_velocity':
                       data = var[indxTime,:,:,0]   #hvel dimensions : [time,node,lev,2]
                    elif par == 'y_sea_water_velocity':
                       data = var[indxTime,:,:,1] #hvel dimensions : [time,node,lev,2]
                    logger.debug('reading 3D velocity data from unstructured reader %s' % (par))

                    # convert 3D data matrix to one column array and define corresponding data coordinates [x,y,z]
                    # (+ update variables dictionary with 3d coords if needed)
                    data,variables = self.convert_3d_to_array(indxTime,data,variables)

            elif (par in ['land_binary_mask']) & (self.use_model_landmask) :
                dry_elem = self.dataset.variables[self.variable_mapping[par]][indxTime,:] # dry_elem =1 if dry, 0 if wet, for each element face
                # find indices of nodes making up the dry elements
                if len(self.dataset['SCHISM_hgrid_face_nodes'].shape) == 3:
                    # when using open_mfdataset, the variable is expanded along time dimension, hence use of indxTime shape = (nb_time,nb_elem,4)
                    node_id = np.ravel(self.dataset['SCHISM_hgrid_face_nodes'][indxTime,dry_elem.values.astype('bool'),:]) # id of nodes making up each element 
                else:
                    node_id = np.ravel(self.dataset['SCHISM_hgrid_face_nodes'][dry_elem.values.astype('bool'),:]) # id of nodes making up each elements shape = (nb_elem,4)
                # remove nans, and keep unique indices
                node_id = np.unique(node_id[~np.isnan(node_id)]) - 1 # add <-1> to convert node index rather than absolute number
                data = 0.0*self.x # build array, same size as node
                data[node_id.astype(int)] = 1.0 # set dry nodes to one to be used as 'land_binary_mask'

            variables[par] = data # save all data slice to dictionary with key 'par'
            # Store coordinates of returned points
            #  >> done in previous steps here, line 380 and in convert_3d_to_array()
            variables[par] = np.asarray(variables[par])

        self.use_subset = False # Functionnal now  - need to check results are consistent.
        if self.use_subset: # for testing subsetting data before computin KD-trees
            variables = self.clip_reader_data(variables_dict = variables, x_particle = x,y_particle = y, requested_variables = requested_variables) # clip reader data to particle cloud coverage 
            # update the 2D KDtree (will be used to initialize the ReaderBlockUnstruct)
            self.reader_KDtree = cKDTree(np.vstack((variables['x'],variables['y'])).T) 
            # the 3D KDtree in updated within ReaderBlockUnstruct()
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
            vertical_levels = self.dataset.variables['zcor'][id_time,:,:]
            # depth are negative down consistent with convention used in OpenDrift 
            # if using the netCDF4 library, vertical_levels is masked array where "masked" levels are those below seabed  (= 9.9692100e+36)
            # if using the xarray library, vertical_levels is nan for levels are those below seabed

            # convert to masked array to be consistent with what netCDF4 lib returns
            vertical_levels = np.ma.array(vertical_levels, mask = np.isnan(vertical_levels.data)) 
            data = np.asarray(data)
            # vertical_levels.mask = np.isnan(vertical_levels.data) # masked using nan's when using xarray
        except:
            logger.debug('no vertical level information present in file ''zcor'' ... stopping')
            raise ValueError('variable ''zcor'' must be present in netcdf file to be able to use 3D currents')
        # flatten 3D data 
        data = np.ravel(data[~vertical_levels.mask])

        # add corresponding 3D coordinates to the 'variable_dict' which is eventually passed to get_variables_interpolated()
        # Note: these 3d coordinates will be the same for all 3d data (current,salt/temp etc..), and will change at each reader time step
        # They are saved as ['x_3d','y_3d','z_3d'] rather than ['x','y','z'] as it would break things when both 2d and 3d data are requested.
        if 'z_3d' not in variable_dict.keys():
            # now make a long 3-column array [lon,lat,z] [n_nodes x  3] at which 'data' is defined
            # tile lon/lat data so that array shape match vertical_levels shape
            x_tiled = np.tile(self.x,(self.nb_levels,1)).T # dimensions [node,vertical_levels]
            y_tiled = np.tile(self.y,(self.nb_levels,1)).T
            # lon_tiled = np.tile(self.lon,(self.nb_levels,1)).T # dimensions [node,vertical_levels]
            # lat_tiled = np.tile(self.lat,(self.nb_levels,1)).T
            # arrays are tiled so that lon_tiled[0,:] = return same value i.e. same lon/lat for all z levels 
            if x_tiled.shape != vertical_levels.shape:
                import pdb;pdb.set_trace()
            # convert to masked array consistent with vertical_levels
            x_tiled_ma = np.ma.array(x_tiled, mask = vertical_levels.mask) 
            y_tiled_ma = np.ma.array(y_tiled, mask = vertical_levels.mask)
            # flatten arrays and add to dictionary
            variable_dict['x_3d'] = np.ravel(x_tiled_ma[~vertical_levels.mask]) 
            variable_dict['y_3d'] = np.ravel(y_tiled_ma[~vertical_levels.mask])
            variable_dict['z_3d'] = np.ravel(vertical_levels[~vertical_levels.mask])

        return data,variable_dict

    def clip_reader_data(self,variables_dict, x_particle,y_particle,requested_variables=None):
        # clip "variables_dict" to current particle cloud extents [x_particle,y_particle] (+ buffer)
        # 
        # # Find a subset of mesh nodes that include the particle 
        # cloud. The frame should not be made too small to make sure
        # particles remain inside the frame over the time 
        # that interpolator is used (i.e. depends on model timestep)
        # 
        # 
        #  returns updated "variables" dict and data array
        deg2meters = 1.1e5 # approx
        buffer = .1  # degrees around given positions ~10km
        buffer = .05  # degrees around given positions ~ 5km - reader seems quite sensitive to that buffer...
        if self.xmax <= 360 :
            pass # latlong reference system - keep same buffer value
        else:
            buffer = buffer * deg2meters# projected reference system in meters
        self.subset_frame = [x_particle.min() - buffer,
                             x_particle.max() + buffer, 
                             y_particle.min() - buffer,
                             y_particle.max() + buffer]
        logger.debug('Spatial frame used for schism reader %s (buffer = %s) ' % (str(self.subset_frame) , buffer))
        # find lon,lat of reader that are within the particle cloud
        self.id_frame = np.where((self.lon >= self.subset_frame[0]) &
                     (self.lon <= self.subset_frame[1]) &
                     (self.lat >= self.subset_frame[2]) &
                     (self.lat <= self.subset_frame[3]))[0]
        logger.debug('Using %s ' % (int(100*self.id_frame.shape[0]/self.lon.shape[0])) + '%' + ' of native nodes' )
        # print(variables_dict['x'].shape)
        # print( self.id_frame.max())

        if False:
            import matplotlib.pyplot as plt
            plt.ion()
            plt.plot(x_particle,y_particle,'.')
            box = np.array([(self.subset_frame[0],self.subset_frame[2]),\
                            (self.subset_frame[1],self.subset_frame[2]),\
                            (self.subset_frame[1],self.subset_frame[3]),\
                            (self.subset_frame[0],self.subset_frame[3]),\
                            (self.subset_frame[0],self.subset_frame[2])])
            plt.plot(box[:,0],box[:,1])
            plt.title('frame used in clip_reader_data()')
            import pdb;pdb.set_trace()
            # plt.close()

        variables_dict['x'] = variables_dict['x'][self.id_frame]
        variables_dict['y'] = variables_dict['y'][self.id_frame]
        variables_dict['z'] = variables_dict['z'][self.id_frame]

        if 'x_3d' in variables_dict:
            ID= np.where((variables_dict['x_3d'] >= self.subset_frame[0]) &
                 (variables_dict['x_3d'] <= self.subset_frame[1]) &
                 (variables_dict['y_3d'] >= self.subset_frame[2]) &
                 (variables_dict['y_3d']<= self.subset_frame[3]))[0]  
            variables_dict['x_3d'] = variables_dict['x_3d'][ID]
            variables_dict['y_3d'] = variables_dict['y_3d'][ID]
            variables_dict['z_3d'] = variables_dict['z_3d'][ID]

        for par in requested_variables:
            if 'x_3d' not in variables_dict:
                # 2D data
                variables_dict[par] = variables_dict[par][self.id_frame]
            else: # can be either 2D or 3D data, in case requested_variables includes both 2D and 3D data
                if variables_dict[par].shape[0] == self.lon.shape[0] : # this "par" is 2D data
                    variables_dict[par] = variables_dict[par][self.id_frame]
                else:
                    # 3D data converted to array
                    variables_dict[par] = variables_dict[par][ID]

        return variables_dict
        
    def set_convolution_kernel(self, convolve):
        """Set a convolution kernel or kernel size (of array of ones) used by `get_variables` on read variables."""
        self.convolve = convolve

    def __convolve_block__(self, env):
        """
        Convolve arrays with a kernel, if reader.convolve is set
        """
        if self.convolve is not None:
            from scipy import ndimage
            N = self.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel / kernel.sum()
            else:
                kernel = N
            logger.debug('Convolving variables with kernel: %s' % kernel)
            for variable in env:
                if variable in ['x', 'y', 'z', 'time']:
                    pass
                else:
                    if env[variable].ndim == 2:
                        env[variable] = ndimage.convolve(env[variable],
                                                         kernel,
                                                         mode='nearest')
                    elif env[variable].ndim == 3:
                        env[variable] = ndimage.convolve(env[variable],
                                                         kernel[:, :, None],
                                                         mode='nearest')
        return env

    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   reader_x, reader_y, z):

        """
        This method _must_ be implemented by every reader. Usually by
        subclassing one of the reader types (e.g.
        :class:`structured.StructuredReader`).

        Arguments are in _native projection_ of reader.

        .. seealso:

            * :meth:`get_variables_interpolated_xy`.
            * :meth:`get_variables_interpolated`.
        """

        """
           Here, this function overloads the _get_variables_interpolated_() methods
           available in unstructured.py.  (which currently doesnt make use of blocks)

           The _get_variables_interpolated_() from structured.py uses regularly gridded 
           data "ReaderBlock" extracted from the netcdf files (which may possibly be "cached" 
           for speed improvements - see code for more detail).

           This function follows a similar approach but is instead using the native 
           high-resolution SCHISM data stored in "ReaderBlockUnstruct" which are used to 
           interpolate data in space and time. 

           The function returns environment data 'env' interpolated at particle positions [x,y] 

        """

        # block = False # legacy stuff 

        # Find reader time_before/time_after
        time_nearest, time_before, time_after, i1, i2, i3 = \
            self.nearest_time(time)
        logger.debug('Reader time:\n\t\t%s (before)\n\t\t%s (after)' %
                      (time_before, time_after))
        # For variables which are not time dependent, we do not care about time
        static_variables = ['sea_floor_depth_below_sea_level', 'land_binary_mask']

        if time == time_before or all(v in static_variables for v in variables):
            time_after = None

        if profiles is not None:
            # If profiles are requested for any parameters, we
            # add two fake points at the end of array to make sure that the
            # requested block has the depth range required for profiles
            mx = np.append(reader_x, [reader_x[-1], reader_x[-1]])
            my = np.append(reader_y, [reader_y[-1], reader_y[-1]])
            mz = np.append(z, [profiles_depth[0], profiles_depth[1]])
        else:
            mx = reader_x
            my = reader_y
            mz = z

        # z = z.copy()[ind_covered]  # Send values and not reference
        # EDIT:commented out: we need to keep the full array so that shapes are consistent
        
        # start interpolation procedure
        # 
        # general idea is to create a  new "ReaderBlockUnstruct" class that will be called instead of
        # the regular "ReaderBlock" as in basereader.py
        # This allows re-using almost the same code as structured.py for the block_before/block_after
        
        self.timer_start('preparing')        
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
            reader_data_dict = \
                 self.__convolve_block__(
                 self.get_variables(blockvariables_before, time_before,
                                    mx, my, mz))          
            # now use reader_data_dict to initialize 
            # a ReaderBlockUnstruct
            logger.debug('initialize ReaderBlockUnstruct var_block_before')
            self.var_block_before[blockvars_before] = \
                ReaderBlockUnstruct(reader_data_dict,
                    KDtree = self.reader_KDtree,
                    interpolation_horizontal=self.interpolation)
            try:
                len_z = len(self.var_block_before[blockvars_before].z)
            except:
                len_z = 1
            logger.debug(('Fetched env-block (size %ix%ix%i) ' +
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
                reader_data_dict = self.__convolve_block__(
                    self.get_variables(blockvariables_after, time_after, mx,
                                       my, mz))
                self.timer_start('preparing')
                logger.debug('initialize ReaderBlockUnstruct var_block_after')
                self.var_block_after[blockvars_after] = \
                    ReaderBlockUnstruct(
                        reader_data_dict,
                        KDtree = self.reader_KDtree,
                        interpolation_horizontal=self.interpolation)
                try:
                    len_z = len(self.var_block_after[blockvars_after].z)
                except:
                    len_z = 1

                logger.debug(('Fetched env-block (size %ix%ix%i) ' +
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
            logger.warning('Data block from %s not large enough to '
                            'cover element positions within timestep. '
                            'Buffer size (%s) must be increased.' %
                            (self.name, str(self.buffer)))
        self.timer_end('preparing') 

        ############################################################
        # Interpolate before/after blocks onto particles in space
        ############################################################
        self.timer_start('interpolation')
        logger.debug('Interpolating before (%s) in space  (%s)' %
                      (block_before.time, self.interpolation))
        env_before, env_profiles_before = block_before.interpolate(
                reader_x, reader_y, z, variables,
                profiles, profiles_depth)

        if (time_after is not None) and (time_before != time):
            logger.debug('Interpolating after (%s) in space  (%s)' %
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
            logger.debug(('Interpolating before (%s, weight %.2f) and'
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
                        logger.warning('Invalid values found for ' + var)
                        logger.warning(env[var][invalid])
                        logger.warning('(allowed range: [%s, %s])' %
                                        (standard_names[var]['valid_min'],
                                         standard_names[var]['valid_max']))
                        logger.warning('Replacing with NaN')
                        env[var][invalid] = np.nan
            # Interpolating vertical profiles in time
            if profiles is not None:
                env_profiles = {}
                logger.info('Interpolating profiles in time')
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
            logger.debug('No time interpolation needed - right on time.')
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
        # the masking, rotation etc.. is done in variables.py get_variables_interpolated_xy()

        # apply log profile if we are interpolating 2D data for ['x_sea_water_velocity','y_sea_water_velocity']
        # not using for now
        if False :
            self.apply_logarithmic_current_profile(env,z)
         
        # additional on-land checks using shore_landmask (if present)
        if 'land_binary_mask' in env.keys() and self.use_model_landmask and hasattr(self,'shore_file'):
            logger.debug('Updating land_binary_mask using shoreline landmask <%s> ' % self.shore_file)
            lon_tmp,lat_tmp = self.xy2lonlat(reader_x,reader_y)
            on_shore_landmask = shapely.vectorized.contains(self.shore_landmask, lon_tmp, lat_tmp)
            # update the 'land_binary_mask' accounting for shoreline landmask
            env['land_binary_mask'] = np.maximum(env['land_binary_mask'],on_shore_landmask.astype(float))

        # make sure dry points have zero velocities which is not always the case
        # we could also look at using depth and thresholds to flag other dry points ?
        if 'land_binary_mask' in env.keys():
            logger.debug('Setting [x_sea_water_velocity,y_sea_water_velocity] to zero at dry points')
            env['x_sea_water_velocity'][env['land_binary_mask'].astype('bool')] = 0
            env['y_sea_water_velocity'][env['land_binary_mask'].astype('bool')] = 0

        return env, env_profiles

    def covers_positions_xy(self, x, y, z=0):
        """
        Check which points are within boundary of mesh.

        Wrapper function of covers_positions() from unstructured.py which is called in 
        get_variables_interpolated_xy() function from variables.py 

        It returns indices of in-mesh points, and in-mesh point coordinates rather than a boolean array (inside/outside) 

        Within get_variables_interpolated_xy() from variables.py, data is queried for these in-mesh points only and the 
        full array (incl. out of mesh positions) is re-generated with correct masking 

        """
        ind_covered = np.where(self.covers_positions(x, y))[0]
        return ind_covered ,x[ind_covered], y[ind_covered]

    def apply_logarithmic_current_profile(self,env,z):
        if not self.use_3d and 'sea_floor_depth_below_sea_level' in self.variables and 'x_sea_water_velocity' in self.variables :
            log_profile_factor = self.logarithmic_current_profile(z,env['sea_floor_depth_below_sea_level'])
            logger.debug('Applying logarithmic current profile to 2D current data [x_sea_water_velocity,y_sea_water_velocity] %s <= factor <=%s' % (np.min(log_profile_factor), np.max(log_profile_factor) ))
            env['x_sea_water_velocity'] = log_profile_factor * env['x_sea_water_velocity']
            env['y_sea_water_velocity'] = log_profile_factor * env['y_sea_water_velocity']
            if False:
                import matplotlib.pyplot as plt
                plt.ion()
                plt.plot(z/env['sea_floor_depth_below_sea_level'],log_profile_factor,'.')
                import pdb;pdb.set_trace()
                plt.close()

    def logarithmic_current_profile(self, particle_z, total_depth):
        ''' 
        Extrapolation of depth-averaged currents to any vertical 
        level of the water column assuming a logarithmic profile


        Inputs :
            particle_z : vertical position of particle in water column (negative down as used in Opendrift)
            total_depth : total water depth at particle position (positive down)
            z0 : roughness length, in meters (default, z0 = 0.001m )

        Returns : 
            Factors to be apply to interpolated raw depth-averaged currents

        Reference :
            Van Rijn, 1993. Principles of Sediment Transport in Rivers,
            Estuaries and Coastal Seas

        '''

        # Opendrift convention : particle_z is 0 at the surface, negative down
        # 
        # The particle_z we need is the height of particle above seabed (positive)
        part_z_above_seabed = np.abs(total_depth) + particle_z 
        # note : taking the absolute value enbsure we have positive down depth (though it make any depth<0 wrong..but log profile probably not critical at these points anyway?)
        # if we are sure that total_depth is positive down then we should just use >> part_z_above_seabed = total_depth + particle_z 
        if not hasattr(self,'z0'): 
            self.z0 = 0.001 # typical value for sandy seabed
        log_fac = ( np.log(part_z_above_seabed / self.z0) ) / ( np.log(np.abs(total_depth)/self.z0)-1 ) # total_depth must be positive, hence the abs()
        log_fac[np.where(part_z_above_seabed<=0)] = 1.0 # do not change velocity value
        return log_fac

        shore = np.loadtxt(
            polygon_file)  # nan-delimited closed polygons for land and islands
    
    def load_shoreline_landmask(self):
        # load shoreline polygon is added by user,used for additional on-land checks
        #  must be in wgs84
        # see reader_landmask_custom.py
        from shapely.geometry import Polygon, MultiPolygon, asPolygon
        import shapely

        shore = np.loadtxt(
            self.shore_file)  # nan-delimited closed polygons for land and islands

        # Loop through polygon to build MultiPolygon object
        #
        # make sure that start and end lines are [nan,nan] as well
        if not np.isnan(shore[0, 0]):
            shore = np.vstack(([np.nan, np.nan], shore))
        if not np.isnan(shore[-1, 0]):
            shore = np.vstack((shore, [np.nan, np.nan]))
        id_nans = np.where(np.isnan(shore[:, 0]))[0]
        poly = []
        for cnt, _id_i in enumerate(id_nans[:-1]):
            # The shapely.geometry.asShape() family of functions can be used to wrap Numpy coordinate arrays
            # https://shapely.readthedocs.io/en/latest/manual.html
            poly.append(
                asPolygon(shore[id_nans[cnt] + 1:id_nans[cnt + 1] - 1, :]))
        # We can pass multiple Polygon -objects into our MultiPolygon as a list
        landmask = MultiPolygon(poly)
        landmask = shapely.prepared.prep(landmask)
        return landmask

    def plot_mesh(self, variable=None, vmin=None, vmax=None,
             filename=None, title=None, buffer=1, lscale='auto',plot_time = None):
        """Plot geographical coverage of reader."""

        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from opendrift_landmask_data import Landmask
        import matplotlib.tri as mtri
        fig = plt.figure()

        if plot_time is None : 
            plot_time = self.start_time


        #####################################
        # In Dev - plot unstructured mesh
        # 
        # To do 
        #  - plot in correct cartopy projection, consistent with basereader.py plot()
        #  - plot a given timestep for flows
        #  - incorporate that to animation as in basemodel.py etc...
        # 
        #####################################
        if False:
            import pdb;pdb.set_trace()
            import matplotlib.tri as mtri
            from mpl_toolkits.mplot3d import Axes3D
            X=self.dataset.variables['SCHISM_hgrid_node_x'][:]
            Y=self.dataset.variables['SCHISM_hgrid_node_y'][:]
            Z=self.dataset.variables['depth'][:]
            face=self.dataset.variables['SCHISM_hgrid_face_nodes'][:,0:3]-1
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

        if variable is None: # generic plot to check extents
            boundary = Polygon(list(zip(xsp, ysp)), alpha=0.5, ec='k', fc='b',
                               zorder=100)
            ax.add_patch(boundary)
            # add nodes
            rx = np.array([self.xmin, self.xmax])
            ry = np.array([self.ymin, self.ymax])
            data = self.get_variables('sea_floor_depth_below_sea_level', plot_time,rx, ry, block=True)
            rlon, rlat = self.xy2lonlat(data['x'], data['y']) # node coordinates
            ax.plot(rlon,rlat, '.',transform=ccrs.PlateCarree())
            # set extents
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
            if len(variable) == 1: # simple scalar variable
                variable = variable[0]
                rx = np.array([self.xmin, self.xmax])
                ry = np.array([self.ymin, self.ymax])
                data = self.get_variables(variable, plot_time,rx, ry, block=True)
                rlon, rlat = self.xy2lonlat(data['x'], data['y']) # node coordinates
                data[variable] = np.ma.masked_invalid(data[variable])
                # ax.plot(rlon,rlat, '.',transform=ccrs.PlateCarree())
                face=self.dataset.variables['SCHISM_hgrid_face_nodes'][:,0:3]-1
                # build triangulation
                triang =mtri.Triangulation(rlon,rlat, triangles=face, mask=None)
                # plot the variable 
                ax.tricontourf(triang, data[variable],vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
                # add the triangular mesh - too slow
                # ax.triplot(triang, transform=ccrs.PlateCarree()) 

            elif len(variable)>1 and all('water' in var for var in variable):  # vector field variable = ['x_sea_water_velocity','y_sea_water_velocity']
                rx = np.array([self.xmin, self.xmax])
                ry = np.array([self.ymin, self.ymax])
                data = self.get_variables(variable, plot_time,rx, ry, block=True)
                rlon, rlat = self.xy2lonlat(data['x'], data['y']) # node coordinates
                # data[variable] = np.ma.masked_invalid(data[variable])
                # ax.plot(rlon,rlat, '.',transform=ccrs.PlateCarree())
                face=self.dataset.variables['SCHISM_hgrid_face_nodes'][:,0:3]-1
                # build triangulation
                triang =mtri.Triangulation(rlon,rlat, triangles=face, mask=None)
                if False:
                    # plot the variable 
                    ax.tricontourf(triang, data[variable],vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
                    # add the triangular mesh - too slow
                    ax.triplot(triang, transform=ccrs.PlateCarree()) 

                ax.quiver(rlon,rlat,data['x_sea_water_velocity'],data['y_sea_water_velocity'],np.sqrt(data['x_sea_water_velocity']**2 + data['y_sea_water_velocity']**2), transform=ccrs.PlateCarree() )

            ax.set_extent([xsp.min()-.5, xsp.max()+.5, ysp.min()-.5, ysp.max()+.5], crs=sp)

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.ion()
            plt.show()
            # import pdb;pdb.set_trace()
            # 
            # variable = ['sea_floor_depth_below_sea_level','x_sea_water_velocity','y_sea_water_velocity']
            # data = self.get_variables(variable, self.start_time,rx, ry, block=True) # where variable = ['x_sea_water_velocity','y_sea_water_velocity']
            # plt.quiver(rlon,rlat,data['x_sea_water_velocity']) #> check what's going here
            # >>> try to add quiver to check flow fields and how they are interpolated over time
            # >>> also consider adding the possibiltyy to overlay quiver+particle motions for the reader_netCDF_CF_generic structured

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

           KDtree : for nearest-neighbor search (initialized using SCHISM nodes in reader's _init_() )
                    This is read from reader object, so that it is not recomputed every time

    """
    logger = logging.getLogger('opendrift')  # using common logger

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
            # probably not valid ..since z are different for each point in SCHISM..rather than fixed
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
        
        logger.debug('saving reader''s 2D (horizontal) KDtree to ReaderBlockUnstruct')
        self.block_KDtree = KDtree # KDtree input to function = one computed during reader's __init__()

        # If we eventually use subset of nodes rather than full mesh, we'll need to re-compute the 2D KDtree
        # as well, instead of re-using the "full" one available from reader's init
        # logger.debug('Compute time-varying KDtree for 2D nearest-neighbor search') 
        # self.block_KDtree = cKDTree(np.vstack((self.x,self.y)).T)  # KDtree input to function = one computed during reader's __init__()

        if hasattr(self,'z_3d'):
            # we need to compute a new KDtree for that time step using vertical coordinates at that time step
            logger.debug('Compute time-varying KDtree for 3D nearest-neighbor search (i.e using ''zcor'') ')
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
            # will not change without the KDtree being recomputedplt

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
                logger.warning('Ensemble data currently not extrapolated towards seafloor')
            elif self.data_dict[var].ndim == 3:
                filled = fill_NaN_towards_seafloor(self.data_dict[var])
                if filled is True:
                    filled_variables.add(var)
                
        if len(filled_variables) > 0:
            logger.debug('Filled NaN-values toward seafloor for :'
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
            logger.debug('Nearest interpolation will be used '
                          'for landmask, and %s for other variables'
                          % interpolation_horizontal)

    def _initialize_interpolator(self, x, y, z=None):
        logger.debug('Initialising interpolator.')
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
                # here we need to make a choice on when to flag as land 
                if False :
                    nearest = True
                    self.interpolator2d_nearest = Nearest2DInterpolator(self.x, self.y, x, y)
                # if closest node is dry then we assume particle is on land
                nb_closest_nodes = 1
                #2D KDtree
                dist,i=self.block_KDtree.query(np.vstack((x,y)).T,nb_closest_nodes, n_jobs=-1) #quick nearest-neighbor lookup
                data_interpolated = data[i] # we keep closest value

            # ensemble data
            if type(data) is list:
                num_ensembles = len(data)
                logger.debug('Interpolating %i ensembles for %s' % (num_ensembles, varname))
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
                    ##############################
                    # PLOT CHECKS
                    if False:
                        import matplotlib.pyplot as plt
                        fig = plt.figure()
                        ax = fig.add_subplot(111, projection='3d')
                        ax.scatter(self.x_3d[i[0]].data.tolist(),self.y_3d[i[0]].data.tolist(),self.z_3d[i[0]].data.tolist(),c='r', marker='o')
                        ax.scatter(x[0:1],y[0:1],z[0:1],c='g', marker='o')
                        plt.ion()
                        plt.show()
                    ##############################3

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

        if len(indices) != len(x):
            import matplotlib.pyplot as plt
            plt.ion()
            plt.plot(x,y,'r.')
            box = np.array([(self.x.min(),self.y.min()),\
            (self.x.max(),self.y.min()),\
            (self.x.max(),self.y.max()),\
            (self.x.min(),self.y.max()),\
            (self.x.min(),self.y.min())])
            plt.plot(box[:,0],box[:,1],'k--')
            plt.title('Increase buffer distance around particle cloud')
            import pdb;pdb.set_trace()
            plt.close()

        if len(indices) == len(x):
            return True
        else:
            return False


    # def covers_positions(self, lon, lat, z=0):
    #     """Return indices of input points covered by reader.
        
    #     For an unstructured reader, this is done by checking that if particle positions are within the convex hull
    #     of the mesh nodes. This means it is NOT using the true polygon bounding the mesh ..but this is generally good enough
    #     to locate the outer boundary (which is often semi-circular). 
    #     >> The convex hull will not be good for the shorelines though, but this is not critical since coast interaction 
    #     will be handled by either by the landmask (read from file or using global_landmask) 
    #     Data interpolation might be an issue for particles reaching the land, and not actually flagged as out-of-bounds...
    #     To Check...

    #     Better alternatives could be: 
    #         - get the true mesh polygon from netCDF files..doesnt seem to be available.
    #         - compute an alpha-shape of mesh nodes (i.e. enveloppe) and use that as polygon. better than convexhull but not perfect 
    #                 (shp = alphaShape(double(x),double(y),10000); work well in matlab for example)
    #         - workout that outer mesh polygon from a cloud of points...maybe using a more evolved trimesh package ?

    #     """        
    #     # weird bug in which x,y become 1e30 sometimes...
    #     # need to print stuff to make it work...

    #     # Calculate x,y coordinates from lon,lat
    #     print(lon.max()) # if I comment this...[x,y] may become 1e+30..have no idea why
    #                      # it runs fine if I have that print statement ....
    #     # print(lat.max())
    #     x, y = self.lonlat2xy(lon, lat)      
    #     # Only checking vertical coverage if zmin, zmax is defined
    #     zmin = -np.inf
    #     zmax = np.inf
    #     if hasattr(self, 'zmin') and self.zmin is not None:
    #         zmin = self.zmin
    #     if hasattr(self, 'zmax') and self.zmax is not None:
    #         zmax = self.zmax

    #     # sometimes x,y will are = 1e30 at this stage..and have no idea why
    #     # recomputing them below seems to fix the issue        
    #     # x, y = self.lonlat2xy(lon, lat)

    #     if self.global_coverage():
    #         pass
    #         # unlikely to be used for SCHISM domain
    #     else:
    #         # Option 1 : use the Path object (matplotlib) defined in __init__() : self.hull_path
    #         # in_hull =self.hull_path.contains_points(np.vstack([x,y]).T)
    #         # Option 2 : use the Prepared Polygon object
    #         in_hull = UnstructuredReader.covers_positions(self,x, y, z)
    #         indices = np.where(in_hull) 
    #     try:
    #         return indices, x[indices], y[indices]
    #     except:
    #         return indices, x, y    