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
# This reader supports native output netcdf files from the WindWaveModelIII 
# used in the coupled SCHISM-WWM III wave-flow model. 
# 
# The interpolation of variable to particle positions is done 
# using a cKDtree-approach defined from the 2D nodes as the SCHISM reader
# 
# Author: Simon Weppe. MetOcean Solution, MetService New Zealand
# 
# Roland, A., Zhang, Y.J., Wang, H.V., Meng, Y., Teng, Y.C., Maderich, V., Brovchenko, I.,
# Dutour-Sikiric, M., Zanke, U., 2012. A fully coupled 3D wave-current interaction
# model on unstructured grids. J. Geophys. Res. 117, C00J33. http://dx.doi.org/10.
# 1029/2012JC007952.
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
##########################################################################


import logging
logger = logging.getLogger(__name__)

import numpy as np
from datetime import datetime
from future.utils import iteritems
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree #cython-based KDtree for quick nearest-neighbor search
import pyproj
from opendrift.readers.basereader import BaseReader, UnstructuredReader
from opendrift.readers.reader_schism_native import Reader as SchismReader
from opendrift.readers.reader_schism_native import ReaderBlockUnstruct 
from opendrift.readers.basereader.consts import *


try:
    import xarray as xr
    has_xarray = True
except:
    has_xarray = False


class Reader(SchismReader): # >> inherits from SCHISM reader

    def __init__(self, filename=None, name=None, proj4=None):
        """Initialise reader_netCDF_CF_unstructured_SCHISM

        Args:
            filename:  name of WWM netcdf file (can have wildcards)
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

        # Default interpolation method, see function interpolate_block()
        self.interpolation = 'linearNDFast'
        self.convolve = None  # Convolution kernel or kernel size
        
        # [name_used_in_schism : equivalent_CF_name]
        WWM_mapping = {
            'depth': 'sea_floor_depth_below_sea_level',
            'WLEV' : 'sea_surface_height',
            'HS'    : 'sea_surface_wave_significant_height',
            'STOKESBAROX' : 'sea_surface_wave_stokes_drift_x_velocity',
            'STOKESBAROY' : 'sea_surface_wave_stokes_drift_y_velocity',
            'TPPD' : 'sea_surface_wave_period_at_variance_spectral_density_maximum'}
            # TPPD 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
            # diffusivity
            # viscosity

            # More variables available :
            # pick names from CF tables here : http://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
            # 
            # HS    'Significant wave height'
            # DM    'Mean wave direction'
            # DSPR   'Directional spreading'
            # TPPD   'Discrete peak wave period'
            # CPP    'peak wave speed'
            # CGPP   'peak group velocity'
            # PEAKDSPR 'Peak directional spreading'
            # ORBITAL   'orbital velocity'
            # CD    'drag coefficient'
            # WATLEV  'water level'
            # DEPDT  'water level change'
            # STOKESBAROX    'barotropic Stokes velocity in X direction'
            # STOKESBAROY    'barotropic Stokes velocity in Y direction'
            # RSXX    'barotropic Stress potential Sxx'
            # RSXY   'barotropic Stress potential Sxy'
            # RSYY  'barotropic Stress potential Syy'
         
        # >> ONGOING WORK
        # load correct coordinates lon,lat
        # coordinate/time variables names to update
        # check interpolation is done correctly


        self.return_block = True

        try:
            # Open file, check that everything is ok
            logger.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                logger.info('Opening files with MFdataset')
                # self.dataset = MFdataset(filename)
                if has_xarray:
                    self.dataset = xr.open_mfdataset(filename,chunks={'ocean_time': 1})
                else:
                    self.dataset = MFDataset(filename,aggdim='ocean_time')
                # in case of issues with file ordering consider inputting an explicit filelist 
                # to reader, for example:
                # 
                # ordered_filelist_1 = []
                # ordered_filelist_1.extend(glob.glob(data_path + 'schism_marl2008*_00z_3D.nc'))# month block 1
                # ordered_filelist_1.sort()
            else:
                logger.info('Opening file with dataset')
                # self.dataset = dataset(filename, 'r')
                if has_xarray:
                    self.dataset = xr.open_dataset(filename,chunks={'ocean_time': 1})
                else:
                    self.dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

        # Define projection of input data
        if proj4 is not None: #  user has provided a projection apriori
            self.proj4 = proj4
        else:                 # no input assumes latlon
            logger.debug('Lon and lat are 1D arrays, assuming latlong projection: proj4 = ''+proj=latlong''')
            self.proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' #'+proj=latlong'
                        
        logger.debug('Finding coordinate variables.')
        # Find x, y and z coordinates
        for var_name in self.dataset.variables:

            var = self.dataset.variables[var_name]
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
            else:
                standard_name = var_name # no standard_name attribute in files in current version, use variable name
            if 'long_name' in attributes:
                long_name = att_dict['long_name']
            if 'axis' in attributes:
                axis = att_dict['axis']
            if 'units' in attributes:
                units = att_dict['units']
            if '_CoordinateAxisType' in attributes:
                CoordinateAxisType = att_dict['_CoordinateAxisType']

            if standard_name == 'lon':
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
            if standard_name == 'lat':
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

                # this is actually levels here, rather than depth 
                # it should now  set self.z using bathy data 

                # if 'positive' not in var.attrs or \
                #         var.__dict__['positive'] == 'up':
                #     self.z = var_data
                # else:
                #     self.z = -var_data
            if standard_name == 'ocean_time' or axis == 'T' or var_name == 'time':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                # Read and store time coverage (of this particular file)
                time = var_data
                time_units = units
                # self.times = num2date(time, time_units)
                if has_xarray:
                    # convert from numpy.datetime64 to datetime
                    self.times = [datetime.utcfromtimestamp((OT -
                        np.datetime64('1970-01-01T00:00:00Z')
                            ) / np.timedelta64(1, 's')) for OT in time]
                else:
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

        self.x = x
        self.y = y
        
        # compute CKDtree of (static) 2D nodes using _build_ckdtree_() from unstructured.py
        logger.debug('Building CKDtree of static 2D nodes for nearest-neighbor search')
        # self.reader_KDtree = self.build_ckdtree(self.lon,self.lat)
        self.reader_KDtree = UnstructuredReader._build_ckdtree_(self,self.x,self.y) 

        # build convex hull of points for particle-in-mesh checks using _build_boundary_polygon_() from unstructured.py
        logger.debug('Building convex hull of nodes for particle''s in-mesh checks')
        # self.hull_path = self.build_boundary_path(x, y)
        self.boundary = UnstructuredReader._build_boundary_polygon_(self,self.x,self.y)
        # self.hull_path = self.build_boundary_path(self.x,self.y) # matplotlib Path object

        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.dataset.variables:
            if var_name in [self.xname, self.yname]: #'depth'
                continue  # Skip coordinate variables
            var = self.dataset.variables[var_name]

            if has_xarray:
                attributes = var.attrs
                att_dict = var.attrs
            else:
                attributes = var.ncattrs()
                att_dict = var.__dict__

            if var_name in WWM_mapping:     
                # build the mapping between WWM and Opendrfit variable names         
                self.variable_mapping[WWM_mapping[var_name]] = str(var_name) 
                   
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
        super(SchismReader, self).__init__()
        
        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present


    def plot(self, variable=None, vmin=None, vmax=None,
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
            if len(variable) == 1: # simple scalar variable - plot colored surface
                variable = variable[0]
                rx = np.array([self.xmin, self.xmax])
                ry = np.array([self.ymin, self.ymax])
                data = self.get_variables(variable, plot_time,rx, ry, block=True)
                rlon, rlat = self.xy2lonlat(data['x'], data['y']) # node coordinates
                data[variable] = np.ma.masked_invalid(data[variable])
                # ax.plot(rlon,rlat, '.',transform=ccrs.PlateCarree())
                face=self.dataset.variables['ele'][:,0:3]-1
                # build triangulation
                triang =mtri.Triangulation(rlon,rlat, triangles=face, mask=None)
                # plot the variable 
                ax.tricontourf(triang, data[variable].filled(fill_value = 0.0),vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
                # add the triangular mesh - too slow
                # ax.triplot(triang, transform=ccrs.PlateCarree()) 

            elif len(variable) == 2:  # (u,v) components of a vector field - plot colored quiver
                rx = np.array([self.xmin, self.xmax])
                ry = np.array([self.ymin, self.ymax])
                data = self.get_variables(variable, plot_time,rx, ry, block=True)
                rlon, rlat = self.xy2lonlat(data['x'], data['y']) # node coordinates
                # data[variable] = np.ma.masked_invalid(data[variable])
                # ax.plot(rlon,rlat, '.',transform=ccrs.PlateCarree())
                face=self.dataset.variables['ele'][:,0:3]-1
                # build triangulation
                triang =mtri.Triangulation(rlon,rlat, triangles=face, mask=None)
                if False:
                    # plot the variable 
                    ax.tricontourf(triang, data[variable],vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
                    # add the triangular mesh - too slow
                    ax.triplot(triang, transform=ccrs.PlateCarree()) 

                ax.quiver(rlon,rlat,data[variable[0]],data[variable[1]],np.sqrt(data[variable[0]]**2 + data[variable[1]]**2), transform=ccrs.PlateCarree() )

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
            import pdb;pdb.set_trace()