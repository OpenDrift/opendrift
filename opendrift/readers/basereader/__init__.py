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
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# Copyright 2020, Gaute Hope, MET Norway

import sys
import logging
logger = logging.getLogger(__name__)
from abc import abstractmethod, ABCMeta
from datetime import datetime, timedelta

import numpy as np
import pyproj
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from .structured import StructuredReader
from .unstructured import UnstructuredReader
from .continuous import ContinuousReader
from .variables import Variables
from .consts import *
from ..operators.ops import Combine, Filter

from opendrift.readers.interpolation import ReaderBlock

class BaseReader(Variables, Combine, Filter):
    """
    An abstract reader. Implementors provide a method to read data and specify how it is interpolated.

    This class inherits :class:`.variables.Variables` which inherits :class:`.variables.ReaderDomain`. `ReaderDomain` is responsible for the extent and domain of the reader, including checking for out-of-bounds and projection conversion. `Variables` is responsible for returning interpolated data at the requests positions or profiles. Apart from coercing the returned data into the right type for :py:mod:`opendrift.models.basemodel`, it defines the abstract interface to :meth:`.variables.Variables._get_variables_interpolated_` which reader-implementations must provide (_usually_ through one of the main reader-types, see: :py:mod:`opendrift.readers`).

    .. seealso::

        :py:mod:`opendrift.readers`.

        :py:mod:`.variables`.
    """

    __metaclass__ = ABCMeta

    verticalbuffer = 1  # To be overridden by application as needed

    # Mapping variable names, e.g. from east-north to x-y, temporarily
    # presuming coordinate system then is lon-lat for equivalence
    variable_aliases = {
        'sea_water_potential_temperature': 'sea_water_temperature',
        'x_wind_10m': 'x_wind',
        'y_wind_10m': 'y_wind',
        'surface_wind_speed': 'wind_speed',
        'surface_wind_direction': 'wind_from_direction',
        'wind_direction': 'wind_from_direction',
        'sea_water_x_velocity': 'x_sea_water_velocity',
        'sea_water_y_velocity': 'y_sea_water_velocity',
        'baroclinic_x_sea_water_velocity': 'x_sea_water_velocity',
        'baroclinic_y_sea_water_velocity': 'y_sea_water_velocity',
        'x_sea_ice_velocity': 'sea_ice_x_velocity',
        'y_sea_ice_velocity': 'sea_ice_y_velocity',
        'barotropic_sea_water_x_velocity': 'sea_ice_x_velocity',
        'barotropic_sea_water_y_velocity': 'sea_ice_y_velocity',
        'salinity_vertical_diffusion_coefficient' : 'ocean_vertical_diffusivity',
        'ocean_mixed_layer_thickness_defined_by_sigma_theta': 'ocean_mixed_layer_thickness',
        'sea_floor_depth_below_sea_surface' : 'sea_floor_depth_below_sea_level',
        'sea_floor_depth_below_geoid' : 'sea_floor_depth_below_sea_level',
        'mass_concentration_of_suspended_matter_in_sea_water' : 'spm'
        }

    xy2eastnorth_mapping = {
        'x_sea_water_velocity': ['eastward_sea_water_velocity', 'surface_eastward_sea_water_velocity',
                                 'baroclinic_eastward_sea_water_velocity',
                                 'eastward_current_velocity', 'eastward_tidal_current',
                                 'eastward_ekman_current_velocity', 'eastward_geostrophic_current_velocity',
                                 'eastward_eulerian_current_velocity', 'surface_geostrophic_eastward_sea_water_velocity',
                                 'surface_geostrophic_eastward_sea_water_velocity_assuming_sea_level_for_geoid',
                                 'surface_eastward_geostrophic_sea_water_velocity_assuming_sea_level_for_geoid'],
        'y_sea_water_velocity': ['northward_sea_water_velocity', 'surface_northward_sea_water_velocity',
                                 'baroclinic_northward_sea_water_velocity',
                                 'northward_current_velocity', 'northward_tidal_current',
                                 'northward_ekman_current_velocity', 'northward_geostrophic_current_velocity',
                                 'northward_eulerian_current_velocity', 'surface_geostrophic_northward_sea_water_velocity',
                                 'surface_geostrophic_northward_sea_water_velocity_assuming_sea_level_for_geoid',
                                 'surface_northward_geostrophic_sea_water_velocity_assuming_sea_level_for_geoid'],
        'x_wind': 'eastward_wind', 'y_wind': 'northward_wind',
        'sea_surface_wave_stokes_drift_x_velocity': 'eastward_surface_stokes_drift',
        'sea_surface_wave_stokes_drift_y_velocity': 'northward_surface_stokes_drift'}

    def __init__(self):
        """Common constructor for all readers"""
        super().__init__()

        self.number_of_fails = 0  # Readers may be quanrantined after a number of fails (config setting)

        self.always_valid = False  # Set to True if a single field should
                                   # be valid at all times

        self.is_lazy = False  # Generally False

        # Set projection for coordinate transformations
        if self.proj is None:
            if self.proj4 is not None:
                try:
                    self.proj = pyproj.Proj(self.proj4)
                except Exception as ex:
                    logger.exception(ex)
                    # Workaround for proj-issue with zero flattening:
                    # https://github.com/OSGeo/proj.4/issues/1191
                    origproj4 = self.proj4
                    self.proj4 = self.proj4.replace('+e=0.0', '')
                    self.proj4 = self.proj4.replace('+e=0', '')
                    self.proj4 = self.proj4.replace('+f=0.0', '')
                    self.proj4 = self.proj4.replace('+f=0', '')
                    if origproj4 != self.proj4:
                        logger.warning('Removing flattening parameter from proj4; %s -> %s' % (origproj4, self.proj4))
                    self.proj = pyproj.Proj(self.proj4)
            else:
                raise Exception("No projection or CRS information available")

        # Check if there are holes in time domain
        if self.start_time is not None and self.time_step is not None:# and len(self.times) > 1:
            self.expected_time_steps = (
                self.end_time - self.start_time).total_seconds() / (
                self.time_step.total_seconds()) + 1
            if self.times is not None:
                self.missing_time_steps = self.expected_time_steps - \
                    len(self.times)
            else:
                self.missing_time_steps = 0
            self.actual_time_steps = self.expected_time_steps - \
                self.missing_time_steps

        # Making sure start_time is datetime, and not cftime object
        if self.start_time is not None:
             self.start_time = datetime(self.start_time.year, self.start_time.month,
                                   self.start_time.day, self.start_time.hour,
                                   self.start_time.minute, self.start_time.second)

        # Calculate shape (size) of domain
        try:
            numx = (self.xmax - self.xmin)/self.delta_x + 1
            numy = (self.ymax - self.ymin)/self.delta_y + 1
            self.shape = (int(numx), int(numy))
        except:
            self.shape = None

        self.set_buffer_size(max_speed = 5.)

        # Check if there are east/north-oriented vectors
        for var in self.variables:
            for xvar, eastnorthvar in self.xy2eastnorth_mapping.items():
                if xvar in self.variables:
                    continue  # We have both x/y and east/north components
                if var in eastnorthvar:
                    logger.info('Variable %s will be rotated from %s' % (xvar, var))
                    self.variables.append(xvar)
                    if not hasattr(self, 'rotate_mapping'):
                        self.rotate_mapping = {}
                    self.rotate_mapping[xvar] = var

        # Adding variables which may be derived from existing ones
        for m in self.environment_mappings:
            em = self.environment_mappings[m]
            if em['active'] is False:
                logger.debug('Variable mapping: %s -> %s is not activated' % (em['input'], em['output']))
                continue
            self.activate_environment_mapping(m)

    def y_is_north(self):
        if (self.proj.crs.is_geographic and 'ob_tran' not in self.proj4) or '+proj=merc' in self.proj.srs:
            return True
        else:
            return False

    def index_of_closest_z(self, requested_z):
        """Return (internal) index of z closest to requested z.

        Thickness of layers (of ocean model) are not assumed to be constant.
        """
        ind_z = [np.abs(np.subtract.outer(
            self.z, requested_z)).argmin(0)]
        return ind_z, self.z[ind_z]

    def indices_min_max_z(self, z):
        """
        Return min and max indices of internal vertical dimension,
        covering the requested vertical positions.
        Needed when block is requested (True).

        Arguments:
            z: ndarray of floats, in meters
        """
        minIndex = (self.z <= z.min()).argmin() - 1
        maxIndex = (self.z >= z.max()).argmax()
        return minIndex, maxIndex

    def performance(self):
        '''Report the time spent on various tasks'''
        outStr = ''
        if hasattr(self, 'timing'):
            for cat, time in self.timing.items():
                time = str(time)[0:str(time).find('.') + 2]
                outStr += '%10s  %s\n' % (time, cat)
        return outStr

    def clip_boundary_pixels(self, numpix):
        '''Trim some (potentially bad) pixels along boundary'''
        logger.info('Trimming %i pixels from boundary' % numpix)
        self.xmin = self.xmin+numpix*self.delta_x
        self.xmax = self.xmax-numpix*self.delta_x
        self.ymin = self.ymin+numpix*self.delta_y
        self.ymax = self.ymax-numpix*self.delta_y
        self.shape = tuple([self.shape[0]-2*numpix,
                            self.shape[1]-2*numpix])
        self.clipped = numpix

    def plot(self, variable=None, vmin=None, vmax=None, time=None,
             filename=None, title=None, buffer=1, lscale='auto',
             cmap = None, cbar_label = None):
        """Plot geographical coverage of reader."""

        fig = plt.figure()

        if self.global_coverage():
            lonmin=-180
            lonmax=180
            latmin=-89
            latmax=89
        else:
            corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                     [self.ymax, self.ymin, self.ymax, self.ymin])
            lonmin = np.min(corners[0]) - buffer*2
            lonmax = np.max(corners[0]) + buffer*2
            latmin = np.min(corners[1]) - buffer
            latmax = np.max(corners[1]) + buffer

        # Initialise map
        if not self.global_coverage():
            # Stereographic projection centred on domain, if small domain
            x0 = (self.xmin + self.xmax) / 2
            y0 = (self.ymin + self.ymax) / 2
            lon0, lat0 = self.xy2lonlat(x0, y0)
            lon0 = lon0[0]
            lat0 = lat0[0]
            sp = ccrs.Stereographic(central_longitude=lon0, central_latitude=lat0)
            latmax = np.maximum(latmax, lat0)
            latmin = np.minimum(latmin, lat0)
        else:
            # Global map if reader domain is large
            sp = ccrs.Mercator()

        ax = fig.add_subplot(1, 1, 1, projection=sp)

        if lscale == 'auto':  # Custom lscale - this should be generalized to Basemodel also
            s = cfeature.AdaptiveScaler('coarse',
                (('low', 100), ('intermediate', 20), ('high', 5), ('full', 1)))
            lscale = s.scale_from_extent([lonmin, lonmax, latmin, latmax])

        # GSHHS coastlines
        f = cfeature.GSHHSFeature(scale=lscale, levels=[1])
        f._geometries_cache = {}
        ax.add_geometries(
            #f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
            f.geometries(),
            ccrs.PlateCarree(),
            facecolor=cfeature.COLORS['land'],
            edgecolor='black')

        gl = ax.gridlines(ccrs.PlateCarree())
        gl.top_labels = False

        # Get boundary
        if self.global_coverage():
            lon = np.array([-180, 180, 180, -180, -180])
            lat = np.array([-89, -89, 89, 89, -89])
        else:
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
            if self.global_coverage():
                ax.set_global()
            else:
                buf = (xsp.max()-xsp.min())*.1  # Some whitespace around polygon
                ax.set_extent([xsp.min()-buf, xsp.max()+buf, ysp.min()-buf, ysp.max()+buf], crs=sp)
        if title is None:
            plt.title(self.name)
        else:
            plt.title(title)
        plt.xlabel('Time coverage: %s to %s' %
                   (self.start_time, self.end_time))

        if time is None:
            time = self.start_time
        if variable is not None:
            rx = np.array([self.xmin, self.xmax])
            ry = np.array([self.ymin, self.ymax])
            if variable in self.derived_variables:
                data = self.get_variables(self.derived_variables[variable], time, rx, ry)
                self.__calculate_derived_environment_variables__(data)
            else:
                data = self.get_variables(variable, time, rx, ry)
            rx, ry = np.meshgrid(data['x'], data['y'])
            rx = np.float32(rx)
            ry = np.float32(ry)
            rlon, rlat = self.xy2lonlat(rx, ry)
            data[variable] = np.ma.masked_invalid(data[variable])

            if data[variable].ndim > 2:
                logger.warning('Ensemble data, plotting only first member')
                data[variable] = data[variable][0,:,:]
            if self.global_coverage():
                mappable = ax.pcolormesh(rlon, rlat, data[variable], vmin=vmin, vmax=vmax,
                                         transform=ccrs.PlateCarree())
            else:
                p = sp.transform_points(ccrs.PlateCarree(), rlon, rlat)
                mapx = p[:,:,0]
                mapy = p[:,:,1]
                mappable = ax.pcolormesh(mapx, mapy, data[variable], vmin=vmin, vmax=vmax,  cmap = cmap)

            cbar = fig.colorbar(mappable, orientation='horizontal', pad=.05, aspect=30, shrink=.4)
            if cbar_label:
                cbar.set_label(cbar_label)
            else:
                cbar.set_label(variable)

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

    def get_timeseries_at_position(self, lon, lat, variables=None,
                                   start_time=None, end_time=None, times=None):
        """ Get timeseries of variables from this reader at given position.
        """

        if times is None:
            if start_time is None:
                start_time = self.start_time
            if end_time is None:
                end_time = self.end_time
            times = [t for t in self.times if t >= start_time and t<= end_time]

        if variables is None:
            variables = self.variables

        # if len(self.covers_positions(lon=lon, lat=lat)[0]) == 0:
        #     return None

        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        if len(lon) == 1:
            lon = lon[0]*np.ones(len(times))
            lat = lat[0]*np.ones(len(times))

        data = {'time': times}
        for var in variables:
            data[var] = np.zeros(len(times))

        for i, time in enumerate(times):
            closest_time = min(self.times, key=lambda d: abs(d - time))
            d = self.get_variables_interpolated(
                lon=np.atleast_1d(lon[i]), lat=np.atleast_1d(lat[i]), z=np.atleast_1d(0),
                time=closest_time, variables=variables, rotate_to_proj='+proj=latlong')[0]
            for var in variables:
                data[var][i] = d[var][0]

        return data

    def shift_start_time(self, start_time):
        """Shift the time coverage of reader to match given start_time"""
        if self.start_time is None:
            logger.warning('Reader has no time coverage, cannot shift times')
            return
        shift_time = self.start_time - start_time
        self.start_time = start_time
        self.end_time = self.end_time - shift_time
        self.times = [t - shift_time for t in self.times]
