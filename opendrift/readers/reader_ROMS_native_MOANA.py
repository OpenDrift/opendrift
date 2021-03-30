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

########################################################################
# Specific reader_ROMS_native for Moana data
# https://www.moanaproject.org/
# http://thredds.moanaproject.org:8080/thredds/catalog/moana/ocean/NZB/v1.9/raw_3D/catalog.html
# 
# same as reader_ROMS_native.py but handles files with longitude >180, and has some sligtly different variables
# name mapping
# 
# updated 27/11/2020 to reflect changes in "official reader_ROMS_native.py"
# updated 15/1/2021 following significant changes on the basereader; this moana-specific reader now inherits from reader_ROMS_native
# 
# 
# Developed by S.Weppe MetOcean Solutions / MetServices NZ
########################################################################

from bisect import bisect_left, bisect_right
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

import numpy as np
from netCDF4 import num2date
import xarray as xr

from opendrift.readers.basereader import BaseReader, vector_pairs_xy, StructuredReader
from opendrift.readers import reader_ROMS_native

from opendrift.readers.roppy import depth


class Reader(reader_ROMS_native.Reader):
    def __init__(self, filename=None, name=None, gridfile=None):

        # Run constructor of parent Reader class ( reader_ROMS_native.Reader )
        super(Reader, self).__init__(filename=filename, name=name, gridfile=gridfile)

        # Add attribute to flag use of longitude>180 deg
        ##########################################################################################################
        # Modification : 
        # Opendrift internally uses longitudes -180<longitude<180, which may cause issues when
        # data is provided for 0<longitude<360 on netcdf files. This is the case for the ROMS Moana data
        # Here we add a flag to specify if that is indeed the case
        # search 'has_lon_0_360' for next adjustments
        ##########################################################################################################
        if  (self.lon[:]>=0).all() and (self.lon[:]<=360).all() : 
        # reader longitude using convention :  0<lon<360, or has longitude between 0 and +180 only
            self.has_lon_0_360 = True
        else:
            self.has_lon_0_360 = False
        ##########################################################################################################

    def covers_positions_xy(self, x, y, z=0):
        """
        Return indices of input points covered by reader.

        Arguments in native projection of reader.
        """
        z = np.atleast_1d(z)

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax)
                               & (z >= self.zmin) & (z <= self.zmax))[0]
        else:
            # there used a to be a if here when self.has_lon_0_360 = True.
            # may not be required anymore since it is in native projection
            indices = np.where((x >= self.xmin) & (x <= self.xmax)
                               & (y >= self.ymin) & (y <= self.ymax)
                               & (z >= self.zmin) & (z <= self.zmax))[0]

        try:
            return indices, x[indices], y[indices]
        except Exception as ex:
            logger.exception(ex)
            return indices, x, y

    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader."""

        # overloads method in variables.py
                
        if not self.has_lon_0_360 : # simple case - no overlap of 180W, use original code
            x, y = self.lonlat2xy(lon, lat)
            return self.covers_positions_xy(x, y, z)

        elif self.has_lon_0_360: #reader longitude using convention :  0<lon<360, or has longitude between 0 and +180 only
            lon360 = self.to_longitude_0_360(lon) # convert opendrift longitude from ([-180,180]) convention to [0,360] convention
            x, y = self.lonlat2xy(lon360, lat) # use the lon0_360 to get the x,y coordinates that will be input to covers_positions_xy()
            return self.covers_positions_xy(x, y, z) # use the 0<lon<360

    
    # def get_variables( ...)
    # we can remove get_variables() as it will be re-used from reader_ROMS_native.Reader

    def get_variables_interpolated(self,
                                   variables,
                                   profiles=None,
                                   profiles_depth=None,
                                   time=None,
                                   lon=None,
                                   lat=None,
                                   z=None,
                                   rotate_to_proj=None):
        """
        `get_variables_interpolated` is the main interface to
        :class:`opendrift.basemodel.OpenDriftSimulation`, and is responsible
        for returning variables at the correct positions.

        Readers should implement :meth:`_get_variables_interpolated_`.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.

            profiles: List of variable names that should be returned for the range in `profiles_depth`.

            profiles_depth: A range [z-start, z-end] for which to return values for profile-variables. The exact z-depth are given by the reader and returned as `z` variable in `env_profiles`.

            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).

            lon: N/A

            lat: N/A

            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)

            block: bool, see return below

            rotate_to_proj: N/A

          Returns:

            (env, env_profiles)

            Interpolated variables at x, y and z. `env` contains values at a fixed depth (`z`), while `env_profiles` contains depth-profiles in the range `profile_depth` for the variables listed in `profiles` for each element (in `x`, `y`). The exact depth is determined by the reader and specified in
            `env_profiles['z']`. Thus variables in `env_profiles` are not interpolated in z-direction.

        .. seealso::

            :meth:`get_variables_interpolated_xy`.

        """

        # overloads method with same name in variables.py

        # x, y = self.lonlat2xy(lon, lat)
        
        if not self.has_lon_0_360 : # simple case - no overlap of 180W, use original code
            x, y = self.lonlat2xy(lon, lat)
        elif self.has_lon_0_360 :
            lon360 = self.to_longitude_0_360(lon) # convert opendrift longitude from ([-180,180]) convention to [0,360] convention
            x, y = self.lonlat2xy(lon360, lat) # use the lon0_360 to get the correct x,y coordinates that will be input to get_variables_interpolated_xy()

        env, env_profiles = self.get_variables_interpolated_xy(
            variables, profiles, profiles_depth, time, x, y, z, rotate_to_proj)

        return env, env_profiles

    # def get_variables(self, requested_variables, time=None,
    #                   x=None, y=None, z=None):
    #     '''
    #     Just re-uses the one from reader_ROMS_native now
    #     '''

    def to_longitude_0_360(self,lon180):
        '''
        converts  longitude to different convention
        -180<longitude<180 toward 0<longitude<360

        '''

        lon360 = lon180.copy() # convert lon to [0-360] convention
        lon360[np.where(lon360<0)] = 360+ lon360[np.where(lon360<0)]
        return lon360

def rotate_vectors_angle(u, v, radians):
    u2 = u*np.cos(radians) - v*np.sin(radians)
    v2 = u*np.sin(radians) + v*np.cos(radians)
    return u2, v2
