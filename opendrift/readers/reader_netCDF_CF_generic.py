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

from datetime import datetime
import pyproj
import numpy as np
from netCDF4 import num2date
import logging
logger = logging.getLogger(__name__)

from opendrift.readers.basereader import BaseReader, StructuredReader
import xarray as xr

def proj_from_CF_dict(c):

    # This method should be extended to other projections:
    # https://cfconventions.org/wkt-proj-4.html
    if not 'grid_mapping_name' in c:
        raise ValueError('grid_mapping not given in dictionary')
    gm = c['grid_mapping_name']

    if 'earth_radius' in c:
        earth_radius = c['earth_radius']
    else:
        earth_radius = 6371000.
    if gm == 'polar_stereographic':
        lon_0 = 0  # default, but dangerous
        for l0 in ['longitude_of_projection_origin',
                   'longitude_of_central_meridian',
                   'straight_vertical_longitude_from_pole']:
            if l0 in c:
                lon_0 = c[l0]
        if 'latitude_of_origin' in c:
            lat_ts = c['latitude_of_origin']
        else:
            lat_ts = c['latitude_of_projection_origin']
        if 'false_easting' in c:
            x0 = c['false_easting']
        else:
            x0 = 0  # dangerous?
        if 'false_northing' in c:
            y0 = c['false_northing']
        else:
            y0 = 0  # dangerous: is there a better default?
        if 'scale_factor_at_projection_origin' in c:
            k0 = c['scale_factor_at_projection_origin']
        else:
            k0 = 1.0
        proj4 = ('+proj={!s} +lat_0={!s} +lon_0={!s} +lat_ts={!s} '
                 '+k_0={!s} +x_0={!s} +y_0={!s} +units=m +a={!s} '
                 '+no_defs'.format('stere',
                                   c['latitude_of_projection_origin'],
                                   lon_0, lat_ts, k0, x0, y0, earth_radius)
                 )

    elif gm == 'rotated_latitude_longitude':
        proj4 = '+proj=ob_tran +o_proj=longlat +lon_0=%s +o_lat_p=%s +R=%s +no_defs' % (
                c['grid_north_pole_longitude']-180, c['grid_north_pole_latitude'], earth_radius)

    proj = pyproj.Proj(proj4)

    return proj4, proj


class Reader(StructuredReader, BaseReader):
    """
    A reader for `CF-compliant <https://cfconventions.org/>`_ netCDF files. It can take a single file, or a file pattern.

    Args:
        :param filename: A single netCDF file, or a pattern of files. The
                         netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, requiered.

        :param name: Name of reader
        :type name: string, optional

        :param proj4: PROJ.4 string describing projection of data.
        :type proj4: string, optional

    Example:

    .. code::

       from opendrift.readers.reader_netCDF_CF_generic import Reader
       r = Reader("arome_subset_16Nov2015.nc")

    Several files can be specified by using a pattern:

    .. code::

       from opendrift.readers.reader_netCDF_CF_generic import Reader
       r = Reader("*.nc")

    or an OPeNDAP URL can be used:

    .. code::

       from opendrift.readers.reader_netCDF_CF_generic import Reader
       r = Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')


    """

    def __init__(self, filename=None, name=None, proj4=None, standard_name_mapping={}, ensemble_member=None):
        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            logger.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                logger.info('Opening files with MFDataset')
                self.Dataset = xr.open_mfdataset(filename, data_vars='minimal', coords='minimal',
                                                 chunks={'time': 1}, decode_times=False)
            elif ensemble_member is not None:
                self.Dataset = xr.open_dataset(filename, decode_times=False).isel(ensemble_member=ensemble_member)
            else:
                self.Dataset = xr.open_dataset(filename, decode_times=False)
        except Exception as e:
            raise ValueError(e)

        # NB: check below might not be waterproof
        if 'ocean_time' in self.Dataset.dims and 'eta_u' in self.Dataset.dims and \
                'eta_rho' in self.Dataset.dims:
            raise ValueError('This seems to be a ROMS native file, should use ROMS native reader instead')

        logger.debug('Finding coordinate variables.')
        if proj4 is not None:  # If user has provided a projection apriori
            self.proj4 = proj4
        # Find x, y and z coordinates
        lon_var_name = None
        lat_var_name = None
        self.unitfactor = 1
        self.realizations = None
        for var_name in self.Dataset.variables:
            var = self.Dataset.variables[var_name]

            if self.proj4 is None:
                if 'proj4' in var.attrs:
                    self.proj4 = str(var.attrs['proj4'])
                elif 'proj4_string' in var.attrs:
                    self.proj4 = str(var.attrs['proj4_string'])
                elif 'grid_mapping_name' in var.attrs:
                    logger.debug(
                        ('Parsing CF grid mapping dictionary:'
                        ' ' + str(var.attrs)))
                    try:
                        self.proj4, proj =\
                            proj_from_CF_dict(var.attrs)
                    except:
                        logger.info('Could not parse CF grid_mapping')

            standard_name = var.attrs['standard_name'] if 'standard_name' in var.attrs else ''
            long_name = var.attrs['long_name'] if 'long_name' in var.attrs else ''
            axis = var.attrs['axis'] if 'axis' in var.attrs else ''
            units = var.attrs['units'] if 'units' in var.attrs else ''
            CoordinateAxisType = var.attrs['CoordinateAxisType'] if 'CoordinateAxisType' in var.attrs else ''
            if standard_name == 'longitude' or \
                    CoordinateAxisType.lower() == 'lon' or \
                    long_name.lower() == 'longitude' or \
                    var_name.lower() in ['longitude', 'lon']:
                lon_var_name = var_name
            if standard_name == 'latitude' or \
                    CoordinateAxisType.lower() == 'lat' or \
                    long_name.lower() == 'latitude' or \
                    var_name.lower() in ['latitude', 'lat']:
                lat_var_name = var_name
            if (axis == 'X' or standard_name == 'projection_x_coordinate') \
                    and var.ndim == 1:
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    self.unitfactor = 1000
                elif units == '100  km':
                    self.unitfactor = 100000
                var_data = var.values
                x = var_data*self.unitfactor
            if (axis == 'Y' or standard_name == 'projection_y_coordinate') \
                    and var.ndim == 1:
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    self.unitfactor = 1000
                elif units == '100  km':
                    self.unitfactor = 100000
                var_data = var.values
                y = var_data*self.unitfactor
            if standard_name == 'depth' or axis == 'Z':
                var_data = var.values
                if var_data.ndim == 1:
                    if 'positive' not in var.attrs or \
                            var.attrs['positive'] == 'up':
                        self.z = var_data
                    else:
                        self.z = -var_data
            if standard_name == 'time' or axis == 'T' or var_name in ['time', 'vtime']:
                # Read and store time coverage (of this particular file)
                var_data = var.values
                time = var_data
                time_units = units

                if isinstance(time[0], np.bytes_):
                    # This hack is probably only necessary for CERSAT/GELOBCURRENT
                    time = [t.decode('ascii') for t in time]
                    self.times = [datetime.fromisoformat(t.replace('Z', '')) for t in time]
                elif time.ndim == 2:
                    self.times = [datetime.fromisoformat(''.join(t).replace('Z', '')) for t in time.astype(str)]
                else:
                    if 'calendar' in var.attrs:
                        calendar = var.attrs['calendar']
                    else:
                        calendar = 'standard'
                    self.times = num2date(time, time_units, calendar=calendar)
                self.start_time = self.times[0]
                self.end_time = self.times[-1]
                if len(self.times) > 1:
                    self.time_step = self.times[1] - self.times[0]
                else:
                    self.time_step = None
            if standard_name == 'realization':
                if ensemble_member == None:
                    var_data = var.values
                    self.realizations = var_data
                    logger.debug('%i ensemble members available'
                                % len(self.realizations))

        # Temporary workaround for Barents EPS model
        if self.realizations is None and 'ensemble_member' in self.Dataset.dims:
            self.realizations = np.arange(self.Dataset.dims['ensemble_member'])

        #########################
        # Summary of geolocation
        #########################
        if 'x' in locals():
            if x[1]-x[0] == 1 and y[1]-y[0] == 1:
                logger.info('deltaX and deltaY are 1, interpreting dataset as unprojected')
                projected = False
            else:
                projected = True
        if 'x' not in locals() or projected is False:  # No x/y-coordinates were detected
            if lon_var_name is None:
                raise ValueError('No geospatial coordinates were detected, cannot geolocate dataset')
            # We load lon and lat arrays into memory
            lon_var = self.Dataset.variables[lon_var_name]
            lat_var = self.Dataset.variables[lat_var_name]
            if lon_var.ndim == 1:
                logger.debug('Lon and lat are 1D arrays - using as projection coordinates')
                x = lon_var.data
                y = lat_var.data
                self.xname = lon_var_name
                self.yname = lat_var_name
                if self.proj4 is None:
                    self.proj4 = '+proj=latlong'
            elif lon_var.ndim == 2:
                logger.debug('Lon and lat are 2D arrays - dataset is unprojected')
                self.lon = lon_var.data
                self.lat = lat_var.data
                self.projected = False
            elif lon_var.ndim == 3:
                logger.debug('Lon lat are 3D arrays, reading first time')
                self.lon = lon_var[0,:,:].data
                self.lat = lat_var[0,:,:].data
                self.projected = False
        else:
            if self.proj4 is None:
                logger.info('Grid coordinates are detected, but proj4 string not given: assuming latlong')
                self.proj4 = '+proj=latlong'

        if 'x' in locals() and 'y' in locals():
            self.numx = len(x)
            self.numy = len(y)
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

        if self.proj4 is not None and 'latlong' in self.proj4 and self.xmax is not None and self.xmax > 360:
            logger.info('Longitudes > 360 degrees, subtracting 360')
            self.xmin -= 360
            self.xmax -= 360
            self.x -= 360

        ##########################################
        # Find all variables having standard_name
        ##########################################
        self.variable_mapping = {}
        skipvars = []
        for var_name in self.Dataset.variables:
            var = self.Dataset.variables[var_name]
            if var.ndim < 2:
                continue
            if var_name in standard_name_mapping:
                # User may specify mapping if standard_name is missing, or to override existing
                standard_name = standard_name_mapping[var_name]
                self.variable_mapping[standard_name] = str(var_name)
            elif 'standard_name' in var.attrs:
                standard_name = str(var.attrs['standard_name'])
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)
            else:
                skipvars.append(var_name)

        if len(skipvars) > 0:
            logger.debug('Skipped variables without standard_name: %s' % skipvars)

        self.variables = list(self.variable_mapping.keys())

        # Run constructor of parent Reader class
        super().__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None,
                      indrealization=None):

        requested_variables, time, x, y, z, _outside = self.check_arguments(
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
            if self.realizations is not None:
                indrealization = range(len(self.realizations))
            else:
                indrealization = None

        # Find indices corresponding to requested x and y
        if hasattr(self, 'clipped'):
            clipped = self.clipped
        else: clipped = 0

        if self.global_coverage():
            if self.lon_range() == '0to360':
                x = np.mod(x, 360)  # Shift x/lons to 0-360 
            elif self.lon_range() == '-180to180':
                x = np.mod(x + 180, 360) - 180 # Shift x/lons to -180-180
        indx = np.floor(np.abs(x-self.x[0])/self.delta_x-clipped).astype(int) + clipped
        indy = np.floor(np.abs(y-self.y[0])/self.delta_y-clipped).astype(int) + clipped
        buffer = self.buffer  # Adding buffer, to cover also future positions of elements
        indy = np.arange(np.max([0, indy.min()-buffer]),
                         np.min([indy.max()+buffer, self.numy]))
        indx = np.arange(indx.min()-buffer, indx.max()+buffer+1)

        if self.global_coverage() and indx.min() < 0 and indx.max() > 0 and indx.max() < self.numx:
            logger.debug('Requested data block is not continuous in file'+
                          ', must read two blocks and concatenate.')
            indx_left = indx[indx<0] + self.numx  # Shift to positive indices
            indx_right = indx[indx>=0]
            if indx_right.max() >= indx_left.min():  # Avoid overlap
                indx_right = np.arange(indx_right.min(), indx_left.min())
            continuous = False
        else:
            continuous = True
            indx = np.arange(np.max([0, indx.min()]),
                             np.min([indx.max(), self.numx]))

        variables = {}

        for par in requested_variables:
            if hasattr(self, 'rotate_mapping') and par in self.rotate_mapping:
                logger.debug('Using %s to retrieve %s' %
                    (self.rotate_mapping[par], par))
                if par not in self.variable_mapping:
                    self.variable_mapping[par] = \
                        self.variable_mapping[
                            self.rotate_mapping[par]]
            var = self.Dataset.variables[self.variable_mapping[par]]

            ensemble_dim = None
            if continuous is True:
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

            variables[par] = np.asarray(variables[par])

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par],
                                         ndmin=2, mask=False)
            # Mask extreme values which might have slipped through
            with np.errstate(invalid='ignore'):
                variables[par] = np.ma.masked_outside(
                    variables[par], -30000, 30000)

            # Ensemble blocks are split into lists
            if ensemble_dim is not None:
                num_ensembles = variables[par].shape[ensemble_dim]
                logger.debug('Num ensembles: %i ' % num_ensembles)
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
        if self.projected is True:
            variables['x'] = \
                self.Dataset.variables[self.xname][indx]*self.unitfactor
            variables['y'] = \
                self.Dataset.variables[self.yname][indy]*self.unitfactor
        else:
            variables['x'] = indx
            variables['y'] = indy
        variables['x'] = np.asarray(variables['x'], dtype=np.float32)
        variables['y'] = np.asarray(variables['y'], dtype=np.float32)

        variables['time'] = nearestTime

        # Rotate any east/north vectors if necessary
        if hasattr(self, 'rotate_mapping'):
            if self.y_is_north() is True:
                logger.debug('North is up, no rotation necessary')
            else:
                self.rotate_variable_dict(variables)

        if hasattr(self, 'shift_x'):
            # "hidden feature": if reader.shift_x and reader.shift_y are defined,
            # the returned fields are shifted this many meters in the x- and y directions
            # E.g. reader.shift_x=10000 gives a shift 10 km eastwards (if x is east direction)
            if self.proj.crs.is_geographic:  # meters to degrees
                shift_y = (self.shift_y/111000)
                shift_x = (self.shift_x/111000)*np.cos(np.radians(variables['y']))
                logger.info('Shifting x between %s and %s' % (shift_x.min(), shift_x.max()))
                logger.info('Shifting y with %s m' % shift_y)
            else:
                shift_x = self.shift_x
                shift_y = self.shift_y
                logger.info('Shifting x with %s m' % shift_x)
                logger.info('Shifting y with %s m' % shift_y)
            variables['x'] += shift_x
            variables['y'] += shift_y

        return variables
