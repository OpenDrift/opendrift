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
from netCDF4 import Dataset, MFDataset, num2date

from opendrift.readers.basereader import BaseReader

try:
    import xarray as xr
    has_xarray = True
except:
    has_xarray = False

def proj_from_CF_dict(c):

    # This method should be extended to other projections:
    # https://cfconventions.org/wkt-proj-4.html
    if not 'grid_mapping_name' in c:
        raise ValueError('grid_mapping not given in dictionary')
    gm = c['grid_mapping_name']

    if gm == 'polar_stereographic':
        if 'earth_radius' in c:
            earth_radius = c['earth_radius']
        else:
            earth_radius = 6371000.
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

    proj = pyproj.Proj(proj4)

    return proj4, proj


class Reader(BaseReader):
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

    def __init__(self, filename=None, name=None, proj4=None, standard_name_mapping={}):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            self.logger.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                self.logger.info('Opening files with MFDataset')
                if has_xarray:
                    self.Dataset = xr.open_mfdataset(filename, concat_dim='time', combine='nested',
                                                     decode_times=False)
                else:
                    try:
                        self.Dataset = MFDataset(filename, aggdim='time')
                    except Exception as e:
                        self.logger.warning('Aggdim *time* does not exist')
                        self.Dataset = MFDataset(filename)
            else:
                self.logger.info('Opening file with Dataset')
                if has_xarray:
                    self.Dataset = xr.open_dataset(filename, decode_times=False)
                else:
                    self.Dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

        self.logger.debug('Finding coordinate variables.')
        if proj4 is not None:  # If user has provided a projection apriori
            self.proj4 = proj4
        # Find x, y and z coordinates
        for var_name in self.Dataset.variables:
            self.logger.debug('Parsing variable: ' +  var_name)
            var = self.Dataset.variables[var_name]
            #if var.ndim > 1:
            #    continue  # Coordinates must be 1D-array
            if has_xarray:
                attributes = var.attrs
                att_dict = var.attrs
            else:
                attributes = var.ncattrs()
                att_dict = var.__dict__
            standard_name = ''
            long_name = ''
            axis = ''
            units = ''
            CoordinateAxisType = ''
            if not hasattr(self, 'proj4'):
                for att in attributes:
                    if 'proj4' in att:
                        if has_xarray:
                            self.proj4 = str(att_dict[att])
                        else:
                            self.proj4 = str(var.__getattr__(att))
                    else:
                        if 'grid_mapping_name' in att:
                            mapping_dict = att_dict
                            self.logger.debug(
                                ('Parsing CF grid mapping dictionary:'
                                ' ' + str(mapping_dict)))
                            try:
                                self.proj4, proj =\
                                    proj_from_CF_dict(mapping_dict)
                            except:
                                self.logger.info('Could not parse CF grid_mapping')

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
            # has_xarray checks in each case below to avoid loading
            # data if it isn't a coord
            # is there a better way??
            if standard_name == 'longitude' or \
                    CoordinateAxisType == 'Lon' or \
                    long_name.lower() == 'longitude':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                self.lon = var_data
                lon_var_name = var_name
            if standard_name == 'latitude' or \
                    CoordinateAxisType == 'Lat' or \
                    long_name.lower() == 'latitude':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                self.lat = var_data
                lat_var_name = var_name
            if axis == 'X' or \
                    standard_name == 'projection_x_coordinate':
                self.xname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                elif units == '100  km':
                    unitfactor = 100000
                else:
                    unitfactor = 1
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                x = var_data*unitfactor
                self.numx = var_data.shape[0]
            if axis == 'Y' or \
                    standard_name == 'projection_y_coordinate':
                self.yname = var_name
                # Fix for units; should ideally use udunits package
                if units == 'km':
                    unitfactor = 1000
                elif units == '100  km':
                    unitfactor = 100000
                else:
                    unitfactor = 1
                self.unitfactor = unitfactor
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                y = var_data*unitfactor
                self.numy = var_data.shape[0]
            if standard_name == 'depth' or axis == 'Z':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                if var_data.ndim == 1:
                    if 'positive' not in attributes or \
                            att_dict['positive'] == 'up':
                        self.z = var_data
                    else:
                        self.z = -var_data
            if standard_name == 'time' or axis == 'T' or var_name in ['time', 'vtime']:
                # Read and store time coverage (of this particular file)
                if has_xarray:
                    var_data = var.values
                else:
                    var.set_auto_mask(False)  # workaround
                    var_data = var[:]
                time = var_data
                time_units = units

                if isinstance(time[0], np.bytes_):
                    # This hack is probably only necessary for CERSAT/GELOBCURRENT
                    time = [t.decode('ascii') for t in time]
                    self.times = [datetime.fromisoformat(t.replace('Z', '')) for t in time]
                elif time.ndim == 2:
                    self.times = [datetime.fromisoformat(''.join(t).replace('Z', '')) for t in time.astype(str)]
                else:
                    self.times = num2date(time, time_units)
                #if has_xarray:
                #    self.times = [datetime.utcfromtimestamp((OT -
                #        np.datetime64('1970-01-01T00:00:00Z')
                #            ) / np.timedelta64(1, 's')) for OT in time]
                #else:
                #    self.times = num2date(time, time_units)
                self.start_time = self.times[0]
                self.end_time = self.times[-1]
                if len(self.times) > 1:
                    self.time_step = self.times[1] - self.times[0]
                else:
                    self.time_step = None
            if standard_name == 'realization':
                if has_xarray:
                    var_data = var.values
                else:
                    var_data = var[:]
                self.realizations = var_data
                self.logger.debug('%i ensemble members available'
                              % len(self.realizations))

        if 'x' not in locals():
            if self.lon.ndim == 1:
                x = self.lon[:]
                self.xname = lon_var_name
                self.numx = len(x)
            #else:
            #    raise ValueError('Did not find x-coordinate variable')
        if 'y' not in locals():
            if self.lat.ndim == 1:
                y = self.lat[:]
                self.yname = lat_var_name
                self.numy = len(y)
            #else:
            #    raise ValueError('Did not find y-coordinate variable')

        if not hasattr(self, 'unitfactor'):
            self.unitfactor = 1
        if 'x' in locals() and 'y' in locals():
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
        else:
            if hasattr(self, 'lon') and hasattr(self, 'lat'):
                self.logger.info('No projection found, using lon/lat arrays')
                self.xname = lon_var_name
                self.yname = lat_var_name
            else:
                raise ValueError('Neither x/y-coordinates or lon/lat arrays found')

        if not hasattr(self, 'proj4'):
            if self.lon.ndim == 1:
                self.logger.debug('Lon and lat are 1D arrays, assuming latong projection')
                self.proj4 = '+proj=latlong'
            elif self.lon.ndim == 2:
                self.logger.debug('Reading lon lat 2D arrays, since projection is not given')
                self.lon = self.lon[:]
                self.lat = self.lat[:]
                self.projected = False

        if hasattr(self, 'proj4') and 'latlong' in self.proj4 and hasattr(self, 'xmax') and self.xmax > 360:
            self.logger.info('Longitudes > 360 degrees, subtracting 360')
            self.xmin -= 360
            self.xmax -= 360
            self.x -= 360
            self.x -= 360

        # Find all variables having standard_name
        self.variable_mapping = {}
        for var_name in self.Dataset.variables:
            if var_name in [self.xname, self.yname, 'depth']:
                continue  # Skip coordinate variables
            var = self.Dataset.variables[var_name]
            if has_xarray:
                attributes = var.attrs
                att_dict = var.attrs
            else:
                attributes = var.ncattrs()
                att_dict = var.__dict__
            if 'standard_name' in attributes:
                standard_name = str(att_dict['standard_name'])
                if standard_name in self.variable_aliases:  # Mapping if needed
                    standard_name = self.variable_aliases[standard_name]
                self.variable_mapping[standard_name] = str(var_name)
            elif var_name in standard_name_mapping:
                # User may specify mapping if standard_name is missing
                standard_name = standard_name_mapping[var_name]
                self.variable_mapping[standard_name] = str(var_name)

        self.variables = list(self.variable_mapping.keys())

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
        if hasattr(self, 'clipped'):
            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped
        # If x or y coordinates are decreasing, we need to flip
        if hasattr(self, 'x'):
            if self.x[0] > self.x[-1]:
                indx = len(self.x) - indx
            if self.y[0] > self.y[-1]:
                indy = len(self.y) - indy
        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            if self.global_coverage():
                #indx = np.arange(indx.min()-buffer, indx.max()+buffer)
                if indx.min() < 0:  # Primitive fix for crossing 0-meridian
                    indx = indx + self.numx
                indx = np.arange(np.max([0, indx.min()-buffer]),
                                 np.min([indx.max()+buffer, self.numx]))
            else:
                indx = np.arange(np.max([0, indx.min()-buffer]),
                                 np.min([indx.max()+buffer, self.numx]))
            indy = np.arange(np.max([0, indy.min()-buffer]),
                             np.min([indy.max()+buffer, self.numy]))
        else:
            indx[outside] = 0  # To be masked later
            indy[outside] = 0
        if indx.min() <= 0 and indx.max() >= self.numx:
            indx = np.arange(0, self.numx)

        variables = {}

        if indx.min() < 0 and indx.max() > 0:
            self.logger.debug('Requested data block is not continous in file'+
                          ', must read two blocks and concatenate.')
            indx_left = indx[indx<0] + self.numx  # Shift to positive indices
            indx_right = indx[indx>=0]
            continous = False
        else:
            continous = True
        for par in requested_variables:
            if hasattr(self, 'rotate_mapping') and par in self.rotate_mapping:
                self.logger.debug('Using %s to retrieve %s' %
                    (self.rotate_mapping[par], par))
                if par not in self.variable_mapping:
                    self.variable_mapping[par] = \
                        self.variable_mapping[
                            self.rotate_mapping[par]]
            var = self.Dataset.variables[self.variable_mapping[par]]

            ensemble_dim = None
            if continous is True:
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

            if has_xarray is True:
                variables[par] = np.asarray(variables[par])

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
            with np.errstate(invalid='ignore'):
                variables[par] = np.ma.masked_outside(
                    variables[par], -30000, 30000)

            # Ensemble blocks are split into lists
            if ensemble_dim is not None:
                num_ensembles = variables[par].shape[ensemble_dim]
                self.logger.debug('Num ensembles: %i ' % num_ensembles)
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
            if self.projected is True:
                variables['x'] = \
                    self.Dataset.variables[self.xname][indx]*self.unitfactor
                variables['y'] = \
                    self.Dataset.variables[self.yname][indy]*self.unitfactor
            else:
                variables['x'] = indx
                variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y
        if has_xarray is True:
            variables['x'] = np.asarray(variables['x'])
            variables['y'] = np.asarray(variables['y'])
        if self.global_coverage():
            # TODO: this should be checked
            if self.xmax + self.delta_x >= 360 and variables['x'].max() > 180:
                variables['x'] -= 360
                    
        variables['time'] = nearestTime

        # Rotate any east/north vectors if necessary
        if hasattr(self, 'rotate_mapping'):
            if self.y_is_north() is True:
                self.logger.debug('North is up, no rotation necessary')
            else:
                self.rotate_variable_dict(variables)

        return variables
