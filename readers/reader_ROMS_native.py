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
from datetime import datetime, timedelta
from bisect import bisect_left, bisect_right

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date

from reader import Reader
from roppy import depth


class Reader(Reader):

    # Map ROMS variable names to CF standard_name
    ROMS_variable_mapping = {
        'mask_psi': 'land_binary_mask',
        'h': 'sea_floor_depth',
        'zeta': 'sea_surface_height',
        'u': 'x_sea_water_velocity',
        'v': 'y_sea_water_velocity',
        'temp': 'sea_water_temperature',
        'salt': 'sea_water_salinity',
        'uice': 'sea_ice_x_velocity',
        'vice': 'sea_ice_y_velocity',
        'aice': 'sea_ice_area_fraction',
        'hice': 'sea_ice_thickness',
        'Uwind': 'x_wind',
        'Vwind': 'y_wind'}

    zbuffer = 1  # Vertical buffer of block around elements

    # z-levels to which sigma-layers may be interpolated
    zlevels = [0, -.5, -1, -3, -5, -10, -25, -50, -75, -100, -150, -200,
               -250, -300, -400, -500, -600, -700, -800, -900, -1000, -1500,
               -2000, -2500, -3000, -3500, -4000, -4500, -5000, -5500, -6000,
               -6500, -7000, -7500, -8000]

    def __init__(self, filename=None, name=None, gridfile=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        if name is None:
            self.name = filename
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            logging.info('Opening dataset: ' + filename)
            self.Dataset = MFDataset(filename, 'r')
        except:
            try:
                logging.info('Could not open with MFDataset, '
                             'trying with Dataset:')
                self.Dataset = Dataset(filename, 'r')
            except:
                raise ValueError('Could not open ' + filename +
                                 ' with netCDF4 library')

        if 's_rho' not in self.Dataset.variables:
            dimensions = 2
        else:
            dimensions = 3

        if dimensions == 3:
            # Read sigma-coordinate values
            try:
                self.sigma = self.Dataset.variables['s_rho'][:]
            except:
                num_sigma = len(self.Dataset.dimensions['s_rho'])
                logging.warning('s_rho not available in dataset, constructing from'
                                ' number of layers (%s).' % num_sigma)
                self.sigma = (np.arange(num_sigma)+.5-num_sigma)/num_sigma

            # Read sigma-coordinate transform parameters
            self.Cs_r = self.Dataset.variables['Cs_r'][:]
            self.hc = self.Dataset.variables['hc'][:]

            self.num_layers = len(self.sigma)
        else:
            self.num_layers = 1
            self.ROMS_variable_mapping['ubar'] = 'x_sea_water_velocity'
            self.ROMS_variable_mapping['vbar'] = 'y_sea_water_velocity'
            del self.ROMS_variable_mapping['u']
            del self.ROMS_variable_mapping['v']

        if 'lat_rho' in self.Dataset.variables:
            # Horizontal oordinates and directions
            self.lat = self.Dataset.variables['lat_rho'][:]
            self.lon = self.Dataset.variables['lon_rho'][:]
        else:
            if gridfile is None:
                raise ValueError(filename + ' does not contain lon/lat '
                                 'arrays, please supply a grid-file '
                                 '"gridfile=<grid_file>"')
            else:
                gf = Dataset(gridfile)
                self.lat = gf.variables['lat_rho'][:]
                self.lon = gf.variables['lon_rho'][:]

        # Get time coverage
        ocean_time = self.Dataset.variables['ocean_time']
        time_units = ocean_time.__dict__['units']
        self.times = num2date(ocean_time[:], time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
        else:
            self.time_step = None

        # x and y are rows and columns for unprojected datasets
        self.xmin = 0.
        self.xmax = np.float(len(self.Dataset.dimensions['xi_rho'])) - 1
        self.delta_x = 1.
        self.ymin = 0.
        self.ymax = np.float(len(self.Dataset.dimensions['eta_rho'])) - 1
        self.delta_y = 1.

        self.name = 'roms native'

        # Find all variables having standard_name
        self.variables = []
        for var_name in self.Dataset.variables:
            if var_name in self.ROMS_variable_mapping.keys():
                var = self.Dataset.variables[var_name]
                self.variables.append(self.ROMS_variable_mapping[var_name])

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        variables = {}

        # Find horizontal indices corresponding to requested x and y
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int)
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int)
        indx_el = indx
        indy_el = indy
        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            # Avoiding the last pixel in each dimension, since there are
            # several grids which are shifted (rho, u, v, psi)
            indx = np.arange(np.max([0, indx.min()-buffer]),
                             np.min([indx.max()+buffer, self.lon.shape[1]-1]))
            indy = np.arange(np.max([0, indy.min()-buffer]),
                             np.min([indy.max()+buffer, self.lon.shape[0]-1]))
        else:
            indx[outside[0]] = 0  # To be masked later
            indy[outside[0]] = 0

        # Find depth levels covering all elements
        if z.min() == 0:
            indz = self.num_layers - 1  # surface layer
            variables['z'] = 0

        else:
            # Find the range of indices covering given z-values
            if not hasattr(self, 'sea_floor_depth'):
                logging.debug('Reading sea floor depth...')
                self.sea_floor_depth = self.Dataset.variables['h'][:]
            indxgrid, indygrid = np.meshgrid(indx, indy)
            H = self.sea_floor_depth[indygrid, indxgrid]
            z_rho = depth.sdepth(H, self.hc, self.Cs_r)
            # Element indices must be relative to extracted subset
            indx_el = indx_el - indx.min()
            indy_el = indy_el - indy.min()

            # Loop to find the layers covering the requested z-values
            indz_min = 0
            for i in range(self.num_layers):
                if np.min(z-z_rho[i, indy_el, indx_el]) > 0:
                    indz_min = i
                if np.max(z-z_rho[i, indy_el, indx_el]) > 0:
                    indz_max = i
            indz = range(np.maximum(0, indz_min-self.zbuffer),
                         np.minimum(self.num_layers,
                                     indz_max + 1 + self.zbuffer))
            z_rho = z_rho[indz, :, :]
            # Determine the z-levels to which to interpolate
            zi1 = np.maximum(0, bisect_left(-np.array(self.zlevels),
                                            -z.max()) - self.zbuffer)
            zi2 = np.minimum(len(self.zlevels),
                             bisect_right(-np.array(self.zlevels),
                                          -z.min()) + self.zbuffer)
            variables['z'] = np.array(self.zlevels[zi1:zi2])

        for par in requested_variables:
            varname = [name for name, cf in
                       self.ROMS_variable_mapping.items() if cf == par]
            var = self.Dataset.variables[varname[0]]
            var.set_auto_maskandscale(True)

            if var.ndim == 2:
                variables[par] = var[indy, indx]
            elif var.ndim == 3:
                variables[par] = var[indxTime, indy, indx]
            elif var.ndim == 4:
                variables[par] = var[indxTime, indz, indy, indx]
                # Regrid from sigma to z levels
                if len(np.atleast_1d(indz)) > 1:
                    logging.debug('sigma to z for ' + varname[0])
                    try:  # Converting masked values to nan, since
                          # not handled by multi_zslice
                        variables[par][variables[par].mask] = np.nan
                    except:
                        pass
                    variables[par] = depth.multi_zslice(variables[par],
                                                    z_rho, variables['z'])
                    # Re-adding mask removed by multi_zslice:
                    variables[par] = np.ma.masked_invalid(variables[par])
            else:
                raise Exception('Wrong dimension of variable: '
                                + self.variable_mapping[par])

            # If 2D array is returned due to the fancy slicing methods
            # of netcdf-python, we need to take the diagonal
            if variables[par].ndim > 1 and block is False:
                variables[par] = variables[par].diagonal()

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par], ndmin=2, mask=False)
            if block is False:
                variables[par].mask[outside[0]] = True

        if 'land_binary_mask' in variables.keys():
            variables['land_binary_mask'] = 1 - variables['land_binary_mask']

        if block is True:
            variables['x'] = indx
            variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['time'] = nearestTime

        if 'x_sea_water_velocity' or 'sea_ice_x_velocity' in variables.keys():
            # We must rotate current vectors
            if not hasattr(self, 'angle_xi_east'):
                logging.debug('Reading angle between xi and east...')
                self.angle_xi_east = self.Dataset.variables['angle'][:]
            rad = self.angle_xi_east[np.meshgrid(indy, indx)].T
            if 'x_sea_water_velocity' in variables.keys():
                variables['x_sea_water_velocity'], \
                    variables['y_sea_water_velocity'] = rotate_vectors(
                        variables['x_sea_water_velocity'],
                        variables['y_sea_water_velocity'], rad)
            if 'sea_ice_x_velocity' in variables.keys():
                variables['sea_ice_x_velocity'], \
                    variables['sea_ice_y_velocity'] = rotate_vectors(
                        variables['sea_ice_x_velocity'],
                        variables['sea_ice_y_velocity'], rad)

        return variables


def rotate_vectors(u, v, radians):
    u = u*np.cos(radians) - v*np.sin(radians)
    v = u*np.sin(radians) + v*np.cos(radians)
    return u, v
