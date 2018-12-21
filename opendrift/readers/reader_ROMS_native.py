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
from bisect import bisect_left, bisect_right
from datetime import datetime

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
try:
    import xarray as xr
    has_xarray = True
except:
    has_xarray = False
has_xarray = False  # Temporary disabled

from opendrift.readers.basereader import BaseReader, vector_pairs_xy
from opendrift.readers.roppy import depth


class Reader(BaseReader):

    def __init__(self, filename=None, name=None, gridfile=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        # Map ROMS variable names to CF standard_name
        self.ROMS_variable_mapping = {
            # Removing (temoprarily) land_binary_mask from ROMS-variables,
            # as this leads to trouble with linearNDFast interpolation
            'mask_rho': 'land_binary_mask',
            'mask_psi': 'land_binary_mask',
            'h': 'sea_floor_depth_below_sea_level',
            'zeta': 'sea_surface_height',
            'u': 'x_sea_water_velocity',
            'v': 'y_sea_water_velocity',
            'w': 'upward_sea_water_velocity',
            'temp': 'sea_water_temperature',
            'salt': 'sea_water_salinity',
            'uice': 'sea_ice_x_velocity',
            'vice': 'sea_ice_y_velocity',
            'aice': 'sea_ice_area_fraction',
            'hice': 'sea_ice_thickness',
            'gls': 'turbulent_generic_length_scale',
            'tke': 'turbulent_kinetic_energy',
            'AKs': 'ocean_vertical_diffusivity',
            'sustr': 'surface_downward_x_stress',
            'svstr': 'surface_downward_y_stress',
            'Uwind': 'x_wind',
            'Vwind': 'y_wind'}

        # z-levels to which sigma-layers may be interpolated
        self.zlevels = np.array([
            0, -.5, -1, -3, -5, -10, -25, -50, -75, -100, -150, -200,
            -250, -300, -400, -500, -600, -700, -800, -900, -1000, -1500,
            -2000, -2500, -3000, -3500, -4000, -4500, -5000, -5500, -6000,
            -6500, -7000, -7500, -8000])

        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            logging.info('Opening dataset: ' + filestr)
            if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                logging.info('Opening files with MFDataset')
                if has_xarray is True:
                    self.Dataset = xr.open_mfdataset(filename)
                else:
                    self.Dataset = MFDataset(filename)
            else:
                logging.info('Opening file with Dataset')
                if has_xarray is True:
                    self.Dataset = xr.open_dataset(filename)
                else:
                    self.Dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

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
                logging.warning(
                    's_rho not available in dataset, constructing from'
                    ' number of layers (%s).' % num_sigma)
                self.sigma = (np.arange(num_sigma)+.5-num_sigma)/num_sigma

            # Read sigma-coordinate transform parameters
            try:
                self.Dataset.variables['Cs_r'].set_auto_mask(False)
            except:
                pass
            self.Cs_r = self.Dataset.variables['Cs_r'][:]
            try:
                self.hc = self.Dataset.variables['hc'][:]
            except:
                if has_xarray is True:
                    self.hc = self.Dataset.variables['hc'].data  # scalar
                else:
                    self.hc = self.Dataset.variables['hc'][0]

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

        try:  # Check for GLS parameters (diffusivity)
            self.gls_parameters = {}
            for gls_param in ['gls_cmu0', 'gls_p', 'gls_m', 'gls_n']:
                self.gls_parameters[gls_param] = \
                    self.Dataset.variables[gls_param][()]
            logging.info('Read GLS parameters from file.')
        except Exception as e:
            logging.info(e)
            logging.info('Did not find complete set of GLS parameters')

        # Get time coverage
        try:
            ocean_time = self.Dataset.variables['ocean_time']
        except:
            ocean_time = self.Dataset.variables['time']
        if has_xarray:
            self.times = [datetime.utcfromtimestamp((OT -
                          np.datetime64('1970-01-01T00:00:00Z')
                            ) / np.timedelta64(1, 's'))
                          for OT in ocean_time.data]
        else:
            time_units = ocean_time.__dict__['units']
            if time_units == 'second':
                logging.info('Ocean time given as seconds relative to start '
                             'Setting artifical start time of 1 Jan 2000.')
                time_units = 'seconds since 2000-01-01 00:00:00'
            self.times = num2date(ocean_time[:], time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
        else:
            self.time_step = None

        # x and y are rows and columns for unprojected datasets
        self.xmin = 0.
        self.delta_x = 1.
        self.ymin = 0.
        self.delta_y = 1.
        if has_xarray:
            self.xmax = self.Dataset['xi_rho'].shape[0] - 1.
            self.ymax = self.Dataset['eta_rho'].shape[0] - 1.
            self.lon = self.lon.data  # Extract, could be avoided downstream
            self.lat = self.lat.data
            self.sigma = self.sigma.data
        else:
            self.xmax = np.float(len(self.Dataset.dimensions['xi_rho'])) - 1
            self.ymax = np.float(len(self.Dataset.dimensions['eta_rho'])) - 1

        self.name = 'roms native'

        self.precalculate_s2z_coefficients = True

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

        # If one vector component is requested, but not the other
        # we must add the other for correct rotation
        for vector_pair in vector_pairs_xy:
            if (vector_pair[0] in requested_variables and 
                vector_pair[1] not in requested_variables):
                requested_variables.extend([vector_pair[1]])
            if (vector_pair[1] in requested_variables and 
                vector_pair[0] not in requested_variables):
                requested_variables.extend([vector_pair[0]])

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        variables = {}

        if z is None:
            z = np.atleast_1d(0)

        # Find horizontal indices corresponding to requested x and y
        if hasattr(self, 'clipped'):
            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped
        indx[outside] = 0  # To be masked later
        indy[outside] = 0
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

        # Find depth levels covering all elements
        if z.min() == 0 or not hasattr(self, 'hc'):
            indz = self.num_layers - 1  # surface layer
            variables['z'] = 0

        else:
            # Find the range of indices covering given z-values
            if not hasattr(self, 'sea_floor_depth_below_sea_level'):
                logging.debug('Reading sea floor depth...')
                self.sea_floor_depth_below_sea_level = \
                    self.Dataset.variables['h'][:]

                Htot = self.sea_floor_depth_below_sea_level
                self.z_rho_tot = depth.sdepth(Htot, self.hc, self.Cs_r)

            if has_xarray is False:
                indxgrid, indygrid = np.meshgrid(indx, indy)
                H = self.sea_floor_depth_below_sea_level[indygrid, indxgrid]
            else:
                H = self.sea_floor_depth_below_sea_level[indy, indx]
            z_rho = depth.sdepth(H, self.hc, self.Cs_r)
            # Element indices must be relative to extracted subset
            indx_el = indx_el - indx.min()
            indy_el = indy_el - indy.min()

            # Loop to find the layers covering the requested z-values
            indz_min = 0
            indz_max = self.num_layers
            for i in range(self.num_layers):
                if np.min(z-z_rho[i, indy_el, indx_el]) > 0:
                    indz_min = i
                if np.max(z-z_rho[i, indy_el, indx_el]) > 0:
                    indz_max = i
            indz = range(np.maximum(0, indz_min-self.verticalbuffer),
                         np.minimum(self.num_layers,
                                    indz_max + 1 + self.verticalbuffer))
            z_rho = z_rho[indz, :, :]
            # Determine the z-levels to which to interpolate
            zi1 = np.maximum(0, bisect_left(-np.array(self.zlevels),
                                            -z.max()) - self.verticalbuffer)
            zi2 = np.minimum(len(self.zlevels),
                             bisect_right(-np.array(self.zlevels),
                                          -z.min()) + self.verticalbuffer)
            variables['z'] = np.array(self.zlevels[zi1:zi2])

        #read_masks = {}  # To store maskes for various grids
        for par in requested_variables:
            varname = [name for name, cf in
                       self.ROMS_variable_mapping.items() if cf == par]
            var = self.Dataset.variables[varname[0]]

            if has_xarray is not True:
                # Automatic masking may lead to trouble for ROMS files
                # with valid_min/max, _Fill_value or missing_value
                # https://github.com/Unidata/netcdf4-python/issues/703
                var.set_auto_maskandscale(False)

                try:
                    FillValue = getattr(var, '_FillValue')
                except:
                    FillValue = None
                try:
                    scale = getattr(var, 'scale_factor')
                except:
                    scale = 1
                try:
                    offset = getattr(var, 'add_offset')
                except:
                    offset = 0

            if var.ndim == 2:
                variables[par] = var[indy, indx]
            elif var.ndim == 3:
                variables[par] = var[indxTime, indy, indx]
            elif var.ndim == 4:
                variables[par] = var[indxTime, indz, indy, indx]
            else:
                raise Exception('Wrong dimension of variable: ' +
                                self.variable_mapping[par])

            if has_xarray is False:
                # Manual scaling, offsetting and masking due to issue with ROMS files
                logging.debug('Manually masking %s, FillValue %s, scale %s, offset %s' % 
                    (par, FillValue, scale, offset))
                if FillValue is not None:
                    if var.dtype != FillValue.dtype:
                        mask = variables[par] == 0
                        if not 'already_warned' in locals():
                            logging.warning('Data type of variable (%s) and _FillValue (%s) is not the same. Masking 0-values instead' % (var.dtype, FillValue.dtype))
                            already_warned = True
                    else:
                        logging.warning('Masking ' + str(FillValue))
                        mask = variables[par] == FillValue
                variables[par] = variables[par]*scale + offset
                if FillValue is not None:
                    variables[par][mask] = np.nan

            if var.ndim == 4:
                # Regrid from sigma to z levels
                if len(np.atleast_1d(indz)) > 1:
                    logging.debug('sigma to z for ' + varname[0])
                    if self.precalculate_s2z_coefficients is True:
                        M = self.sea_floor_depth_below_sea_level.shape[0]
                        N = self.sea_floor_depth_below_sea_level.shape[1]
                        O = len(self.z_rho_tot)
                        if not hasattr(self, 's2z_A'):
                            logging.info('Calculating sigma2z-coefficients for whole domain')
                            starttime = datetime.now()
                            dummyvar = np.ones((O, M, N))
                            dummy, self.s2z_total = depth.multi_zslice(dummyvar, self.z_rho_tot, self.zlevels)
                            # Store arrays/coefficients
                            self.s2z_A = self.s2z_total[0].reshape(len(self.zlevels), M, N)
                            self.s2z_C = self.s2z_total[1].reshape(len(self.zlevels), M, N)
                            #self.s2z_I = self.s2z_total[2].reshape(M, N)
                            self.s2z_kmax = self.s2z_total[3]
                            del self.s2z_total  # Free memory
                            logging.info('Time: ' + str(datetime.now() - starttime))
                        if 'A' not in locals():
                            logging.info('Re-using sigma2z-coefficients')
                            # Select relevant subset of full arrays
                            zle = np.arange(zi1, zi2)  # The relevant depth levels
                            A = self.s2z_A.copy()  # Awkward subsetting to prevent losing one dimension
                            A = A[:,:,indx]
                            A = A[:,indy,:]
                            A = A[zle,:,:]
                            C = self.s2z_C.copy()
                            C = C[:,:,indx]
                            C = C[:,indy,:]
                            C = C[zle,:,:]
                            C = C - C.max() + variables[par].shape[0] - 1
                            C[C<1] = 1
                            A = A.reshape(len(zle), len(indx)*len(indy))
                            C = C.reshape(len(zle), len(indx)*len(indy))
                            I = np.arange(len(indx)*len(indy))
                            ## Check
                            #dummyvar2, s2z = depth.multi_zslice(
                            #    variables[par].copy(), z_rho.copy(), variables['z'])
                            #print len(zle), variables[par].shape, 'zle, varshape'
                            #Ac,Cc,Ic,kmaxc = s2z
                            #print C, 'C'
                            #print Cc, 'Cc'
                            #print C.shape, Cc.shape
                            #if C.max() != Cc.max():
                            #    print 'WARNING!!'
                            #    import sys; sys.exit('stop')
                            kmax = len(zle)  # Must be checked. Or number of sigma-layers?
                    if 'A' not in locals():
                        logging.info('Calculating new sigma2z-coefficients')
                        variables[par], s2z = depth.multi_zslice(
                            variables[par], z_rho, variables['z'])
                        A,C,I,kmax = s2z
                        # Reshaping to compare with subset of full array
                        #zle = np.arange(zi1, zi2)
                        #A = A.reshape(len(zle), len(indx), len(indy))
                        #C = C.reshape(len(zle), len(indx), len(indy))
                        #I = I.reshape(len(indx), len(indy))
                    else:
                        logging.info('Applying sigma2z-coefficients')
                        # Re-using sigma2z koefficients:
                        F = np.asarray(variables[par])
                        Fshape = F.shape
                        N = F.shape[0]
                        M = F.size // N
                        F = F.reshape((N, M))
                        R = (1-A)*F[(C-1, I)]+A*F[(C, I)]
                        variables[par] = R.reshape((kmax,) + Fshape[1:])

                    # Nan in input to multi_zslice gives extreme values in output
                    variables[par][variables[par]>1e+9] = np.nan

            # If 2D array is returned due to the fancy slicing methods
            # of netcdf-python, we need to take the diagonal
            if variables[par].ndim > 1 and block is False:
                variables[par] = variables[par].diagonal()

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par], ndmin=2, mask=False)
            if block is False:
                variables[par].mask[outside] = True

            # Skipping de-staggering, as it leads to invalid values at later interpolation
            #if block is True:
            #    # Unstagger grid for vectors
            #    logging.debug('Unstaggering ' + par)
            #    if 'eta_v' in var.dimensions:
            #        variables[par] = np.ma.array(variables[par],
            #                            mask=variables[par].mask)
            #        variables[par][variables[par].mask] = 0
            #        if variables[par].ndim == 2:
            #            variables[par] = \
            #                (variables[par][0:-1,0:-1] +
            #                variables[par][0:-1,1::])/2
            #        elif variables[par].ndim == 3:
            #            variables[par] = \
            #                (variables[par][:,0:-1,0:-1] +
            #                variables[par][:,0:-1,1::])/2
            #        variables[par] = np.ma.masked_where(variables[par]==0,
            #                                            variables[par])
            #    elif 'eta_u' in var.dimensions:
            #        variables[par] = np.ma.array(variables[par],
            #                            mask=variables[par].mask)
            #        variables[par][variables[par].mask] = 0
            #        if variables[par].ndim == 2:
            #            variables[par] = \
            #                (variables[par][0:-1,0:-1] +
            #                 variables[par][1::,0:-1])/2
            #        elif variables[par].ndim == 3:
            #            variables[par] = \
            #                (variables[par][:,0:-1,0:-1] +
            #                 variables[par][:,1::,0:-1])/2
            #        variables[par] = np.ma.masked_where(variables[par]==0,
            #                                            variables[par])
            #    else:
            #        if variables[par].ndim == 2:
            #            variables[par] = variables[par][1::, 1::]
            #        elif variables[par].ndim == 3:
            #            variables[par] = variables[par][:,1::, 1::]

        if block is True:
            # TODO: should be midpoints, but angle array below needs integer
            #indx = indx[0:-1]  # Only if de-staggering has been performed
            #indy = indy[1::]
            variables['x'] = indx
            variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['x'] = variables['x'].astype(np.float)
        variables['y'] = variables['y'].astype(np.float)
        variables['time'] = nearestTime

        if 'x_sea_water_velocity' or 'sea_ice_x_velocity' \
                or 'x_wind' in variables.keys():
            # We must rotate current vectors
            if not hasattr(self, 'angle_xi_east'):
                logging.debug('Reading angle between xi and east...')
                self.angle_xi_east = self.Dataset.variables['angle'][:]
            if has_xarray is False:
                rad = self.angle_xi_east[np.meshgrid(indy, indx)].T
            else:
                rad = self.angle_xi_east[indy, indx].T
            if 'x_sea_water_velocity' in variables.keys():
                variables['x_sea_water_velocity'], \
                    variables['y_sea_water_velocity'] = rotate_vectors_angle(
                        variables['x_sea_water_velocity'],
                        variables['y_sea_water_velocity'], rad)
            if 'sea_ice_x_velocity' in variables.keys():
                variables['sea_ice_x_velocity'], \
                    variables['sea_ice_y_velocity'] = rotate_vectors_angle(
                        variables['sea_ice_x_velocity'],
                        variables['sea_ice_y_velocity'], rad)
            if 'x_wind' in variables.keys():
                variables['x_wind'], \
                    variables['y_wind'] = rotate_vectors_angle(
                        variables['x_wind'],
                        variables['y_wind'], rad)

        if 'land_binary_mask' in requested_variables:
            variables['land_binary_mask'] = \
                1 - variables['land_binary_mask']

        # Masking NaN
        for var in requested_variables:
            variables[var] = np.ma.masked_invalid(variables[var])
        
        return variables


def rotate_vectors_angle(u, v, radians):
    u2 = u*np.cos(radians) - v*np.sin(radians)
    v2 = u*np.sin(radians) + v*np.cos(radians)
    return u2, v2
