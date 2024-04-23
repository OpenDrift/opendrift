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

from bisect import bisect_left, bisect_right
from datetime import datetime
import logging
logger = logging.getLogger(__name__)

import numpy as np
from netCDF4 import num2date
import xarray as xr

from opendrift.readers.basereader import BaseReader, vector_pairs_xy, StructuredReader
from opendrift.readers.roppy import depth


class Reader(BaseReader, StructuredReader):
    """
    A reader for ROMS Output files. It can take a single file, a file pattern, a URL or an xarray Dataset.

    Args:
        :param filename: A single netCDF file, a pattern of files, or a xr.Dataset. The
                         netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, xr.Dataset (required).

        :param name: Name of reader
        :type name: string, optional

        :param save_interpolator: Whether or not to save the interpolator that goes from lon/lat to x/y (calculated in structured.py)
        :type save_interpolator: bool

        :param interpolator_filename: If save_interpolator is True, user can input this string to control where interpolator is saved.
        :type interpolator_filename: Path, str, optional

    Example:

    .. code::

       from opendrift.readers.reader_ROMS_native import Reader
       r = Reader("roms.nc")

    Several files can be specified by using a pattern:

    .. code::

       from opendrift.readers.reader_ROMS_native import Reader
       r = Reader("*.nc")

    An OPeNDAP URL can be used:

    .. code::

       from opendrift.readers.reader_ROMS_native import Reader
       r = Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')

    A xr.Dataset can be used:

    .. code::

        from opendrift.readers.reader_ROMS_native import Reader
        ds = xr.open_dataset(filename, decode_times=False)
        r = Reader(ds)
    """

    def __init__(self, filename=None, name=None, gridfile=None, standard_name_mapping={},
                 save_interpolator=False, interpolator_filename=None):
        
        self._mask_rho = None
        self._mask_u = None
        self._mask_v = None
        self._zeta = None
        self._angle = None
        self.land_binary_mask = None
        self.sea_floor_depth_below_sea_level = None
        self.z_rho_tot = None
        self.s2z_A = None

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        # Map ROMS variable names to CF standard_name
        self.ROMS_variable_mapping = {
            # Removing (temoprarily) land_binary_mask from ROMS-variables,
            # as this leads to trouble with linearNDFast interpolation
            'mask_rho': 'land_binary_mask',
            # 'mask_psi': 'land_binary_mask',  # don't want two variables mapping together - raises error now
            'h': 'sea_floor_depth_below_sea_level',
            'zeta': 'sea_surface_height',
            'u': 'x_sea_water_velocity',
            'v': 'y_sea_water_velocity',
            #'u_eastward': 'x_sea_water_velocity',  # these are wrognly rotated below
            #'v_northward': 'y_sea_water_velocity',
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
            'tair': 'air_temperature',
            'wspd': 'wind_speed',
            'uwnd': 'x_wind',
            'vwnd': 'y_wind',
            'uwind': 'x_wind',
            'vwind': 'y_wind',
            'Uwind': 'x_wind',
            'Vwind': 'y_wind',
            }

        # Add user provided variable mappings
        self.ROMS_variable_mapping.update(standard_name_mapping)

        # z-levels to which sigma-layers may be interpolated
        self.zlevels = np.array([
            0, -.5, -1, -3, -5, -10, -25, -50, -75, -100, -150, -200,
            -250, -300, -400, -500, -600, -700, -800, -900, -1000, -1500,
            -2000, -2500, -3000, -3500, -4000, -4500, -5000, -5500, -6000,
            -6500, -7000, -7500, -8000])

        gls_param = ['gls_cmu0', 'gls_p', 'gls_m', 'gls_n']
        
        self.name = name or 'roms native'

        if isinstance(filename, xr.Dataset):
            self.Dataset = filename
        else:

            filestr = str(filename)

            try:
                # Open file, check that everything is ok
                logger.info('Opening dataset: ' + filestr)
                if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
                    logger.info('Opening files with MFDataset')
                    def drop_non_essential_vars_pop(ds):
                        dropvars = [v for v in ds.variables if v not in
                                    list(self.ROMS_variable_mapping.keys()) + gls_param +
                                    ['ocean_time', 'time', 'bulk_time', 's_rho',
                                     'Cs_r', 'hc', 'angle', 'Vtransform']
                                    and v[0:3] not in ['lon', 'lat', 'mas']]
                        logger.debug('Dropping variables: %s' % dropvars)
                        ds = ds.drop_vars(dropvars)
                        return ds
                    self.Dataset = xr.open_mfdataset(filename,
                        chunks={'ocean_time': 1}, compat='override', decode_times=False,
                        preprocess=drop_non_essential_vars_pop,
                        data_vars='minimal', coords='minimal')
                else:
                    logger.info('Opening file with Dataset')
                    self.Dataset = xr.open_dataset(filename, decode_times=False)
            except Exception as e:
                raise ValueError(e)

        # this is an opporunity to save interpolators to pickle to save sim time
        self.save_interpolator = save_interpolator
        self.interpolator_filename = interpolator_filename or f'{self.name}_interpolators'

        if gridfile is not None:  # Merging gridfile dataset with main dataset
            gf = xr.open_dataset(gridfile)
            self.Dataset = xr.merge([self.Dataset, gf])

        if 'Vtransform' in self.Dataset.variables:
            self.Vtransform = self.Dataset.variables['Vtransform'].data  # scalar
        else:
            logger.warning('Vtransform not found, using 1')
            self.Vtransform = 1
        self.Vtransform = np.asarray(self.Vtransform)

        if 's_rho' not in self.Dataset.variables:
            dimensions = 2
        else:
            dimensions = 3

        if dimensions == 3:
            # Read sigma-coordinate values
            try:
                self.sigma = self.Dataset.variables['s_rho'][:]
                self.sigma = self.sigma.data  # Extract data
            except:
                num_sigma = len(self.Dataset.dimensions['s_rho'])
                logger.warning(
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
                self.hc = self.Dataset.variables['hc'].data  # scalar
            else:
                self.hc = None

            self.num_layers = len(self.sigma)
        else:
            logger.warning("2D dataset, so deleting u and v from ROMS_variable_mapping")
            self.num_layers = 1
            self.ROMS_variable_mapping['ubar'] = 'x_sea_water_velocity'
            self.ROMS_variable_mapping['vbar'] = 'y_sea_water_velocity'
            del self.ROMS_variable_mapping['u']
            del self.ROMS_variable_mapping['v']

        if 'lat_rho' in self.Dataset.variables:
            # Horizontal oordinates and directions
            self.lat = self.Dataset.variables['lat_rho'][:]
            self.lon = self.Dataset.variables['lon_rho'][:]
            self.lon = self.lon.data  # Extract, could be avoided downstream
            self.lat = self.lat.data
            if self.lat.ndim == 1:
                self.lon, self.lat = np.meshgrid(self.lon, self.lat)
                # self.angle_xi_east = 0  # this was moved to the angle property
        else:
            raise ValueError(filename + ' does not contain lon/lat '
                             'arrays, please supply a grid-file: "gridfile=<grid_file>"')

        for var in list(self.ROMS_variable_mapping):  # Remove unused variables
            if var not in self.Dataset.variables:
                del self.ROMS_variable_mapping[var]

        try:  # Check for GLS parameters (diffusivity)
            self.gls_parameters = {}
            for gls_par in gls_param:
                self.gls_parameters[gls_par] = \
                    self.Dataset.variables[gls_par][()]
            logger.info('Read GLS parameters from file.')
        except Exception as e:
            logger.info(e)
            logger.info('Did not find complete set of GLS parameters')

        # Get time coverage
        ocean_time = None
        for tv in ['ocean_time', 'time', 'bulk_time']:
            if tv in self.Dataset.variables:
                ocean_time = self.Dataset.variables[tv]
        if ocean_time is None:
            raise ValueError('Time variable not found in ROMS file')
        time_units = ocean_time.attrs['units']
        if time_units == 'second':
            logger.info('Ocean time given as seconds relative to start '
                         'Setting artifical start time of 1 Jan 2000.')
            time_units = 'seconds since 2000-01-01 00:00:00'
        if time_units == 'days':
            logger.info('Ocean time given as days relative to start '
                         'Setting artifical start time of 1 Jan 2000.')
            time_units = 'days since 2000-01-01 00:00:00'
        self.times = num2date(ocean_time[:], time_units)
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
        else:
            self.time_step = None

        self.precalculate_s2z_coefficients = True

        # Find all variables having standard_name
        self.variables = []
        for var_name in list(self.Dataset.variables):
            var = self.Dataset.variables[var_name]
            if 'standard_name' in var.attrs and var_name not in self.ROMS_variable_mapping.keys():
                self.ROMS_variable_mapping[var_name] = var.attrs['standard_name']
            if var_name in self.ROMS_variable_mapping.keys():
                self.variables.append(self.ROMS_variable_mapping[var_name])

        # A bit hackish solution:
        # If variable names or their standard_name contain "east" or "north", 
        # these should not be rotated from xi-direction to east-direction
        self.do_not_rotate = []
        for var, stdname in self.ROMS_variable_mapping.items():
            if 'east' in var.lower() or 'east' in stdname.lower() or \
                    'north' in var.lower() or 'north' in stdname.lower():
                self.do_not_rotate.append(stdname)
        if len(self.do_not_rotate)>0:
            logger.debug('The following ROMS vectors are considered east-north, and will not be rotated %s' % self.do_not_rotate)

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    @property
    def mask_rho(self):
        """Mask for the rho-points.
        
        Uses wetdry_mask_rho (which should be 3D) if available, otherwise mask_rho (2D).
        If this mask is 2D, read it in this one time and use going forward in simulation. If 3D,
        will read in parts of the mask each loop.
        """
        if self._mask_rho is None:
            if 'wetdry_mask_rho' in self.Dataset.data_vars:
                self._mask_rho = self.Dataset.variables['wetdry_mask_rho']
                logger.info("Using wetdry_mask_rho for mask_rho")
            else:
                # Read landmask for whole domain, for later re-use
                self._mask_rho = self.Dataset.variables['mask_rho'][:]
                # load in once if static mask
                if self._mask_rho.chunks is not None:  # to see if dask array
                    self._mask_rho = self._mask_rho.compute()
                logger.info("Using mask_rho for mask_rho")
        return self._mask_rho

    @property
    def mask_u(self):
        """Mask for the u-points.
        
        Uses wetdry_mask_u (which should be 3D) if available, otherwise mask_u (2D).
        If this mask is 2D, read it in this one time and use going forward in simulation. If 3D,
        will read in parts of the mask each loop.
        """
        if self._mask_u is None:
            if 'wetdry_mask_u' in self.Dataset.data_vars:
                self._mask_u = self.Dataset.variables['wetdry_mask_u']
                logger.info("Using wetdry_mask_u for mask_u")
            else:
                # Read landmask for whole domain, for later re-use
                self._mask_u = self.Dataset.variables['mask_u'][:]
                # load in once if static mask
                if self._mask_u.chunks is not None:  # to see if dask array
                    self._mask_u = self._mask_u.compute()
                logger.info("Using mask_u for mask_u")
        return self._mask_u

    @property
    def mask_v(self):
        """Mask for the v-points.
        
        Uses wetdry_mask_v (which should be 3D) if available, otherwise mask_v (2D).
        If this mask is 2D, read it in this one time and use going forward in simulation. If 3D,
        will read in parts of the mask each loop.
        """
        if self._mask_v is None:
            if 'wetdry_mask_v' in self.Dataset.data_vars:
                self._mask_v = self.Dataset.variables['wetdry_mask_v']
                logger.info("Using wetdry_mask_v for mask_v")
            else:
                # Read landmask for whole domain, for later re-use
                self._mask_v = self.Dataset.variables['mask_v'][:]
                # load in once if static mask
                if self._mask_v.chunks is not None:  # to see if dask array
                    self._mask_v = self._mask_v.compute()
                logger.info("Using mask_v for mask_v")
        return self._mask_v
    
    @property
    def zeta(self):
        """Sea surface height."""
        if self._zeta is None:
            if 'zeta' in self.Dataset.data_vars:
                self._zeta = self.Dataset.variables['zeta']
                logger.info("Using zeta for sea surface height")
            else:
                self._zeta = np.zeros(self.mask_rho.shape)
                logger.info("No zeta found, using 0 array for sea surface height")
        return self._zeta        
    
    @property
    def angle(self):
        """Grid angle if curvilinear."""
        if self._angle is None:
            if 'lat_rho' in self.Dataset.variables and self.lat.ndim == 1:
                self._angle = 0
            elif 'angle' in self.Dataset.data_vars:
                self._angle = self.Dataset.variables['angle']
                logger.info("Using angle from Dataset.")
            # else:
            #     self._angle = 0
            #     logger.warning("No angle found, using 0 integer for angle.")

            if self._angle is None:
                raise ValueError('No angle between xi and east found')
        return self._angle

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, testing=False):
        start_time = datetime.now()
        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        # land_binary_mask should be based on the rho grid
        if 'land_binary_mask' in requested_variables and self.land_binary_mask is None:
            self.land_binary_mask = 1 - self.mask_rho

        if 'sea_floor_depth_below_sea_level' in requested_variables and self.sea_floor_depth_below_sea_level is None:
            self.sea_floor_depth_below_sea_level = self.Dataset.variables['h'][:]

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
        if self.clipped is not None:
            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped
        indx[outside] = 0  # To be masked later
        indy[outside] = 0
        indx_el = indx.copy()
        indy_el = indy.copy()
        # Adding buffer, to cover also future positions of elements
        buffer = self.buffer
        # Avoiding the last pixel in each dimension, since there are
        # several grids which are shifted (rho, u, v, psi)
        indx = np.arange(np.max([0, indx.min()-buffer]),
                            np.min([indx.max()+buffer, self.lon.shape[1]-1]))
        indy = np.arange(np.max([0, indy.min()-buffer]),
                            np.min([indy.max()+buffer, self.lon.shape[0]-1]))

        # define indices
        ixy = (indy,indx)
        itxy = (indxTime,indy,indx)
        
        # use these same indices for all mask subsetting since if one is
        # 3D they all should be
        if self.mask_rho.ndim == 2:
            imask = ixy
        elif self.mask_rho.ndim == 3:
            imask = itxy

        # Find depth levels covering all elements
        if z.min() == 0 or self.hc is None:
            indz = self.num_layers - 1  # surface layer
            variables['z'] = 0

        else:
            # Find the range of indices covering given z-values
            if self.sea_floor_depth_below_sea_level is None:
                logger.debug('Reading sea floor depth...')
                self.sea_floor_depth_below_sea_level = \
                    self.Dataset.variables['h'][:]

            if self.z_rho_tot is None:
                Htot = self.sea_floor_depth_below_sea_level
                zeta = self.zeta[indxTime]
                self.z_rho_tot = depth.sdepth(Htot, zeta, self.hc, self.Cs_r,
                                              Vtransform=self.Vtransform)
                # z_rho is positive relative to mean sea level but z is
                # 0 at the surface.
                # Transform z_rho to match convention of z.
                self.z_rho_tot -= np.asarray(zeta)[np.newaxis]

            H = self.sea_floor_depth_below_sea_level[indy, indx]
            zeta = self.zeta[itxy]
            z_rho = depth.sdepth(H, zeta, self.hc, self.Cs_r,
                                 Vtransform=self.Vtransform)

            # z_rho is positive relative to mean sea level but z is
            # 0 at the surface.
            # Transform z_rho to match convention of z.
            z_rho -= np.asarray(zeta)[np.newaxis]

            # Check for positive values in z_rho
            # if there are any, nan them out since they are above
            # the surface and therefore the cell is dry
            # this should only come up for wet/dry simulations
            if (np.nanmax(z_rho) > 0).any():
                logger.info(f'z_rho had positive values that are now nans.')
                z_rho[z_rho>0] = np.nan

            # Element indices must be relative to extracted subset
            indx_el = np.clip(indx_el - indx.min(), 0, z_rho.shape[2]-1)
            indy_el = np.clip(indy_el - indy.min(), 0, z_rho.shape[1]-1)

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
        
        # define another set of indices
        itzxy = (indxTime, indz, indy, indx)
            
        def get_mask(mask_name, imask, masks_store):
            if mask_name in masks_store:
                mask = masks_store[mask_name]
            else:
                mask = getattr(self, mask_name)[imask]
            return mask, mask_name

        masks_store = {}  # To store masks for various grids
        for par in requested_variables:
            varname = [name for name, cf in
                       self.ROMS_variable_mapping.items() if cf == par]
            if len(varname) > 1:
                raise ValueError("Multiple variables exist with standard name and "
                                 "are present in reader. Either remove the duplicate mapping "
                                 "or remove the variable from the reader."
                                 "Variables: " + str(varname))
            var = self.Dataset.variables[varname[0]]

            if par == 'land_binary_mask':
                variables[par] = self.land_binary_mask[imask]
            elif par == 'sea_floor_depth_below_sea_level':
                variables[par] = self.sea_floor_depth_below_sea_level[ixy]
            elif var.ndim == 2:
                variables[par] = var[ixy]
            elif var.ndim == 3:
                variables[par] = var[itxy]
            elif var.ndim == 4:
                variables[par] = var[itzxy]
            else:
                raise Exception('Wrong dimension of variable: ' +
                                self.ROMS_variable_mapping[par])

            # If Xarray, load in so can be used in future loop iterations too
            variables[par] = np.asarray(variables[par])            

            if par != 'land_binary_mask':

                # make sure that var has matching horizontal dimensions with the mask
                # make sure coord names also match
                if self.mask_rho.shape[-2:] == var.shape[-2:] and self.mask_rho.dims[-2:] == var.dims[-2:]:
                    mask, mask_name = get_mask("mask_rho", imask, masks_store)
                elif self.mask_u.shape[-2:] == var.shape[-2:] and self.mask_u.dims[-2:] == var.dims[-2:]:
                    mask, mask_name = get_mask("mask_u", imask, masks_store)
                elif self.mask_v.shape[-2:] == var.shape[-2:] and self.mask_v.dims[-2:] == var.dims[-2:]:
                    mask, mask_name = get_mask("mask_v", imask, masks_store)
                else:
                    raise Exception('No mask found for ' + par)

                masks_store[mask_name] = np.asarray(mask)
                mask = np.asarray(mask)

                if mask.min() == 0:
                    # Not using the fill value directly allows the mask to have more control
                    # which is necessary when using a wetdry mask that changes in time
                    # since the fill value will not cover all masked locations then.
                    variables[par][...,mask==0] = np.nan

            if var.ndim == 4:
                # Regrid from sigma to z levels
                if len(np.atleast_1d(indz)) > 1:
                    logger.debug('sigma to z for ' + varname[0])
                    if self.precalculate_s2z_coefficients is True:
                        M = self.sea_floor_depth_below_sea_level.shape[0]
                        N = self.sea_floor_depth_below_sea_level.shape[1]
                        O = len(self.z_rho_tot)
                        if self.s2z_A is None:
                            logger.debug('Calculating sigma2z-coefficients for whole domain')
                            starttime = datetime.now()
                            dummyvar = np.ones((O, M, N))
                            dummy, self.s2z_total = depth.multi_zslice(dummyvar, self.z_rho_tot, self.zlevels)
                            # Store arrays/coefficients
                            self.s2z_A = self.s2z_total[0].reshape(len(self.zlevels), M, N)
                            self.s2z_C = self.s2z_total[1].reshape(len(self.zlevels), M, N)
                            #self.s2z_I = self.s2z_total[2].reshape(M, N)
                            self.s2z_kmax = self.s2z_total[3]
                            del self.s2z_total  # Free memory
                            logger.info('Time: ' + str(datetime.now() - starttime))
                        if 'A' not in locals():
                            logger.debug('Re-using sigma2z-coefficients')
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
                        logger.debug('Calculating new sigma2z-coefficients')
                        variables[par], s2z = depth.multi_zslice(
                            variables[par], z_rho, variables['z'])
                        A,C,I,kmax = s2z
                        # Reshaping to compare with subset of full array
                        #zle = np.arange(zi1, zi2)
                        #A = A.reshape(len(zle), len(indx), len(indy))
                        #C = C.reshape(len(zle), len(indx), len(indy))
                        #I = I.reshape(len(indx), len(indy))
                    else:
                        logger.debug('Applying sigma2z-coefficients')
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

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par], ndmin=2, mask=False)

            # Skipping de-staggering, as it leads to invalid values at later interpolation
            #if block is True:
            #    # Unstagger grid for vectors
            #    logger.debug('Unstaggering ' + par)
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

        # TODO: should be midpoints, but angle array below needs integer
        #indx = indx[0:-1]  # Only if de-staggering has been performed
        #indy = indy[1::]
        variables['x'] = indx
        variables['y'] = indy

        variables['x'] = variables['x'].astype(np.float32)
        variables['y'] = variables['y'].astype(np.float32)
        variables['time'] = nearestTime

        if 'x_sea_water_velocity' in variables.keys() and \
            'x_sea_water_velocity' not in self.do_not_rotate:
            if isinstance(self.angle, int):
                rad = self.angle
            else:
                rad = np.ma.asarray(self.angle[indy, indx])
            variables['x_sea_water_velocity'], \
                variables['y_sea_water_velocity'] = rotate_vectors_angle(
                    variables['x_sea_water_velocity'],
                    variables['y_sea_water_velocity'], rad)
            logger.debug('Rotated x_sea_water_velocity and y_sea_water_velocity')
        
        if 'sea_ice_x_velocity' in variables.keys() and \
            'sea_ice_x_velocity' not in self.do_not_rotate:
            if isinstance(self.angle, int):
                rad = self.angle
            else:
                rad = np.ma.asarray(self.angle[indy, indx])
            if 'sea_ice_x_velocity' in variables.keys() and \
                    'sea_ice_x_velocity' not in self.do_not_rotate:
                variables['sea_ice_x_velocity'], \
                    variables['sea_ice_y_velocity'] = rotate_vectors_angle(
                        variables['sea_ice_x_velocity'],
                        variables['sea_ice_y_velocity'], rad)
                logger.debug('Rotated sea_ice_x_velocity and sea_ice_y_velocity')
        
        if 'x_wind' in variables.keys() and 'x_wind' not in self.do_not_rotate:
            if isinstance(self.angle, int):
                rad = self.angle
            else:
                rad = np.ma.asarray(self.angle[indy, indx])
            if 'x_wind' in variables.keys() and \
                    'x_wind' not in self.do_not_rotate:
                variables['x_wind'], \
                    variables['y_wind'] = rotate_vectors_angle(
                        variables['x_wind'],
                        variables['y_wind'], rad)
                logger.debug('Rotated x_wind and y_wind')

        # Masking NaN
        for var in requested_variables:
            variables[var] = np.ma.masked_invalid(variables[var])

        logger.debug('Time for ROMS native reader: ' + str(datetime.now()-start_time))

        if testing:
            return variables, masks_store
        else:
            return variables


def rotate_vectors_angle(u, v, radians):
    u2 = u*np.cos(radians) - v*np.sin(radians)
    v2 = u*np.sin(radians) + v*np.cos(radians)
    return u2, v2
