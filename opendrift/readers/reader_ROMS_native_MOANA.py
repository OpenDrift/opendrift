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
# same as reader_ROMS_native.py but handles longitude >180, and has some sligtly different variables
# name mapping
# 
# updated 27/11/2020 to reflect changes in "official reader_ROMS_native.py"
# 
# Developed by S.Weppe MetOcean Solutions / MetServices NZ
########################################################################

from bisect import bisect_left, bisect_right
from datetime import datetime

import numpy as np
from netCDF4 import num2date
import xarray as xr

from opendrift.readers.basereader import BaseReader, vector_pairs_xy
from opendrift.readers.roppy import depth
from opendrift.readers.interpolation import ReaderBlock


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
            'u_eastward' : 'x_sea_water_velocity', 
            'v_northward' : 'y_sea_water_velocity',
            'u': 'x_sea_water_velocity',
            'v': 'y_sea_water_velocity',
            'u_eastward': 'x_sea_water_velocity',
            'v_northward': 'y_sea_water_velocity',
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

        gls_param = ['gls_cmu0', 'gls_p', 'gls_m', 'gls_n']

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
                def drop_non_essential_vars_pop(ds):
                    dropvars = [v for v in ds.variables if v not in
                                list(self.ROMS_variable_mapping.keys()) + gls_param +
                                ['ocean_time', 's_rho', 'Cs_r', 'hc', 'angle']
                                and v[0:3] not in ['lon', 'lat', 'mas']]
                    self.logger.debug('Dropping variables: %s' % dropvars)
                    ds = ds.drop(dropvars)
                    return ds
                self.Dataset = xr.open_mfdataset(filename,
                    chunks={'ocean_time': 1}, concat_dim='ocean_time',
                    combine='by_coords',
                    compat='override',
                    decode_times=False,
                    preprocess=drop_non_essential_vars_pop,
                    data_vars='minimal', coords='minimal')
            else:
                self.logger.info('Opening file with Dataset')
                self.Dataset = xr.open_dataset(filename, decode_times=False)
        except Exception as e:
            raise ValueError(e)


        if 'Vtransform' in self.Dataset.variables:
            self.Vtransform = self.Dataset.variables['Vtransform'].data  # scalar
        else:
            self.logger.warning('Vtransform not found, using 1')
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
            except:
                num_sigma = len(self.Dataset.dimensions['s_rho'])
                self.logger.warning(
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

            self.num_layers = len(self.sigma)
        else:
            self.num_layers = 1
            self.ROMS_variable_mapping['ubar'] = 'x_sea_water_velocity'
            self.ROMS_variable_mapping['vbar'] = 'y_sea_water_velocity'
            del self.ROMS_variable_mapping['u']
            del self.ROMS_variable_mapping['v']

        if 'u_eastward' in self.Dataset.variables:
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

        try:  # Check for GLS parameters (diffusivity)
            self.gls_parameters = {}
            for gls_par in gls_param:
                self.gls_parameters[gls_par] = \
                    self.Dataset.variables[gls_par][()]
            self.logger.info('Read GLS parameters from file.')
        except Exception as e:
            self.logger.info(e)
            self.logger.info('Did not find complete set of GLS parameters')

        # Get time coverage
        try:
            ocean_time = self.Dataset.variables['ocean_time']
        except:
            ocean_time = self.Dataset.variables['time']
        time_units = ocean_time.attrs['units']
        if time_units == 'second':
            self.logger.info('Ocean time given as seconds relative to start '
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
        self.xmax = self.Dataset['xi_rho'].shape[0] - 1.
        self.ymax = self.Dataset['eta_rho'].shape[0] - 1.
        self.lon = self.lon.data  # Extract, could be avoided downstream
        self.lat = self.lat.data
        self.sigma = self.sigma.data

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

    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader.

           Overloads covers_positions() from basereader.py
           to account for netcdf files where 0<longitude<360

        """

        z = np.atleast_1d(z)
        # Calculate x,y coordinates from lon,lat
        x, y = self.lonlat2xy(lon, lat)

        # Only checking vertical coverage if zmin, zmax is defined
        zmin = -np.inf
        zmax = np.inf
        if hasattr(self, 'zmin') and self.zmin is not None:
            zmin = self.zmin
        if hasattr(self, 'zmax') and self.zmax is not None:
            zmax = self.zmax

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        else:
            # at this stage 'lon' follows the convention -180<lon<180
            # need to account for reader(s) that may use 0<lon<360 instead
            if not self.has_lon_0_360 : # simple case - no overlap of 180W, use original code
                indices = np.where((x >= self.xmin) & (x <= self.xmax) &
                                   (y >= self.ymin) & (y <= self.ymax) &
                                   (z >= zmin) & (z <= zmax))[0]
            elif self.has_lon_0_360: #reader longitude using convention :  0<lon<360, or has longitude between 0 and +180 only
                x360 = self.to_longitude_0_360(x) # convert opendrift x ([-180,180]) to [0,360] convention
                indices = np.where((x360 >= self.xmin) & (x360 <= self.xmax) &
                                   (y >= self.ymin) & (y <= self.ymax) &
                                   (z >= zmin) & (z <= zmax))[0]
        try:
            return indices, x[indices], y[indices]
        except:
            return indices, x, y


    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):

        start_time = datetime.now()
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
        ##########################################################################################################
        # handling cases when longitude > 180 
        # indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped # original code
        if self.has_lon_0_360 :
            indx = np.floor((self.to_longitude_0_360(x)-self.xmin)/self.delta_x).astype(int) 
        else:
            indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped # original code
        ##########################################################################################################
        # indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped
        indx[outside] = 0  # To be masked later
        indy[outside] = 0
        indx_el = indx.copy()
        indy_el = indy.copy()
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
                self.logger.debug('Reading sea floor depth...')
                self.sea_floor_depth_below_sea_level = \
                    self.Dataset.variables['h'][:]

                Htot = self.sea_floor_depth_below_sea_level
                self.z_rho_tot = depth.sdepth(Htot, self.hc, self.Cs_r,
                                              Vtransform=self.Vtransform)

            H = self.sea_floor_depth_below_sea_level[indy, indx]
            z_rho = depth.sdepth(H, self.hc, self.Cs_r,
                                 Vtransform=self.Vtransform)
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

        #read_masks = {}  # To store maskes for various grids
        mask_values = {}
        for par in requested_variables:
            varname = [name for name, cf in
                       self.ROMS_variable_mapping.items() if cf == par]
            var = self.Dataset.variables[varname[0]]

            if par == 'land_binary_mask':
                if not hasattr(self, 'land_binary_mask'):
                    # Read landmask for whole domain, for later re-use
                    self.land_binary_mask = \
                        1 - self.Dataset.variables['mask_rho'][:]
                variables[par] = self.land_binary_mask[indy, indx]
            elif var.ndim == 2:
                variables[par] = var[indy, indx]
            elif var.ndim == 3:
                variables[par] = var[indxTime, indy, indx]
            elif var.ndim == 4:
                variables[par] = var[indxTime, indz, indy, indx]
            else:
                raise Exception('Wrong dimension of variable: ' +
                                self.variable_mapping[par])

            variables[par] = np.asarray(variables[par])  # If Xarray
            start = datetime.now()

            if par not in mask_values:
                indxgrid = indx
                indygrid = indy
                if par == 'x_sea_water_velocity':
                    if not hasattr(self, 'mask_u'):
                        if 'mask_u' in self.Dataset.variables:
                            self.mask_u = self.Dataset.variables['mask_u'][:]
                        else:
                            self.mask_u = self.Dataset.variables['mask_rho'][:]
                    mask = self.mask_u[indygrid, indxgrid]
                elif par == 'y_sea_water_velocity':
                    if not hasattr(self, 'mask_v'):
                        if 'mask_v' in self.Dataset.variables:
                            self.mask_v = self.Dataset.variables['mask_v'][:]
                        else:
                            self.mask_v = self.Dataset.variables['mask_rho'][:]
                    mask = self.mask_v[indygrid, indxgrid]
                else:
                    if not hasattr(self, 'mask_rho'):
                        # For ROMS-Agrif this must perhaps be mask_psi?
                        self.mask_rho = self.Dataset.variables['mask_rho'][:]
                    mask = self.mask_rho[indygrid, indxgrid]
                mask = np.asarray(mask)
                if mask.min() == 0 and par != 'land_binary_mask':
                    first_mask_point = np.where(mask.ravel()==0)[0][0]
                    if variables[par].ndim == 3:
                        upper = variables[par][0,:,:]
                    else:
                        upper = variables[par]
                    mask_values[par] = upper.ravel()[first_mask_point]
                    variables[par][variables[par]==mask_values[par]] = np.nan

            if var.ndim == 4:
                # Regrid from sigma to z levels
                if len(np.atleast_1d(indz)) > 1:
                    self.logger.debug('sigma to z for ' + varname[0])
                    if self.precalculate_s2z_coefficients is True:
                        M = self.sea_floor_depth_below_sea_level.shape[0]
                        N = self.sea_floor_depth_below_sea_level.shape[1]
                        O = len(self.z_rho_tot)
                        if not hasattr(self, 's2z_A'):
                            self.logger.debug('Calculating sigma2z-coefficients for whole domain')
                            starttime = datetime.now()
                            dummyvar = np.ones((O, M, N))
                            dummy, self.s2z_total = depth.multi_zslice(dummyvar, self.z_rho_tot, self.zlevels)
                            # Store arrays/coefficients
                            self.s2z_A = self.s2z_total[0].reshape(len(self.zlevels), M, N)
                            self.s2z_C = self.s2z_total[1].reshape(len(self.zlevels), M, N)
                            #self.s2z_I = self.s2z_total[2].reshape(M, N)
                            self.s2z_kmax = self.s2z_total[3]
                            del self.s2z_total  # Free memory
                            self.logger.info('Time: ' + str(datetime.now() - starttime))
                        if 'A' not in locals():
                            self.logger.debug('Re-using sigma2z-coefficients')
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
                        self.logger.debug('Calculating new sigma2z-coefficients')
                        variables[par], s2z = depth.multi_zslice(
                            variables[par], z_rho, variables['z'])
                        A,C,I,kmax = s2z
                        # Reshaping to compare with subset of full array
                        #zle = np.arange(zi1, zi2)
                        #A = A.reshape(len(zle), len(indx), len(indy))
                        #C = C.reshape(len(zle), len(indx), len(indy))
                        #I = I.reshape(len(indx), len(indy))
                    else:
                        self.logger.debug('Applying sigma2z-coefficients')
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
            #    self.logger.debug('Unstaggering ' + par)
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
                self.logger.debug('Reading angle between xi and east...')
                self.angle_xi_east = self.Dataset.variables['angle'][:]
            rad = self.angle_xi_east[indy, indx]
            rad = np.ma.asarray(rad)
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

        # Masking NaN
        for var in requested_variables:
            variables[var] = np.ma.masked_invalid(variables[var])
        
        self.logger.debug('Time for ROMS native reader: ' + str(datetime.now()-start_time))

        return variables

    def get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):

        # this will overload the get_variables_interpolated() from basereader.py
        # to allow applying log profile correction, adding readers etc..
        # eventually this may be a nice feature to add the basereader ?

        self.timer_start('total')
        self.timer_start('preparing')

        #############################################################################
        # Changed ;  
        # convert lon,lat positions to fit with reader conventions i.e. 0<lon<360 or -180<lon<180
        # so that interpolation is handled correctly
        # opendrit internally uses -180<lon<180 i.e. at this stage -180<lon<180
        if self.has_lon_0_360:
            lon = self.to_longitude_0_360(lon) 
        #############################################################################

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
        self.logger.debug('Reader time:\n\t\t%s (before)\n\t\t%s (after)' %
                      (time_before, time_after))
        # For variables which are not time dependent, we do not care about time
        static_variables = ['sea_floor_depth_below_sea_level', 'land_binary_mask']
        if time == time_before or all(v in static_variables for v in variables):
            time_after = None

        z = z.copy()[ind_covered]  # Send values and not reference
                                   # to avoid modifications
        if block is False or self.return_block is False:
            # Analytical reader, continous in space and time
            self.timer_end('preparing')
            env_before = self._get_variables(variables, profiles,
                                             profiles_depth,
                                             time,
                                             #time_before,
                                             reader_x, reader_y, z,
                                             block=block)
            self.logger.debug('Fetched env-before')
            self.timer_start('preparing')

        else:
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
                self.timer_start('preparing')
                self.var_block_before[blockvars_before] = \
                    ReaderBlock(reader_data_dict,
                                interpolation_horizontal=self.interpolation)
                try:
                    len_z = len(self.var_block_before[blockvars_before].z)
                except:
                    len_z = 1
                self.logger.debug(('Fetched env-block (size %ix%ix%i) ' +
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
                    self.var_block_after[blockvars_after] = \
                        ReaderBlock(
                            reader_data_dict,
                            interpolation_horizontal=self.interpolation)
                    try:
                        len_z = len(self.var_block_after[blockvars_after].z)
                    except:
                        len_z = 1

                    self.logger.debug(('Fetched env-block (size %ix%ix%i) ' +
                                  'for time after (%s)') %
                                  (len(self.var_block_after[blockvars_after].x),
                                   len(self.var_block_after[blockvars_after].y),
                                   len_z, time_after))
                    block_after = self.var_block_after[blockvars_after]

            if (block_before is not None and block_before.covers_positions(
                reader_x, reader_y) is False) or (\
                block_after is not None and block_after.covers_positions(
                    reader_x, reader_y) is False):
                self.logger.warning('Data block from %s not large enough to '
                                'cover element positions within timestep. '
                                'Buffer size (%s) must be increased.' %
                                (self.name, str(self.buffer)))

            self.timer_end('preparing')
            ############################################################
            # Interpolate before/after blocks onto particles in space
            ############################################################
            self.timer_start('interpolation')
            self.logger.debug('Interpolating before (%s) in space  (%s)' %
                          (block_before.time, self.interpolation))
            env_before, env_profiles_before = block_before.interpolate(
                    reader_x, reader_y, z, variables,
                    profiles, profiles_depth)

            if (time_after is not None) and (time_before != time):
                self.logger.debug('Interpolating after (%s) in space  (%s)' %
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
            self.logger.debug(('Interpolating before (%s, weight %.2f) and'
                           '\n\t\t      after (%s, weight %.2f) in time') %
                          (block_before.time, 1 - weight_after,
                           block_after.time, weight_after))
            env = {}
            for var in variables:
                # Weighting together, and masking invalid entries
                env[var] = np.ma.masked_invalid((env_before[var] *
                                                (1 - weight_after) +
                                                env_after[var] * weight_after))

                # if var in standard_names.keys():
                #     invalid = np.where((env[var] < standard_names[var]['valid_min'])
                #                | (env[var] > standard_names[var]['valid_max']))[0]
                #     if len(invalid) > 0:
                #         self.logger.warning('Invalid values found for ' + var)
                #         self.logger.warning(env[var][invalid])
                #         self.logger.warning('(allowed range: [%s, %s])' %
                #                         (standard_names[var]['valid_min'],
                #                          standard_names[var]['valid_max']))
                #         self.logger.warning('Replacing with NaN')
                #         env[var][invalid] = np.nan
                
            # Interpolating vertical profiles in time
            if profiles is not None:
                env_profiles = {}
                self.logger.info('Interpolating profiles in time')
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
            self.logger.debug('No time interpolation needed - right on time.')
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
                self.logger.debug('Reader SRS is the same as calculation SRS - '
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
            self.logger.debug('Masking %i elements outside coverage' %
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
