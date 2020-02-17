#-*-coding:utf-8-*-
from bisect import bisect_left, bisect_right
from datetime import datetime

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
try:
    import xarray as xr
    has_xarray = True
except:
    has_xarray = False
#has_xarray = False  # Temporary disabled

from opendrift.readers.basereader import BaseReader, vector_pairs_xy
from opendrift.readers.roppy import depth


class Reader(BaseReader):

    def __init__(self, filename=None, name=None, gridfile=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        # Map ECOM variable names to CF standard_name
        self.ECOM_variable_mapping = {
            'time': {'time': 'time'},
            'sigma':'ocean_sigma_coordinate',
            'depth':'depth',
            'ang':'angle_of_rotation_from_east_to_x',
            'elev':'sea_surface_height_above_sea_level',
            'x_wind'{'wu': 'x_wind'},
            'y_wind'{'wv':'y_wind'},
            'patm':'air_pressure_at_sea_level',
            'x_sea_water_velocity'{'u': 'x_sea_water_velocity'},
            'y_sea_water_velocity'{'v':'y_sea_water_velocity'},
            'w':'upward_sea_water_velocity',
            'salt':'sea_water_salinity',
            'temp':'sea_water_temperature',
            #Variaveis do ECOM sem CF-standard-name
            'xpos':'X_coordinate_in_Cartesian system',
            'ypos':'Y_coordinate_in_Cartesian system',
            'zpos':'Z_coordinate_in_Cartesian_system',
            'date':'Date_Time',
            'layer_bnds':'bounds_of_stretched_vertical_coordinate_levels',
            'x':'Corner_longitude',
            'y':'Corner_latitude',
            'h1':'dx_metric',
            'h2':'dy_metric',
            'FSM':'land_binary_mask', #0 ou 1 essa tal Free Surface Mask
            'DUM':'U1_direction_mask',
            'DVM':'V1_direction_mask',
            'lon_bnds':'Vertex_longitude',
            'lat_bnds':'Vertex_latitude',
            'cbc':'Bottom Drag_Coefficient',
            'lon':'longitude',
            'lat':'latitude'}

             # z-levels to which sigma-layers may be interpolated (21 niveis sigma, depth max = 2000m)
        self.zlevels = np.array([0, -5, -10, -25,-30, -50, -75, -100, -150, -200,
            -250, -300, -400, -500, -600, -700, -800, -900, -1000, -1500,
            -2000])



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
                                list(self.ECOM_variable_mapping.keys())  +
                                ['time', 'sigma', 'depth', 'elev', 'ang']]

                    self.logger.debug('Dropping variables: %s' % dropvars)
                    ds = ds.drop(dropvars)
                    return ds

                if has_xarray is True:
                    self.Dataset = xr.open_mfdataset(filename,
                        chunks={'time': 1}, concat_dim='time',
                        preprocess=drop_non_essential_vars_pop,
                        data_vars='minimal', coords='minimal')
                else:
                    self.Dataset = MFDataset(filename)
            else:
                self.logger.info('Opening file with Dataset')
                if has_xarray is True:
                    self.Dataset = xr.open_dataset(filename)
                else:
                    self.Dataset = Dataset(filename, 'r')
        except Exception as e:
            raise ValueError(e)

        if 'sigma' not in self.Dataset.variables:
            dimensions = 2
        else:
            dimensions = 3

        if dimensions == 3:
            # Read sigma-coordinate values
            try:
                self.sigma = self.Dataset.variables['sigma'][:]
            except:
                num_sigma = len(self.Dataset.dimensions['sigma'])
                self.logger.warning(
                    'sigma not available in dataset, constructing from'
                    ' number of layers (%s).' % num_sigma)
                self.sigma = (np.arange(num_sigma)+.5-num_sigma)/num_sigma


                # Read sigma-coordinate transform parameters
            try:
                self.Dataset.variables['sigma'].set_auto_mask(False)
            except:
                pass
            self.sigma = self.Dataset.variables['sigma'][:]
            try:
                self.depth = self.Dataset.variables['depth'][:]
            except:
                if has_xarray is True:
                    self.depth = self.Dataset.variables['depth'].data  # scalar
                else:
                    self.depth = self.Dataset.variables['depth'][0]

            self.num_layers = len(self.sigma)

        else:
            self.num_layers = 1
            self.ECOM_variable_mapping['u'] = 'x_sea_water_velocity'
            self.ECOM_variable_mapping['v'] = 'y_sea_water_velocity'


        if 'lat' in self.Dataset.variables:
            # Horizontal coordinates and directions
            self.lat = self.Dataset.variables['lat'][:]
            self.lon = self.Dataset.variables['lon'][:]
           #self.lat[np.where(self.lat == 0)] = np.nan
           #self.lon[np.where(self.lon == 0)] = np.nan
           #self.lat = self.lat[~np.isnan(self.lat)]
           #self.lon = self.lon[~np.isnan(self.lon)]
        else:
            if gridfile is None:
                raise ValueError(filename + ' does not contain lon/lat '
                                 'arrays, please supply a grid-file '
                                 '"gridfile=<grid_file>"')
            else:
                gf = Dataset(gridfile)
                self.lat = gf.variables['lat'][:]
                self.lon = gf.variables['lon_'][:]

        # Get time coverage

        ocean_time = self.Dataset.variables['time']
        if has_xarray:
            self.times = [datetime.utcfromtimestamp((OT -
                          np.datetime64('1970-01-01T00:00:00Z')
                            ) / np.timedelta64(1, 's'))
                          for OT in ocean_time.data]
        else:
            time_units = ocean_time.__dict__['units']
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
        self.h1 = 1.
        self.ymin = 0.
        self.h2 = 1.
        if has_xarray:
            self.xmax = self.Dataset['xpos'].shape[0] - 1.
            self.ymax = self.Dataset['ypos'].shape[0] - 1.
            self.lon = self.lon.data  # Extract, could be avoided downstream
            self.lat = self.lat.data
            self.sigma = self.sigma.data
        else:
            self.xmax = np.float(len(self.Dataset.dimensions['xpos'])) - 1
            self.ymax = np.float(len(self.Dataset.dimensions['ypos'])) - 1

        self.name = 'ECOM'

        self.precalculate_s2z_coefficients = True

        # Find all variables having standard_name
        self.variables = []
        for var_name in self.Dataset.variables:
            if var_name in self.ECOM_variable_mapping.keys():
                var = self.Dataset.variables[var_name]
                self.variables.append(self.ECOM_variable_mapping[var_name])

        # Run constructor of parent Reader class
        super(Reader, self).__init__()
        '''[Super is used to] return a proxy object that delegates method calls to a parent or sibling class of type.
        This is useful for accessing inherited methods that have been overridden in a class. The search order is same as that used by getattr() except that the type itself is skipped'''

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
            z = np.atleast_1d(0) #Convert inputs to arrays with at least one dimension.


# Find horizontal indices corresponding to requested x and y
        if hasattr(self, 'clipped'): #The hasattr() method returns:True, if object has the given named attribute ; False, if object has no given named attribute

            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.h1).astype(int) + clipped #np.floor organiza em ordem crescente os dados
        indy = np.floor((y-self.ymin)/self.h2).astype(int) + clipped
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

# Find depth levels covering all elements (Qual o análogo de HC - profundidade critica do ROMS no ECOM?)

        if z.min() == 0 or not hasattr(self, 'depth'):
            indz = self.num_layers - 1  # surface layer
            variables['z'] = 0

        else:
            # Find the range of indices covering given z-values
            if not hasattr(self, 'depth'):
                self.logger.debug('Reading sea floor depth...')
                self.depth = \
                    self.Dataset.variables['depth'][:]

                Htot = self.depth
                self.z__tot = depth.depth(Htot, self.depth, self.sigma)

            if has_xarray is False:
                indxgrid, indygrid = np.meshgrid(indx, indy)
                H = self.depth[indygrid, indxgrid]
            else:
                H = self.depth[indy, indx]
            z_rho = depth.depth(H, self.depth, self.sigma)
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
                       self.ECOM_variable_mapping.items() if cf == par]
            var = self.Dataset.variables[varname[0]]

            if par == 'land_binary_mask':
                if not hasattr(self, 'land_binary_mask'):
                    # Read landmask for whole domain, for later re-use
                    self.land_binary_mask = \
                        1 - self.Dataset.variables['FSM'][:]
                if has_xarray is False:
                    indxgrid, indygrid = np.meshgrid(indx, indy)
                    variables[par] = self.land_binary_mask[indygrid, indxgrid]
                else:
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
                if has_xarray is False:
                    indxgrid, indygrid = np.meshgrid(indx, indy)
                else:
                    indxgrid = indx
                    indygrid = indy
                if par == 'x_sea_water_velocity':
                    if not hasattr(self, 'DUM'):
                        if 'DUM' in self.Dataset.variables:
                            self.mask_u = self.Dataset.variables['DUM'][:]
                        else:
                            self.mask_u = self.Dataset.variables['FSM'][:]
                    mask = self.mask_u[indygrid, indxgrid]
                elif par == 'y_sea_water_velocity':
                    if not hasattr(self, 'DVM'):
                        if 'DVM' in self.Dataset.variables:
                            self.mask_v = self.Dataset.variables['DVM'][:]
                        else:
                            self.mask_v = self.Dataset.variables['FSM'][:]
                    mask = self.mask_v[indygrid, indxgrid]
                else:
                    if not hasattr(self, 'FSM'):
                        self.mask_rho = self.Dataset.variables['FSM'][:]
                    mask = self.mask_rho[indygrid, indxgrid]
                if has_xarray is True:
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
                        M = depth.shape[0]
                        N = depth.shape[1]
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

        if 'x_sea_water_velocity' or 'x_wind' in variables.keys():
            # We must rotate current vectors
            if not hasattr(self, 'angle_of_rotation_from_east_to_x'):
                self.logger.debug('Reading angle between xi and east...')
                self.angle_of_rotation_from_east_to_x = self.Dataset.variables['ang'][:]
            if has_xarray is False:
                rad = self.angle_of_rotation_from_east_to_x[tuple(np.meshgrid(indy, indx))].T
            else:
                rad = self.angle_of_rotation_from_east_to_x[indy, indx]
            if 'x_sea_water_velocity' in variables.keys():
                variables['x_sea_water_velocity'], \
                    variables['y_sea_water_velocity'] = rotate_vectors_angle(
                        variables['x_sea_water_velocity'],
                        variables['y_sea_water_velocity'], rad)

            if 'x_wind' in variables.keys():
                variables['x_wind'], \
                    variables['y_wind'] = rotate_vectors_angle(
                        variables['x_wind'],
                        variables['y_wind'], rad)

        # Masking NaN
        for var in requested_variables:
            variables[var] = np.ma.masked_invalid(variables[var])

        self.logger.debug('Time for ECOM reader:' + str(datetime.now()-start_time))

        return variables


def rotate_vectors_angle(u, v, radians):
    u2 = u*np.cos(radians) - v*np.sin(radians)
    v2 = u*np.sin(radians) + v*np.cos(radians)
    return u2, v2
