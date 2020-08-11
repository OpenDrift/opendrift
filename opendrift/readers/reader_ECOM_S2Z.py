#Horizontal reader to ECOM model at SBB
#Based on ROMS_native_reader
#Developed by Arian Dialectaquiz Santos and Danilo Silva from LHiCo - IO -USP (Brazil)



from bisect import bisect_left, bisect_right
from datetime import datetime

import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset, MFDataset, num2date
try:
    import xarray as xr
    has_xarray = True
except:
    has_xarray = False
#has_xarray = False  # Temporary disabled

from opendrift.readers.basereader import BaseReader, vector_pairs_xy
from opendrift.readers.roppy import depth_ECOM



class Reader(BaseReader):

    def __init__(self, filename=None, name=None, gridfile=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        # Map ECOM variable names to CF standard_name
        self.ECOM_variable_mapping = {
            'time': 'time',
            'sigma':'ocean_sigma_coordinate',
            'depth':'sea_floor_depth_below_sea_level',
            'ang':'angle_of_rotation_from_east_to_x',
            'elev':'sea_surface_height_above_sea_level',
            'wu': 'x_wind',
            'wv':'y_wind',
            'patm':'air_pressure_at_sea_level',
            'u': 'x_sea_water_velocity',
            'v':'y_sea_water_velocity',
            'w':'upward_sea_water_velocity',
            'salt':'sea_water_salinity',
            'temp':'sea_water_temperature',
            #Variaveis do ECOM sem CF-standard-name
            'xpos':'X_coordinate_in_Cartesian system',
            'ypos':'Y_coordinate_in_Cartesian system',
            'date':'Date_Time',
            'layer_bnds':'bounds_of_stretched_vertical_coordinate_levels',
            'x':'Corner_longitude',
            'y':'Corner_latitude',
            'h1':'delta_x',
            'h2':'delta_y',
            'FSM':'land_binary_mask', 
            'DUM':'U1_direction_mask',
            'DVM':'V1_direction_mask',
            'lon_bnds':'Vertex_longitude',
            'lat_bnds':'Vertex_latitude',
            'cbc':'Bottom Drag_Coefficient',
            'lon':'lon',
            'lat':'lat'}


             # z-levels to which sigma-layers may be interpolated (21 niveis sigma, depth max = 2000m, +0)
        self.zlevels = np.array([0, -5, -10, -15, -25,-30, -50, -75, -100, -150, -200,
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
                                list(self.ECOM_variable_mapping.keys())]
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
            raise ValueError('e')

        if 'sigma' not in self.Dataset.variables:
            dimensions = 2
        else:
            dimensions = 3

        if dimensions == 3:
            # Read sigma-coordinate values
            try:
                self.sigma = self.Dataset.variables['sigma'][:]
                self.sigma = np.nan_to_num(self.sigma)
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

            self.depth = self.Dataset.variables['depth'][:]

        if 'lat' in self.Dataset.variables:
            # Horizontal coordinates and directions
            self.lat = self.Dataset.variables['lat'][:]
            self.lat =  np.nan_to_num(self.lat)
            self.lon = self.Dataset.variables['lon'][:]            
            self.lon =  np.nan_to_num(self.lon)

           
        else:
            if gridfile is None:
                raise ValueError(filename + ' does not contain lon/lat '
                                 'arrays, please supply a grid-file '
                                 '"gridfile=<grid_file>"')
            else:
                gf = Dataset(gridfile)
                self.lat = gf.variables['lat'][:]
                self.lat =  np.nan_to_num(self.lat)
                self.lon = gf.variables['lon_'][:]
                self.lon =  np.nan_to_num(self.lon)
       


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
        self.xmin = 1
        self.delta_x = 1.

        self.ymin = 1
        self.delta_y = 1.


        if has_xarray:
            self.xmax = self.Dataset['xpos'].shape[0] - 1.
            self.ymax = self.Dataset['ypos'].shape[0] - 1.
            
            self.lon = self.lon.data  # Extract, could be avoided downstream
            self.lat = self.lat.data
            #self.sigma = self.sigma.data

       
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

        #if z is None:
            #z = np.atleast_1d(0) #Convert inputs to arrays with at least one dimension.

# Find horizontal indices corresponding to requested x and y
        if hasattr(self, 'clipped'):
            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped 
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped 


        indx_el = indx.copy()
        indy_el = indy.copy()

        
        def find_nearest(array, value):
            array = np.asarray(array)
                
            idx = (np.abs(array - value)).argmin()

            return idx

        # define an empty list as indz
        if not hasattr(self, 'z') and (z is not None):   
            print ("z ==", z)
            
            dz = self.zlevels   

            indz = [find_nearest(dz, value) for value in z]


        else:
            indz.append(0)

        indz = np.asarray(indz)

        

        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            indx = indx
            indy = indy

        print ("indx ==", indx)
        print ("indy ==", indy)
        print ("indz ==", indz)
##################################################################################################################################
#########################---This must be converted to work with sigma ----############################################
##################################################################################################################################

        # Find depth levels covering all elements 
        sigma = np.asarray(self.sigma)
        H_ND = np.asarray(self.depth)
        
        layers = self.num_layers - 1  # surface layer
        #meters_cell = depth[indy,indx]*sigma_val
        H_shape = H_ND.shape      # Save the shape of H
        H = H_ND.ravel()        # and make H 1D for easy shape maniplation
        L_C = len(sigma)
        outshape = (L_C,) + H_shape
        S = -1.0 + (0.5+np.arange(L_C))/L_C 
        A = (S - sigma)[:, None]
        B = np.outer(sigma, H)
        self.z_rho_tot = (A + B).reshape(outshape)  
        
        #print ("Z_rho_tot ==", self.z_rho_tot)
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

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

            variables[par] = np.asarray(variables[par])  # If Xarrays ######!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#################
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

                print("par_requested ==", par)
##################################################################################################################################
#########################---This must be converted to work with sigma in 3D----############################################
##################################################################################################################################
            
            if var.ndim == 4:
                # Regrid from sigma to z levels
                if len(np.atleast_1d(indz)) >= 1:
                    self.logger.debug('sigma to z for ' + varname[0])
                    if self.precalculate_s2z_coefficients is True: #what it is (view line 216)
                        y_depth = self.depth[0]
                        x_depth = self.depth[1]
                        M = len(y_depth)
                        N = len(x_depth) 
                        O = len(self.z_rho_tot)

                        if not hasattr(self, 's2z_A'):
                            self.logger.debug('Calculating sigma2z-coefficients for whole domain')
                            starttime = datetime.now()
                            dummyvar = np.ones((21, 137, 110))
                            dummy, self.s2z_total = depth_ECOM.multi_zslice(dummyvar, self.z_rho_tot, self.zlevels)
                            # Store arrays/coefficients
                            self.s2z_A = self.s2z_total[0].reshape(len(self.zlevels), 137, 110)
                            self.s2z_C = self.s2z_total[1].reshape(len(self.zlevels), 137, 110)
                            self.s2z_I = self.s2z_total[2].reshape(137, 110)
                            self.s2z_kmax = self.s2z_total[3]
                            del self.s2z_total  # Free memory
                            self.logger.info('Time: ' + str(datetime.now() - starttime))
                                        
                    self.logger.debug('Applying sigma2z-coefficients')
                    # Re-using sigma2z koefficients:
                    S = np.asarray(self.z_rho_tot)
                    Z = self.zlevels
                    kmax = Z.size
                    dummyvar = np.ones((21, 137, 110))
                    F = np.asarray(dummyvar)
                    Fshape = F.shape
                    N = F.shape[0]
                    M = F.size // N
                    I = np.arange(M, dtype= int)
                    F = F.reshape((N, M))
                    Z = Z[:, np.newaxis] + np.zeros((kmax, M))

                    S = S.reshape((N, M))
                    G = np.sum(S[np.newaxis, :, :] < Z[:, np.newaxis, :], axis=1)
                    G = np.asarray(G, dtype = int)
                    G = G.clip(1, N-1)
                    A = (Z - S[(G-1, I)])/(S[(G, I)]-S[(G-1, I)])
                    A = A.clip(0.0, 1.0)

                    R = (1-A)*F[(G-1, I)]+A*F[(G, I)]
                    R = np.nan_to_num(R)
                    interp_coefs = R.reshape((kmax,) + Fshape[1:])
                    variables[par]= var[indxTime, indz, indy, indx]*interp_coefs[indz,indy,indx]
                    variables[par] = np.nan_to_num(variables[par])

                    # Nan in input to multi_zslice gives extreme values in output
                    #variables[par][variables[par]>1e+9] = np.nan                   
                    #variables[par][variables[par]<=1e-9] = np.nan



            # If 2D array is returned due to the fancy slicing methods
            # of netcdf-python, we need to take the diagonal
            if variables[par].ndim > 1 and block is False:
                variables[par] = variables[par].diagonal()

            # Mask values outside domain
            variables[par] = np.ma.array(variables[par], ndmin=2, mask=False)
            if block is False:
                variables[par].mask[outside] = True
        
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

        if block is True:

            variables['x'] = indx
            variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['x'] = variables['x'].astype(np.float)
        variables['y'] = variables['y'].astype(np.float)
        variables['time'] = nearestTime

        if 'x_sea_water_velocity' or 'x_wind' in variables.keys():
 

            if not hasattr(self, 'angle_of_rotation_from_east_to_x'):
                self.logger.debug('Reading angle between xi and east...')
                self.angle_of_rotation_from_east_to_x = self.Dataset.variables['ang'][:]
            if has_xarray is False:
                rad = self.angle_of_rotation_from_east_to_x[tuple(np.meshgrid(indy, indx))].T
            else:
                rad = self.angle_of_rotation_from_east_to_x[indy, indx]
                rad = np.nan_to_num(rad)

            #print ("angle_of_rotation_from_east_to_x ==", rad)

            if 'x_sea_water_velocity' in variables.keys():
                #variables['x_sea_water_velocity'], \
                #    variables['y_sea_water_velocity'] = rotate_vectors_angle(
                #        variables['x_sea_water_velocity'],
                #        variables['y_sea_water_velocity'], rad)
                
                variables['x_sea_water_velocity'] = np.nan_to_num(variables['x_sea_water_velocity'])
                variables['y_sea_water_velocity'] = np.nan_to_num(variables['y_sea_water_velocity'])
              
                print ("u ==", variables['x_sea_water_velocity'])
                print ("v ==", variables['y_sea_water_velocity'])
                print ("w ==", variables['upward_sea_water_velocity'])

            if 'x_wind' in variables.keys():

               # variables['x_wind'], \
               #     variables['y_wind'] = rotate_vectors_angle(
               #         variables['x_wind'],
               #         variables['y_wind'], rad)
               # 
                variables['x_wind'] = np.nan_to_num(variables['x_wind'])
                variables['y_wind'] = np.nan_to_num(variables['y_wind'])



                print ("wu ==", variables['x_wind'])
                print ("wv ==", variables['y_wind'])


        # Masking NaN of the others variables, considering u and v always requested
        for var in requested_variables:
          
            variables[var] = np.nan_to_num(variables[var])

        self.logger.debug('Time for ECOM reader: ' + str(datetime.now()-start_time))

        return variables
        
#This follow function is not to be used at the SBB grid.

def rotate_vectors_angle(x_sea_water_velocity, y_sea_water_velocity, radians):
    x_sea_water_velocity2 = x_sea_water_velocity*np.cos(radians) - y_sea_water_velocity*np.sin(radians)
    y_sea_water_velocity2 = x_sea_water_velocity*np.sin(radians) + y_sea_water_velocity*np.cos(radians)
    return x_sea_water_velocity2, y_sea_water_velocity2
