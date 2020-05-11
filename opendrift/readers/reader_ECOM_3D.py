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
from opendrift.readers.roppy import depth


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
            'xpos':'X_coordinate_in_Cartesian_system',
            'ypos':'Y_coordinate_in_Cartesian_system',
            'z':'Z_coordinate_in_Cartesian_system',
            'date':'Date_Time',
            'layer_bnds':'bounds_of_stretched_vertical_coordinate_levels',
            'x':'Corner_longitude',
            'y':'Corner_latitude',
            'h1':'delta_x',
            'h2':'delta_y',
            'FSM':'land_binary_mask', #0 ou 1 essa tal Free Surface Mask
            'DUM':'U1_direction_mask',
            'DVM':'V1_direction_mask',
            'lon_bnds':'Vertex_longitude',
            'lat_bnds':'Vertex_latitude',
            'cbc':'Bottom Drag_Coefficient',
            'lon':'lon',
            'lat':'lat'}


            
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



        if 'z' not in self.Dataset.variables:
            dimensions = 2
        else:
            dimensions = 3

        if dimensions == 3:
            
            try:
                self.z = self.Dataset.variables['z'][:]
            except:
                num_z = len(self.Dataset.dimensions['z'])
                self.logger.warning(
                    'Z not available in dataset, constructing from'
                    ' number of layers (%s).' % num_z)
                self.z = (np.arange(num_z)+.5-num_z)/num_z
                
                # Read z-coordinate transform parameters
            try:
                self.Dataset.variables['z'].set_auto_mask(False)
            except:
                pass
            self.z = self.Dataset.variables['z'][:]
            try:
                self.depth = self.Dataset.variables['depth'][:]

            except:
                if has_xarray is True:
                    self.depth = self.Dataset.variables['depth'].data  # scalar
                else:
                    self.depth = self.Dataset.variables['depth'][0]

            self.num_layers = len(self.z)


        else:
            self.num_layers = 1
            self.ECOM_variable_mapping['u'] = 'x_sea_water_velocity'
            self.ECOM_variable_mapping['v'] = 'y_sea_water_velocity' 

            self.depth = self.Dataset.variables['depth'][:]


        
        if 'lat' in self.Dataset.variables:
            # Horizontal coordinates and directions
                

            #self.lat = np.nan_to_num(self.Dataset.variables['lat'][:])
            #self.lon = np.nan_to_num(self.Dataset.variables['lon'][:])

            self.lat = self.Dataset.variables['lat'][:]
            self.lon = self.Dataset.variables['lon'][:]

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
        self.delay_time = self.times[1]
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

       
        else:
            self.xmax = np.float(len(self.Dataset.dimensions['xpos'])) - 1
            self.ymax = np.float(len(self.Dataset.dimensions['ypos'])) - 1
          

        self.name = 'ECOM'

       # self.precalculate_s2z_coefficients = True

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
        if hasattr(self, 'z') and (z is not None):   
            print ("z ==", z)
            
            dz = -np.array(self.Dataset['z'])   

            indz = [find_nearest(dz, value) for value in z]


        else:
            indz.append(0)

        indz = np.asarray(indz)

        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            
            indx = indx
            indy = indy
            indz = indz
                

        print ("indz = ", indz)
        print ("indx ==", indx)
        print ("indy ==", indy)


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

            #    variables[par] = np.asarray(variables[par])  # If Xarray
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
                    
        if block is True:
            
            variables['x'] = indx
            variables['y'] = indy
            variables['z'] = indz
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['x'] = variables['x'].astype(np.float)
        variables['y'] = variables['y'].astype(np.float)
        variables['z'] = variables['z'].astype(np.float)

        variables['time'] = nearestTime

        if 'x_sea_water_velocity' or 'x_wind' in variables.keys():
 

            if not hasattr(self, 'angle_of_rotation_from_east_to_x'):
                self.logger.debug('Reading angle between xi and east...')
                self.angle_of_rotation_from_east_to_x = self.Dataset.variables['ang'][:]
            if has_xarray is False:
                rad = self.angle_of_rotation_from_east_to_x[tuple(np.meshgrid(indy, indx))].T
            else:
                rad = self.angle_of_rotation_from_east_to_x[indy, indx]

            print ("angle_of_rotation_from_east_to_x ==", rad)

            if 'x_sea_water_velocity' in variables.keys():
       
                variables['x_sea_water_velocity'], \
                    variables['y_sea_water_velocity'] = rotate_vectors_angle(
                        variables['x_sea_water_velocity'],
                        variables['y_sea_water_velocity'], rad)

                print ("u ==", variables['x_sea_water_velocity'])
                print ("v ==", variables['y_sea_water_velocity'])
                print ("w ==", variables['upward_sea_water_velocity'])

            if 'x_wind' in variables.keys():
                variables['x_wind'], \
                    variables['y_wind'] = rotate_vectors_angle(
                        variables['x_wind'],
                        variables['y_wind'], rad)

                print ("wu ==", variables['x_wind'])
                print ("wv ==", variables['y_wind'])


        # Masking NaN
        for var in requested_variables:
          
            #variables[var] = np.ma.masked_invalid(variables[var])
            variables[var] = np.nan_to_num(variables[var])
        self.logger.debug('Time for ECOM reader: ' + str(datetime.now()-start_time))

        return variables

def rotate_vectors_angle(x_sea_water_velocity, y_sea_water_velocity, radians):
    x_sea_water_velocity2 = x_sea_water_velocity*np.cos(radians) - y_sea_water_velocity*np.sin(radians)
    y_sea_water_velocity2 = x_sea_water_velocity*np.sin(radians) + y_sea_water_velocity*np.cos(radians)
    return x_sea_water_velocity2, y_sea_water_velocity2
