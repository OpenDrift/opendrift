#Reader to wave model:
#Building to work with UMWM


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

class Reader(BaseReader):

    def __init__(self, filename=None, name=None, gridfile=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        # Map waves variable names to CF standard_name
        self.waves_variable_mapping = {
            'time': 'time',
            'x':'X_coordinate_in_Cartesian system',
            'y':'Y_coordinate_in_Cartesian system',
            'lon':'lon',
            'lat':'lat',
            'mwp' :'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
            'dwp': 'sea_surface_wave_period_at_variance_spectral_density_maximum',
            'u_stokes' : 'sea_surface_wave_stokes_drift_x_velocity',
            'v_stokes': 'sea_surface_wave_stokes_drift_y_velocity'}


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
                                list(self.waves_variable_mapping.keys())]
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
            self.xmax = self.Dataset['x'].shape[0] - 1.
            self.ymax = self.Dataset['y'].shape[0] - 1.
            
            self.lon = self.lon.data  # Extract, could be avoided downstream
            self.lat = self.lat.data

        else:
            self.xmax = np.float(len(self.Dataset.dimensions['x'])) - 1
            self.ymax = np.float(len(self.Dataset.dimensions['y'])) - 1

        self.name = 'waves'

        # Find all variables having standard_name
        self.variables = []
        for var_name in self.Dataset.variables:
            if var_name in self.waves_variable_mapping.keys():
                var = self.Dataset.variables[var_name]
                self.variables.append(self.waves_variable_mapping[var_name])
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

        if z is None:
            z = np.atleast_1d(0)

# Find horizontal indices corresponding to requested x and y
        if hasattr(self, 'clipped'):
            clipped = self.clipped
        else: clipped = 0
        indx = np.floor((x-self.xmin)/self.delta_x).astype(int) + clipped 
        indy = np.floor((y-self.ymin)/self.delta_y).astype(int) + clipped 


        indx_el = indx.copy()
        indy_el = indy.copy()

        if block is True:
            # Adding buffer, to cover also future positions of elements
            buffer = self.buffer
            indx = indx
            indy = indy

        print ("indx ==", indx)
        print ("indy ==", indy)

        #Working with the variables:
        #mask_values = {}
        for par in requested_variables:
            varname = [name for name, cf in
                       self.waves_variable_mapping.items() if cf == par]
            var = self.Dataset.variables[varname[0]]

            if var.ndim == 3:
                variables[par] = var[indxTime, indy, indx]
            
            else:
                raise Exception('Wrong dimension of variable: ' +
                                self.variable_mapping[par])

            variables[par] = np.asarray(variables[par])  # If Xarrays ######!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#################
            start = datetime.now()

                            
        if block is True:

            variables['x'] = indx
            variables['y'] = indy
        else:
            variables['x'] = self.xmin + (indx-1)*self.delta_x
            variables['y'] = self.ymin + (indy-1)*self.delta_y

        variables['x'] = variables['x'].astype(np.float)
        variables['y'] = variables['y'].astype(np.float)
        variables['time'] = nearestTime



        if 'sea_surface_wave_stokes_drift_x_velocity' in variables.keys():
                                
                variables['sea_surface_wave_stokes_drift_x_velocity'] = np.nan_to_num(variables['sea_surface_wave_stokes_drift_x_velocity'])
                variables['sea_surface_wave_stokes_drift_y_velocity'] = np.nan_to_num(variables['sea_surface_wave_stokes_drift_y_velocity'])
                variables['sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'] = np.nan_to_num(variables['sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'])
                variables['sea_surface_wave_period_at_variance_spectral_density_maximum'] = np.nan_to_num(variables['sea_surface_wave_period_at_variance_spectral_density_maximum'])
              
                print ("u_stokes ==", variables['sea_surface_wave_stokes_drift_x_velocity'])
                print ("v_stokes ==", variables['sea_surface_wave_stokes_drift_y_velocity'])

                print ("mwp ==", variables['sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'])
                print ("dwp ==", variables['sea_surface_wave_period_at_variance_spectral_density_maximum'])

        # Masking NaN of the others variables, considering u and v always requested
        for var in requested_variables:
          
            variables[var] = np.nan_to_num(variables[var])

        self.logger.debug('Time for Waves reader: ' + str(datetime.now()-start_time))

        return variables

