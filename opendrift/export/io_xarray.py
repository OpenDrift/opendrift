import sys
import os
from datetime import datetime, timedelta
import logging; logging.captureWarnings(True); logger = logging.getLogger(__name__)
import string
from shutil import move

import numpy as np
#from netCDF4 import Dataset, num2date, date2num
import xarray as xr
import dask.array

# Module with functions to export/import trajectory data to/from netCDF file
# Strives to be compliant with netCDF CF-convention on trajectories
# https://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#idp8377728
# https://geo-ide.noaa.gov/wiki/index.php?title=NODC_NetCDF_Trajectory_Template
# http://xarray.pydata.org/en/stable/user-guide/data-structures.html#creating-a-dataset

skip_parameters = ['ID']  # Do not write to file


def init(self, filename, method="zarr"):
    """
    Assemble the data into a xarray structure so that it can be written in any
    format CF compliant (netcdf, zarr,...)
    Args:
        :param filename: A single file (or a pattern of files?).
        :type filename: string, required.
        :param method: A method to write the file
        :type filename: string, optional (netcdf4, zarr), grib, h5netcdf?)

    """
    assert method in ["netcdf4", "scipy", 'h5netcdf', 'zarr'],\
                    "method not recognised"
    assert method == 'zarr', '{} engine not yet supported'.format(method)
    self.outfile_name = filename
    self.xr_dataset = xr.Dataset(
                coords={
                    'trajectory':
                        (['trajectory'], np.arange(self.num_elements_total(),
                        dtype=np.int32),
                        {'cf_role':'trajectory_id', 'units':'1'}),
                    'time':
                        (['time'], self.times),
                            {
                            'units':'seconds since 1970-01-01 00:00:00',
                            'standard_name': 'time',
                            'long_name': 'time',
                            # 'axis': 'T'?
                            }
                        ),
                    'lon': (['trajectory', 'time']),
                    'lat': (['trajectory', 'time'])
                    },
                attrs={
                    'Conventions': 'CF-1.6',
                    'standard_name_vocabulary': 'CF-1.6',
                    'featureType': 'trajectory',
                    'history': 'Created ' + str(datetime.now()),
                    'source': 'Output from simulation with OpenDrift',
                    'model_url': 'https://github.com/OpenDrift/opendrift',
                    'opendrift_class': self.__class__.__name__,
                    'opendrift_module': self.__class__.__module__,
                    'readers': str(self.readers.keys()),
                    'time_coverage_start': str(self.start_time),
                    'time_step_calculation': str(self.time_step),
                    'time_step_output': str(self.time_step_output),
                    },)

    # Write config settings
    for key in self._config:
        value = self.get_config(key)
        if isinstance(value, (bool, type(None))):
            value = str(value)
        self.outfile.attrs['config_' + key]= value

    # Write additionaly metadata attributes, if given
    if hasattr(self, 'metadata_dict'):
        for key, value in self.metadata_dict.items():
            self.xr_dataset.attrs[key]= str(value)

    # Add all element properties as variables
    for prop in self.history.dtype.fields:
        if prop in skip_parameters:
            continue
        # Note: Should use 'f8' if 'f4' is not accurate enough,
        #       at expense of larger files
        try:
            dtype = self.history.dtype[prop]
        except:
            dtype = 'f4'
        self.xr_dataset[prop]=(
            ('trajectory', 'time'), # coord
            dask.array.zeros(  # empty dask array
                (self.num_elements_total(),len(self.times)),
                dtype= dtype,
                chunks={'time':1}),
            {'coordinates':'lat lon time'})
        for subprop in self.history_metadata[prop].items():
            if subprop[0] not in ['dtype', 'constant', 'default', 'seed']:
                # Apparently axis attribute shall not be given for lon and lat:
                if prop in ['lon', 'lat'] and subprop[0] == 'axis':
                    continue
                self.xr_dataset[prop].attrs[subprop[0]]= subprop[1]

    # writing all the metadata without computing any array value
    self.xr_dataset.to_zarr(self.outfile_name, , compute=False)

def write_buffer(self):
    """
    Write the buffer in the dataset and write it in the right part of the
    dask.array.
    """

    self.buffer_ds = xr.Dataset()
    num_steps_to_export = self.steps_output - self.steps_exported
    slice = slice(self.steps_exported, self.steps_output+1)
    for prop in self.history_metadata:
        if prop in skip_parameters:
            continue
        self.xr_dataset[prop].isel(time=slice) = \
            self.history[prop][:, :num_steps_to_export]
    self.xr_dataset.isel(time=slice).to_zarr(self.outfile_name, region={
            "time": slice})
    logger.info('Wrote %s steps to file %s' % (num_steps_to_export,
                                                self.outfile_name))
    self.history.mask = True  # Reset history array, for new data
    self.steps_exported = self.steps_exported + num_steps_to_export
    self.xr_dataset.attrs['steps_exported'] = self.steps_exported



def close(self):
    """
    Writing final metadata of the dataset before closing the dataset.
    """

    # Write status categories metadata
    status_dtype = self.ElementType.variables['status']['dtype']
    self.xr_dataset['status'].attrs['valid_range'] = np.array(
        (0, len(self.status_categories) - 1)).astype(status_dtype)
    self.xr_dataset['status'].attrs['flag_values'] = \
        np.array(np.arange(len(self.status_categories)), dtype=status_dtype)
    self.xr_dataset['status'].attrs['flag_meanings'] = \
        " ".join(self.status_categories)

    # Write origin_marker definitions
    if 'origin_marker' in self.xr_dataset.data_vars:
        self.xr_dataset['origin_marker'].attrs['flag_values'] = \
            np.array(np.arange(len(self.origin_marker)))
        self.xr_dataset['origin_marker'].attrs['flag_meanings'] = \
            " ".join(self.origin_marker.values())

    # Write final timesteps to file
    self.xr_dataset.attrs['time_coverage_end'] = str(self.time)

    # Write performance data
    self.xr_dataset.attrs['performance'] = self.performance()

    # Write metadata items anew, if any are added during simulation
    # this should be optimised with flags?
    if hasattr(self, 'metadata_dict'):
        for key, value in self.metadata_dict.items():
            self.xr_dataset.attrs[key]= value

    # Write min and max values as variable attributes
    for var in self.history_metadata:
        if var in self.xr_dataset.data_vars:
            self.xr_dataset[var].attrs['minval']= self.minvals[var])
            self.xr_dataset[var].attrs['maxval']= self.maxvals[var])

    # Write bounds metadata
    self.xr_dataset.attrs['geospatial_lat_min'] = self.history['lat'].min()
    self.xr_dataset.attrs['geospatial_lat_max'] = self.history['lat'].max()
    self.xr_dataset.attrs['geospatial_lat_units'] = 'degrees_north'
    self.xr_dataset.attrs['geospatial_lat_resolution'] = 'point'
    self.xr_dataset.attrs['geospatial_lon_min'] = self.history['lon'].min()
    self.xr_dataset.attrs['geospatial_lon_max'] = self.history['lon'].max()
    self.xr_dataset.attrs['geospatial_lon_units'] = 'degrees_east'
    self.xr_dataset.attrs['geospatial_lon_resolution'] = 'point'
    self.xr_dataset.attrs['runtime'] = str(datetime.now() -
                               self.timers['total time'])
    # Finally close file
    self.xr_dataset.close()

# def import_file_xarray(self, filename, chunks):
#
#     import xarray as xr
#     logger.debug('Importing with Xarray from ' + filename)
#     self.ds = xr.open_dataset(filename, chunks=chunks)
#
#     self.steps_output = len(self.ds.time)
#     ts0 = (self.ds.time[0] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
#     self.start_time = datetime.utcfromtimestamp(float(ts0))
#     tse = (self.ds.time[-1] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
#     self.end_time = datetime.utcfromtimestamp(float(tse))
#     if len(self.ds.time) > 1:
#         ts1 = (self.ds.time[1] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
#         self.time_step_output = timedelta(seconds=float(ts1 - ts0))
#     self.time = self.end_time  # Using end time as default
#     self.status_categories = self.ds.status.flag_meanings.split()
#     if 'flag_meanings' in self.ds.origin_marker.attrs:
#         self.origin_marker = [s.replace('_', ' ') for s in self.ds.origin_marker.flag_meanings.split()]
#
#     num_elements = len(self.ds.trajectory)
#     elements=np.arange(num_elements)
#
#     # Data types for variables
#     dtype = np.dtype([(var[0], var[1]['dtype'])
#                       for var in self.ElementType.variables.items()])
#     history_dtype_fields = [
#         (name, self.ElementType.variables[name]['dtype'])
#         for name in self.ElementType.variables]
#     # Add environment variables
#     self.history_metadata = self.ElementType.variables.copy()
#     for env_var in self.required_variables:
#         if env_var in self.ds.variables:
#             history_dtype_fields.append((env_var, np.dtype('float32')))
#             self.history_metadata[env_var] = {}
#     history_dtype = np.dtype(history_dtype_fields)
#
#     # Masking where elements are not been seeded
#     # TODO: mask from deactivation towards end
#     for da in ['lon', 'lat']:
#         self.ds[da] = self.ds[da].where(self.ds.status>=0)
#
#     if 'minval' in self.ds.lon.attrs:
#         self.lonmin = np.float32(self.ds.lon.minval)
#         self.latmin = np.float32(self.ds.lat.minval)
#         self.lonmax = np.float32(self.ds.lon.maxval)
#         self.latmax = np.float32(self.ds.lat.maxval)

# def import_file(self, filename, times=None, elements=None, load_history=True):
#     """Create OpenDrift object from imported file.
#      times: indices of time steps to be imported, must be contineous range.
#      elements: indices of elements to be imported
#     """
#
#     logger.debug('Importing from ' + filename)
#     infile = Dataset(filename, 'r')
#     # 'times' can be used to import subset. Not yet implemented.
#     if times is None and hasattr(infile, 'steps_exported'):
#         self.steps_output = infile.steps_exported
#         times = np.arange(infile.steps_exported)
#     else:
#         #self.steps_output = len(infile.dimensions['time'])
#         self.steps_output = len(times)
#
#     filetime = infile.variables['time'][times]
#     units = infile.variables['time'].units
#     self.start_time = num2date(filetime[0], units)
#     if len(filetime) > 1:
#         self.end_time = num2date(filetime[self.steps_output-1], units)  # Why -1?
#         self.time_step_output = num2date(filetime[1], units) - self.start_time
#     else:
#         self.time_step_output = timedelta(hours=1)
#         self.end_time = self.start_time
#     self.time = self.end_time  # Using end time as default
#     self.status_categories = infile.variables['status'].flag_meanings.split()
#
#     if elements is None:
#         num_elements = len(infile.dimensions['trajectory'])
#         elements=np.arange(num_elements)
#     else:
#         num_elements=len(elements)
#
#     dtype = np.dtype([(var[0], var[1]['dtype'])
#                       for var in self.ElementType.variables.items()])
#
#     history_dtype_fields = [
#         (name, self.ElementType.variables[name]['dtype'])
#         for name in self.ElementType.variables]
#     # Add environment variables
#     self.history_metadata = self.ElementType.variables.copy()
#     for env_var in self.required_variables:
#         if env_var in infile.variables:
#             history_dtype_fields.append((env_var, np.dtype('float32')))
#             self.history_metadata[env_var] = {}
#     history_dtype = np.dtype(history_dtype_fields)
#
#     # Import dataset (history)
#     if load_history is True:
#         self.history = np.ma.array(
#             np.zeros([num_elements, self.steps_output]),
#             dtype=history_dtype, mask=[True])
#         for var in infile.variables:
#             if var in ['time', 'trajectory']:
#                 continue
#             try:
#                 self.history[var] = infile.variables[var][elements, times]
#             except Exception as e:
#                 logger.info(e)
#                 pass
#
#         # Initialise elements from given (or last) state/time
#         firstlast = np.ma.notmasked_edges(self.history['status'], axis=1)
#         index_of_last = firstlast[1][1]
#         kwargs = {}
#         for var in infile.variables:
#             if var in self.ElementType.variables:
#                 kwargs[var] = self.history[var][
#                     np.arange(len(index_of_last)), index_of_last]
#         kwargs['ID'] = np.arange(len(elements)) + 1
#         self.elements = self.ElementType(**kwargs)
#         self.elements_deactivated = self.ElementType()
#     else:
#         self.history = None
#         logger.warning('Not importing history')
#
#     # Remove elements which are scheduled for deactivation
#     self.remove_deactivated_elements()
#
#     # Import and apply config settings
#     attributes = infile.ncattrs()
#     for attr in attributes:
#         if attr.startswith('config_'):
#             value = infile.getncattr(attr)
#             conf_key = attr[7:]
#             if value == 'True':
#                 value = True
#             if value == 'False':
#                 value = False
#             if value == 'None':
#                 value = None
#             try:
#                 self.set_config(conf_key, value)
#                 logger.debug('Setting imported config: %s -> %s' %
#                              (conf_key, value))
#             except:
#                 logger.warning('Could not set config: %s -> %s' %
#                                 (conf_key, value))
#
#     # Import time steps from metadata
#     def timedelta_from_string(timestring):
#         if 'day' in timestring:
#             days = int(timestring.split('day')[0])
#             hs = timestring.split(' ')[-1]
#             th = datetime.strptime(hs, '%H:%M:%S')
#             return timedelta(days=days, hours=th.hour, minutes=th.minute, seconds=th.second)
#         else:
#             t = datetime.strptime(timestring, '%H:%M:%S')
#             return timedelta(
#                 hours=t.hour, minutes=t.minute, seconds=t.second)
#     try:
#         self.time_step = timedelta_from_string(infile.time_step_calculation)
#         self.time_step_output = timedelta_from_string(infile.time_step_output)
#     except Exception as e:
#         logger.warning(e)
#         logger.warning('Could not parse time_steps from netCDF file')
#
#     infile.close()
