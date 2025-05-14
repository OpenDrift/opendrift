import sys
import os
from datetime import datetime, timedelta
import logging; logging.captureWarnings(True); logger = logging.getLogger(__name__)
import string
import shutil
import pandas as pd

import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date, date2num
from opendrift.models.basemodel import Mode

# Module with functions to export/import trajectory data to/from netCDF file.
# Follows netCDF CF-convention on trajectories:
# https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#_multidimensional_array_representation_of_trajectories

def datetime_from_datetime64(dt64):
    t = (dt64 - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    return datetime.utcfromtimestamp(float(t))

def init(self, filename):

    logger.debug('Initialising output netCDF file {self.outfile_name}')
    self.outfile_name = filename
    if os.path.exists(self.outfile_name):
        logger.warning(f'Deleting existing {self.outfile_name}')
        os.remove(self.outfile_name)

def write_buffer(self):

    if not os.path.exists(self.outfile_name):
        logger.debug('Initialising output netCDF file '
                    f'{self.outfile_name} with {self.result.sizes["time"]} timesteps')
        encoding = {'trajectory': {'dtype': self.result.trajectory.attrs['dtype'], '_FillValue': None},
                    'time': {'units': 'seconds since 1970-01-01 00:00:00',  # Expected by some downstream clients
                             'dtype': self.result.time.attrs['dtype'], '_FillValue': None}}
        for varname, var in self.result.variables.items():
            attrs = var.attrs.copy()
            if varname not in ['time', 'trajectory'] and 'dtype' in attrs and issubclass(attrs['dtype'], np.integer):
                FillValue = np.iinfo(attrs['dtype']).max
                encoding[varname] = {'dtype': attrs['dtype'], '_FillValue': FillValue}
            attrs.pop('dtype', None)
            self.result[varname].attrs = attrs
        self.result.to_netcdf(self.outfile_name, unlimited_dims={'time': True},
                              encoding=encoding)
        self._netCDF_encoding = encoding
        return

    self.outfile = Dataset(self.outfile_name, 'a')  # Re-open file at each write
    numtimes = self.outfile['time'].shape[0]

    for varname in self.result.data_vars:
        if 'time' in self.result[varname].dims:
            var = self.outfile.variables[varname]
            var[:, numtimes:numtimes + self.result.sizes['time']] = self.result[varname]
    self.outfile.variables['time'][numtimes:numtimes + self.result.sizes['time']] = \
        date2num(pd.to_datetime(self.result.time).to_pydatetime(),
                 self.outfile['time'].units, self.outfile['time'].calendar)
    self.outfile.sync()  # Flush from memory to disk
    self.outfile.close()  # close file temporarily

    logger.debug(f'Appended {numtimes} steps to file {self.outfile_name}')

def close(self):

    logger.debug('Closing netCDF-file')

    self.outfile = Dataset(self.outfile_name, 'a')
    for var in self.result.data_vars:  # Updating variable attributes, if changed during simulation
        for atn, atv in self.result[var].attrs.items():
            if atn != '_FillValue':
                self.outfile[var].setncattr(atn, atv)
    for atn, atv in self.result.attrs.items():  # Updating global attributes
        self.outfile.setncattr(atn, atv)
    self.outfile.close()  # Finally close file
    logger.debug('Closed netCDF-file')

    # Re-opening (lazy) and storing dataset as self.result
    self.result = xr.open_dataset(self.outfile_name)

    if self.num_elements_scheduled() > 0:
        logger.info(f'Removing {self.num_elements_scheduled()} unseeded elements already written to file')
        seeded_indices = [n for n in np.arange(self.num_elements_total())
                     if n not in self.elements_scheduled.ID]
        self.result = self.result.isel(trajectory=seeded_indices)

    # Finally changing UNLIMITED time dimension to fixed, for CDM compliance.
    compression = {'zlib': True, 'complevel': 6}
    logger.debug(f'Making netCDF file CDM compliant with fixed dimensions, and compressing with {compression}')
    for varname in self.result.data_vars:
        if varname in self._netCDF_encoding:
            self._netCDF_encoding[varname].update(compression)
        else:
            self._netCDF_encoding[varname] = compression
    self.result.to_netcdf(self.outfile_name + '_tmp', unlimited_dims={}, encoding=self._netCDF_encoding)
    self.result.close()  # Closing so that tmp-file can be renamed, thereafter opening lazily again
    shutil.move(self.outfile_name + '_tmp', self.outfile_name)  # Replace original
    self.result = xr.open_dataset(self.outfile_name)

def import_file(self, filename):
    """Create OpenDrift object from imported file.
    """

    logger.debug('Importing from ' + filename)

    self.result = xr.open_dataset(filename)

    self.steps_output = self.result.sizes['time']
    self.start_time = datetime_from_datetime64(self.result.time[0])
    self.time = datetime_from_datetime64(self.result.time[-1])
    self.status_categories = self.result.status.flag_meanings.split()

    num_elements = self.result.sizes['trajectory']
    index_of_first, index_of_last = self.index_of_first_and_last()

    kwargs = {}
    for var in self.result.data_vars:
        if var in self.ElementType.variables:
            last_vals_var = xr.apply_ufunc(
                lambda arr, i: arr[i],
                self.result[var], index_of_last,
                input_core_dims=[['time'], []],
                output_core_dims=[[]],
                vectorize=True, dask='parallelized',
                output_dtypes=[self.result[var].dtype]
            )
            kwargs[var] = last_vals_var.compute()
    kwargs['ID'] = np.arange(num_elements)
    self.elements = self.ElementType(**kwargs)
    self.elements_deactivated = self.ElementType()

    # Remove elements which are scheduled for deactivation
    self.remove_deactivated_elements()

    # Import and apply config settings
    attributes = self.result.attrs
    self.mode = Mode.Config  # To allow setting config
    for key, value in attributes.items():
        if key.startswith('config_'):
            conf_key = key[7:]
            if value == 'True':
                value = True
            if value == 'False':
                value = False
            if value == 'None':
                value = None
            try:
                self.set_config(conf_key, value)
                logger.debug('Setting imported config: %s -> %s' %
                             (conf_key, value))
            except Exception as e:
                logger.warning(e)
                logger.warning('Could not set config: %s -> %s' %
                                (conf_key, value))
    self.mode = Mode.Result

    self.status_categories = self.result.status.flag_meanings.split()
    if 'origin_marker' in self.result.data_vars:
        if 'flag_meanings' in self.result.origin_marker.attrs:
            self.origin_marker = [s.replace('_', ' ') for s in self.result.origin_marker.flag_meanings.split()]


    # Import time steps from metadata
    def timedelta_from_string(timestring):
        if 'day' in timestring:
            days = int(timestring.split('day')[0])
            hs = timestring.split(' ')[-1]
            th = datetime.strptime(hs, '%H:%M:%S')
            return timedelta(days=days, hours=th.hour, minutes=th.minute, seconds=th.second)
        else:
            t = datetime.strptime(timestring, '%H:%M:%S')
            return timedelta(
                hours=t.hour, minutes=t.minute, seconds=t.second)
    try:
        self.time_step = timedelta_from_string(self.result.time_step_calculation)
        self.time_step_output = timedelta_from_string(self.result.time_step_output)
    except Exception as e:
        logger.warning(e)
        logger.warning('Could not parse time_steps from netCDF file')
