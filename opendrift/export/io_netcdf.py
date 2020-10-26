import sys
import os
from datetime import datetime, timedelta
import string
from shutil import move

import numpy as np
from netCDF4 import Dataset, num2date, date2num

# Module with functions to export/import trajectory data to/from netCDF file
# Strives to be compliant with netCDF CF-convention on trajectories
# https://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#idp8377728
# https://geo-ide.noaa.gov/wiki/index.php?title=NODC_NetCDF_Trajectory_Template

skip_parameters = ['ID']  # Do not write to file


def init(self, filename):

    self.outfile_name = filename
    self.outfile = Dataset(filename, 'w')
    self.outfile.createDimension('trajectory', self.num_elements_total())
    self.outfile.createVariable('trajectory', 'i4', ('trajectory',))
    self.outfile.createDimension('time', None)  # Unlimited time dimension
    self.outfile.createVariable('time', 'f8', ('time',))
    # NB: trajectory_id must be changed for future ragged array representation
    self.outfile.variables['trajectory'][:] = \
        np.arange(self.num_elements_total())+1
    self.outfile.variables['trajectory'].cf_role = 'trajectory_id'
    self.outfile.variables['trajectory'].units = '1'

    self.outfile.Conventions = 'CF-1.6'
    self.outfile.standard_name_vocabulary = 'CF-1.6'
    self.outfile.featureType = 'trajectory'
    self.outfile.history = 'Created ' + str(datetime.now())
    self.outfile.source = 'Output from simulation with OpenDrift'
    self.outfile.model_url = 'https://github.com/OpenDrift/opendrift'
    self.outfile.opendrift_class = self.__class__.__name__
    self.outfile.opendrift_module = self.__class__.__module__
    self.outfile.readers = str(self.readers.keys())
    self.outfile.time_coverage_start = str(self.start_time)
    self.outfile.time_step_calculation = str(self.time_step)
    self.outfile.time_step_output = str(self.time_step_output)

    # Time
    self.timeStr = 'seconds since 1970-01-01 00:00:00'
    self.outfile.variables['time'].units = self.timeStr
    self.outfile.variables['time'].standard_name = 'time'
    self.outfile.variables['time'].long_name = 'time'
    # Apparently axis attribute shall not be given for time, lon and lat
    #self.outfile.variables['time'].axis = 'T'

    # Write config settings
    for key in self._config:
        value = self.get_config(key)
        if isinstance(value, (bool, type(None))):
            value = str(value)
        self.outfile.setncattr('config_' + key, value)

    # Write additionaly metadata attributes, if given
    if hasattr(self, 'metadata_dict'):
        for key, value in self.metadata_dict.items():
            self.outfile.setncattr(key, str(value))

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
        var = self.outfile.createVariable(prop, dtype, ('trajectory', 'time'))
        var.setncattr('coordinates', 'lat lon time')
        for subprop in self.history_metadata[prop].items():
            if subprop[0] not in ['dtype', 'constant', 'default']:
                # Apparently axis attribute shall not be given for lon and lat:
                if prop in ['lon', 'lat'] and subprop[0] == 'axis':
                    continue
                var.setncattr(subprop[0], subprop[1])

def write_buffer(self):
    if self.outfile._isopen == 0:
        self.outfile = Dataset(self.outfile_name, 'a')
    num_steps_to_export = self.steps_output - self.steps_exported
    for prop in self.history_metadata:
        if prop in skip_parameters:
            continue
        var = self.outfile.variables[prop]
        var[:, self.steps_exported:self.steps_exported+num_steps_to_export] = \
            self.history[prop][:, 0:num_steps_to_export]

    times = [self.start_time + n*self.time_step_output for n in
             range(self.steps_exported, self.steps_output)]
    self.outfile.variables['time'][
        self.steps_exported:self.steps_exported+len(times)] = \
            date2num(times, self.timeStr)

    # Write status categories metadata
    status_dtype = self.ElementType.variables['status']['dtype']
    self.outfile.variables['status'].valid_range = np.array(
        (0, len(self.status_categories) - 1)).astype(status_dtype)
    self.outfile.variables['status'].flag_values = \
        np.array(np.arange(len(self.status_categories)), dtype=status_dtype)
    self.outfile.variables['status'].flag_meanings = \
        " ".join(self.status_categories)

    self.logger.info('Wrote %s steps to file %s' % (num_steps_to_export,
                                                self.outfile_name))
    self.history.mask = True  # Reset history array, for new data
    self.steps_exported = self.steps_exported + num_steps_to_export
    self.outfile.steps_exported = self.steps_exported
    self.outfile.sync()  # Flush from memory to disk
    self.outfile.close()  # close file temporarily

def close(self):
    self.outfile = Dataset(self.outfile_name, 'a')
    # Write status categories metadata
    status_dtype = self.ElementType.variables['status']['dtype']
    self.outfile.variables['status'].valid_range = np.array(
        (0, len(self.status_categories) - 1)).astype(status_dtype)
    self.outfile.variables['status'].flag_values = \
        np.array(np.arange(len(self.status_categories)), dtype=status_dtype)
    self.outfile.variables['status'].flag_meanings = \
        " ".join(self.status_categories)
    # Write final timesteps to file
    self.outfile.time_coverage_end = str(self.time)

    # Write performance data
    self.outfile.performance = self.performance()

    # Write bounds metadata
    self.outfile.geospatial_lat_min = self.history['lat'].min()
    self.outfile.geospatial_lat_max = self.history['lat'].max()
    self.outfile.geospatial_lat_units = 'degrees_north'
    self.outfile.geospatial_lat_resolution = 'point'
    self.outfile.geospatial_lon_min = self.history['lon'].min()
    self.outfile.geospatial_lon_max = self.history['lon'].max()
    self.outfile.geospatial_lon_units = 'degrees_east'
    self.outfile.geospatial_lon_resolution = 'point'
    self.outfile.runtime = str(datetime.now() -
                               self.timers['total time'])

    self.outfile.close()  # Finally close file

    # Finally changing UNLIMITED time dimension to fixed, for CDM compliance.
    # Fortunately this is quite fast.
    # https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/FeatureDatasets/CFpointImplement.html
    try:
        self.logger.debug('Making netCDF file CDM compliant with fixed dimensions')
        if self.num_elements_scheduled() > 0:
            self.logger.info('Removing %i unseeded elements already written to file' % self.num_elements_scheduled())
            mask = np.ones(self.history.shape[0], dtype=bool)
            mask[self.elements_scheduled.ID-1] = False
        with Dataset(self.outfile_name) as src, \
                Dataset(self.outfile_name + '_tmp', 'w') as dst:
            for name, dimension in src.dimensions.items():
                if name=='trajectory':
                    # Truncate dimension length to  number actually seeded
                    dst.createDimension(name, self.num_elements_activated())
                else:
                    dst.createDimension(name, len(dimension))

            for name, variable in src.variables.items():
                dstVar = dst.createVariable(name, variable.datatype,
                                             variable.dimensions)
                srcVar = src.variables[name]
                # Truncate data to number actually seeded
                if 'trajectory' in variable.dimensions:
                    if self.num_elements_scheduled() > 0:
                        if len(variable.dimensions) == 2:
                            dstVar[:] = srcVar[mask, :] 
                        else:
                            dstVar[:] = srcVar[mask]  # Copy data
                    else:
                        dstVar[:] = srcVar[:]
                else:
                    dstVar[:] = srcVar[:]
                for att in src.variables[name].ncattrs():
                    # Copy variable attributes
                    dstVar.setncattr(att, srcVar.getncattr(att))

            for att in src.ncattrs():  # Copy global attributes
                dst.setncattr(att, src.getncattr(att))

        move(self.outfile_name + '_tmp', self.outfile_name)  # Replace original
    except Exception as me:
        print(me)
        print('Could not convert netCDF file from unlimited to fixed dimension. Could be due to netCDF library incompatibility(?)')
    
def import_file_xarray(self, filename, chunks):

    import xarray as xr
    self.logger.debug('Importing with Xarray from ' + filename)
    self.ds = xr.open_dataset(filename, chunks=chunks)

    self.steps_output = len(self.ds.time)
    ts0 = (self.ds.time[0] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    self.start_time = datetime.utcfromtimestamp(np.float(ts0))
    tse = (self.ds.time[-1] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    self.end_time = datetime.utcfromtimestamp(np.float(tse))
    if len(self.ds.time) > 1:
        ts1 = (self.ds.time[1] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
        self.time_step_output = timedelta(seconds=np.float(ts1 - ts0))
    self.time = self.end_time  # Using end time as default
    self.status_categories = self.ds.status.flag_meanings.split()

    num_elements = len(self.ds.trajectory)
    elements=np.arange(num_elements)

    # Data types for variables
    dtype = np.dtype([(var[0], var[1]['dtype'])
                      for var in self.ElementType.variables.items()])
    history_dtype_fields = [
        (name, self.ElementType.variables[name]['dtype'])
        for name in self.ElementType.variables]
    # Add environment variables
    self.history_metadata = self.ElementType.variables.copy()
    for env_var in self.required_variables:
        if env_var in self.ds.variables:
            history_dtype_fields.append((env_var, np.dtype('float32')))
            self.history_metadata[env_var] = {}
    history_dtype = np.dtype(history_dtype_fields)

    # Masking where elements are not been seeded
    # TODO: mask from deactivation towards end
    for da in ['lon', 'lat']:
        self.ds[da] = self.ds[da].where(self.ds.status>=0)

    # Read some saved parameters
    if os.path.exists(self.analysis_file):
        self.af = xr.open_dataset(self.analysis_file)
        self.lonmin = self.af.lonmin
        self.lonmax = self.af.lonmax
        self.latmin = self.af.latmin
        self.latmax = self.af.latmax

def import_file(self, filename, times=None, elements=None):

    self.logger.debug('Importing from ' + filename)
    infile = Dataset(filename, 'r')
    # 'times' can be used to import subset. Not yet implemented.
    if times is None and hasattr(infile, 'steps_exported'):
        self.steps_output = infile.steps_exported
    else:
        self.steps_output = len(infile.dimensions['time'])
    
    self.start_time = num2date(infile.variables['time'][0],
                               infile.variables['time'].units)
    if len(infile.variables['time']) > 1:
        self.end_time = num2date(
            infile.variables['time'][self.steps_output-1],
            infile.variables['time'].units)
        self.time_step_output = num2date(infile.variables['time'][1],
            infile.variables['time'].units) - self.start_time
    else:
        self.time_step_output = timedelta(hours=1)
        self.end_time = self.start_time
    self.time = self.end_time  # Using end time as default
    self.status_categories = infile.variables['status'].flag_meanings.split()

    if elements is None:
        num_elements = len(infile.dimensions['trajectory'])
        elements=np.arange(num_elements)
    else:
        num_elements=len(elements)

    dtype = np.dtype([(var[0], var[1]['dtype'])
                      for var in self.ElementType.variables.items()])

    history_dtype_fields = [
        (name, self.ElementType.variables[name]['dtype'])
        for name in self.ElementType.variables]
    # Add environment variables
    self.history_metadata = self.ElementType.variables.copy()
    for env_var in self.required_variables:
        if env_var in infile.variables:
            history_dtype_fields.append((env_var, np.dtype('float32')))
            self.history_metadata[env_var] = {}
    history_dtype = np.dtype(history_dtype_fields)

    # Import dataset (history)
    self.history = np.ma.array(
        np.zeros([num_elements, self.steps_output]),
        dtype=history_dtype, mask=[True])
    for var in infile.variables:
        if var in ['time', 'trajectory']:
            continue
        try:
            self.history[var] = infile.variables[var][elements, 0:self.steps_output]
        except Exception as e:
            self.logger.info(e)
            pass

    # Initialise elements from given (or last) state/time
    firstlast = np.ma.notmasked_edges(self.history['status'], axis=1)
    index_of_last = firstlast[1][1]
    kwargs = {}
    for var in infile.variables:
        if var in self.ElementType.variables:
            kwargs[var] = self.history[var][
                np.arange(len(index_of_last)), index_of_last]
    kwargs['ID'] = elements + 1
    self.elements = self.ElementType(**kwargs)
    self.elements_deactivated = self.ElementType()

    # Remove elements which are scheduled for deactivation
    self.remove_deactivated_elements()

    # Import and apply config settings
    attributes = infile.ncattrs()
    for attr in attributes:
        if attr.startswith('config_'):
            value = infile.getncattr(attr)
            conf_key = attr[7:]
            if value == 'True':
                value = True
            if value == 'False':
                value = False
            if value == 'None':
                value = None
            try:
                self.set_config(conf_key, value)
                self.logger.debug('Setting imported config: %s -> %s' %
                             (conf_key, value))
            except:
                self.logger.warning('Could not set config: %s -> %s' %
                                (conf_key, value))

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
        self.time_step = timedelta_from_string(infile.time_step_calculation)
        self.time_step_output = timedelta_from_string(infile.time_step_output)
    except Exception as e:
        self.logger.warning(e)
        self.logger.warning('Could not parse time_steps from netCDF file')

    infile.close()
