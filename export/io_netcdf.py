import sys
import datetime
import string
from shutil import move
import logging

import numpy as np
from netCDF4 import Dataset, num2date, date2num

# Module with functions to export/import trajectory data to/from netCDF file
# Strives to be compliant with netCDF CF-convention on trajectories
# http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#idp8377728
# https://geo-ide.noaa.gov/wiki/index.php?title=NODC_NetCDF_Trajectory_Template

skip_parameters = ['ID']  # Do not write to file


def init(self, filename, times=None):

    self.outfile_name = filename
    self.outfile = Dataset(filename, 'w')
    self.outfile.createDimension('trajectory', len(self.elements.lon))
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
    self.outfile.history = 'Created ' + str(datetime.datetime.now())
    self.outfile.source = 'Output from simulation with OpenDrift'
    self.outfile.model_url = 'https://github.com/knutfrode/opendrift'
    self.outfile.readers = str(self.readers.keys())
    self.outfile.time_coverage_start = str(self.start_time)
    self.outfile.time_step = str(self.time_step)

    # Add all element properties as variables
    for prop in self.elements.variables:
        if prop in skip_parameters:
            continue
        # Note: Should use 'f8' if 'f4' is not accurate enough,
        #       at expense of larger files
        try:
            dtype = self.elements.variables[prop]['dtype']
        except:
            dtype = 'f4'
        var = self.outfile.createVariable(prop, dtype, ('trajectory', 'time'))
        for subprop in self.elements.variables[prop].items():
            if subprop[0] not in ['dtype', 'constant', 'default']:
                # Apparently axis attribute shall not be given for lon and lat:
                if prop in ['lon', 'lat'] and subprop[0] == 'axis':
                    continue
                var.setncattr(subprop[0], subprop[1])
                var.setncattr('coordinates', 'lat lon time')
    # list and number all readers


def write_buffer(self):
    num_steps_to_export = self.steps + 1 - self.steps_exported
    for prop in self.elements.variables:
        if prop in skip_parameters:
            continue
        var = self.outfile.variables[prop]
        var[:, self.steps_exported:self.steps_exported+num_steps_to_export] = \
            self.history[prop][:, 0:num_steps_to_export]

    logging.info('Wrote %s steps to file %s' % (num_steps_to_export,
                                                self.outfile_name))
    #self.history.mask = True  # Reset history array, for new data
    self.steps_exported = self.steps_exported + num_steps_to_export

    # for each element:
    #   - should be possible to select which to write and not
    # - properties, including ID
    # - environment variables, including from which reader
    # - which reader for which


def close(self):

    # Write status categories metadata
    status_dtype = self.ElementType.variables['status']['dtype']
    self.outfile.variables['status'].valid_range = np.array([0]).astype(
        status_dtype), \
        np.array([len(self.status_categories) - 1]).astype(status_dtype)
    self.outfile.variables['status'].flag_values = \
        np.array(np.arange(len(self.status_categories)), dtype=status_dtype)
    self.outfile.variables['status'].flag_meanings = \
        string.join(self.status_categories)
    # Write timesteps to file
    self.outfile.time_coverage_end = str(self.time)
    timeStr = 'seconds since 1970-01-01 00:00:00'
    times = [self.start_time + n*self.time_step for n in range(self.steps + 1)]
    self.outfile.variables['time'][0:len(times)] = date2num(times, timeStr)
    self.outfile.variables['time'].units = timeStr
    self.outfile.variables['time'].standard_name = 'time'
    self.outfile.variables['time'].long_name = 'time'
    # Apparently axis attribute shall not be given for time, lon and lat
    #self.outfile.variables['time'].axis = 'T'

    # Write bounds metadata
    self.outfile.geospatial_lat_min = self.history['lat'].min()
    self.outfile.geospatial_lat_max = self.history['lat'].max()
    self.outfile.geospatial_lat_units = 'degrees_north'
    self.outfile.geospatial_lat_resolution = 'point'
    self.outfile.geospatial_lon_min = self.history['lon'].min()
    self.outfile.geospatial_lon_max = self.history['lon'].max()
    self.outfile.geospatial_lon_units = 'degrees_east'
    self.outfile.geospatial_lon_resolution = 'point'
    self.outfile.runtime = str(self.runtime_environment + self.runtime_model)

    self.outfile.close()  # Finally close file

    # Finally changing UNLIMITED time dimension to fixed, for CDM compliance.
    # Fortunately this is quite fast.
    # http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/FeatureDatasets/CFpointImplement.html
    try:
        with Dataset(self.outfile_name) as src, \
                Dataset(self.outfile_name + '_tmp', 'w') as dst:
            for name, dimension in src.dimensions.iteritems():
                dst.createDimension(name, len(dimension))

            for name, variable in src.variables.iteritems():
                dstVar = dst.createVariable(name, variable.datatype,
                                             variable.dimensions)
                srcVar = src.variables[name]
                dstVar[:] = srcVar[:]  # Copy data
                for att in src.variables[name].ncattrs():
                    # Copy variable attributes
                    dstVar.setncattr(att, srcVar.getncattr(att))

            for att in src.ncattrs():  # Copy global attributes
                dst.setncattr(att, src.getncattr(att))

        move(self.outfile_name + '_tmp', self.outfile_name)  # Replace original
    except:
        print 'Could not convert netCDF file from unlimited to fixed dimension. Could be due to netCDF library incompatibility(?)'
    

def import_file(self, filename, time=None):

    infile = Dataset(filename, 'r')
    self.start_time = num2date(infile.variables['time'][0],
                               infile.variables['time'].units)
    self.end_time = num2date(infile.variables['time'][-1],
                             infile.variables['time'].units)
    self.time_step = num2date(infile.variables['time'][1],
                              infile.variables['time'].units) - self.start_time
    self.time = self.end_time  # Using end time as default
    self.status_categories = infile.variables['status'].flag_meanings.split()

    for var in infile.variables:
        if var not in self.ElementType.variables:
            print var + ' does not exist - adding to element class definition'

    num_elements = len(infile.dimensions['trajectory'])
    num_timesteps = len(infile.dimensions['time'])
    self.steps = num_timesteps - 1  # we do not here count initial state
    dtype = np.dtype([(var[0], var[1]['dtype'])
                      for var in self.ElementType.variables.items()])

    # Import whole dataset (history)
    self.history = np.ma.array(np.zeros([num_elements, num_timesteps]),
                               dtype=dtype, mask=[True])
    for var in infile.variables:
        if var in ['time', 'trajectory']:
            continue
        self.history[var] = infile.variables[var][:, :]

    # Initialise elements from given (or last) state/time
    indx = -1  # using last as default - to be fixed
    lon = self.history['lon'][:, :]
    index_of_last = (self.history['status'] == 0).sum(axis=1)
    index_of_last[index_of_last >= num_timesteps] = num_timesteps - 1
    stat = self.history['status'][np.arange(len(index_of_last)), index_of_last]
    kwargs = {}
    for var in infile.variables:
        if var in ['time', 'trajectory']:
            continue
        kwargs[var] = self.history[var][np.arange(len(index_of_last)),
                                        index_of_last]
    # Import element IDs, which are named 'trajectory' in netCDF CF convention
    kwargs['ID'] = infile.variables['trajectory'][:]
    self.elements = self.ElementType(**kwargs)

    # Remove elements which are scheduled for deactivation
    self.remove_deactivated_elements()
