import sys
import datetime
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
    self.outfile.createVariable('elementID', 'i4', ('trajectory',))
    self.outfile.createDimension('obs', None)  # Unlimited time dimension
    self.outfile.createVariable('time', 'f8', ('obs',))

    self.outfile.Conventions = 'CF-1.6'
    self.outfile.standard_name_vocabulary = 'CF-1.6'
    self.outfile.featureType = 'trajectory'
    self.outfile.history = 'Created ' + str(datetime.datetime.now())
    self.outfile.source = 'Output from simulation with OpenDrift'
    self.outfile.model_url = 'https://github.com/knutfrode/opendrift'
    self.outfile.readers = str(self.readers.keys())
    self.outfile.time_coverage_start = str(self.start_time)
    self.outfile.time_step = str(self.time_step)

    self.outfile.createVariable('crs', 'i4')
    self.outfile.variables['crs'].grid_mapping_name = 'latitude_longitude'
    self.outfile.variables['crs'].epsg_code = 'EPSG:4326'
    self.outfile.variables['crs'].semi_major_axis = 6378137.0
    self.outfile.variables['crs'].inverse_flattening = 298.257223563

    # Add all element properties as variables
    for prop in self.elements.variables:
        if prop in skip_parameters: continue
        # Note: Should use 'f8' if 'f4' is not accurate enough,
        #       at expense of larger files
        var = self.outfile.createVariable(prop, 'f4', ('trajectory', 'obs'))
        for subprop in self.elements.variables[prop].items():
            if subprop[0] not in ['dtype', 'constant']:
                var.setncattr(subprop[0], subprop[1])
    # list and number all readers

def write_buffer(self):
    num_steps_to_export = self.steps + 1 - self.steps_exported
    for prop in self.elements.variables:
        if prop in skip_parameters: continue
        var = self.outfile.variables[prop]
        var[:,self.steps_exported:self.steps_exported+num_steps_to_export] = \
            self.history[prop][:,0:num_steps_to_export]
        
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

    # Write timesteps to file
    self.outfile.time_coverage_end = str(self.time)
    timeStr = 'seconds since 1970-01-01 00:00:00'
    times = [self.start_time + n*self.time_step for n in range(self.steps + 1)]
    self.outfile.variables['time'][:] = date2num(times, timeStr)
    self.outfile.variables['time'].units = timeStr
    self.outfile.variables['time'].axis = 'T'
    self.outfile.close()  # Finally close file

def import_file(self, filename, time=None):
    
    infile = Dataset(filename, 'r')
    self.start_time = num2date(infile.variables['time'][0],
                               infile.variables['time'].units)
    self.end_time = num2date(infile.variables['time'][-1],
                             infile.variables['time'].units)
    self.time_step = num2date(infile.variables['time'][1],
                              infile.variables['time'].units) - self.start_time
    self.time = self.end_time  # Using end time as default

    for var in infile.variables:
        if var not in self.ElementType.variables:
            print var + ' does not exist - adding to element class definition'

    num_elements = len(infile.dimensions['trajectory'])
    num_timesteps = len(infile.dimensions['obs'])
    self.steps = num_timesteps - 1  # we do not here count initial state
    dtype = np.dtype([(var[0],var[1]['dtype'])
                for var in self.ElementType.variables.items()])

    # Import whole dataset (history)
    self.history = np.ma.array(np.zeros([num_elements, num_timesteps]),
                               dtype=dtype, mask=[True])
    for var in infile.variables:
        if var in ['elementID', 'time', 'crs']:
            continue
        self.history[var] = infile.variables[var][:,:]

    # Initialise elements from given (or last) state/time
    indx = -1  # using last as default - to be fixed
    lon = self.history['lon'][:,:]
    index_of_last = (self.history['status']==0).sum(axis=1)
    index_of_last[index_of_last>=num_timesteps] = num_timesteps - 1
    stat = self.history['status'][np.arange(len(index_of_last)), index_of_last]
    kwargs = {}
    for var in infile.variables:
        if var in ['elementID', 'time', 'crs']:
            continue
        kwargs[var] = self.history[var][np.arange(len(index_of_last)),
                                        index_of_last]
    kwargs['ID'] = np.arange(num_elements) + 1  # Should rader be saved/read from file TBD
    self.elements = self.ElementType(**kwargs)

    # Remove elements which are schedules for deactivation
    self.remove_deactivated_elements()
