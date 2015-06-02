import sys
import datetime
import logging
from netCDF4 import Dataset, num2date, date2num

# Module with functions to export/import trajectory data to/from netCDF file
# Strives to be compliant with netCDF CF-convention on trajectories
# http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#idp8377728
# https://geo-ide.noaa.gov/wiki/index.php?title=NODC_NetCDF_Trajectory_Template

skip_parameters = ['ID', 'status']  # Do not write to file

def init(self, filename, times=None):
    
    self.outfile_name = filename
    self.outfile = Dataset(filename, 'w')
    self.outfile.createDimension('elementID', len(self.elements.lon))
    self.outfile.createVariable('elementID', 'i4', ('elementID',))
    self.outfile.createDimension('timeID', None)  # Unlimited time dimension
    self.outfile.createVariable('time', 'f4', ('timeID',))

    self.outfile.Conventions = 'CF-1.6'
    self.outfile.standard_name_vocabulary = 'CF-1.6'
    self.outfile.featureType = 'trajectory'
    self.outfile.history = 'Created ' + str(datetime.datetime.now())
    self.outfile.source = 'Output from simulation with OpenDrift'
    self.outfile.model_url = 'https://github.com/knutfrode/opendrift'
    self.outfile.time_coverage_start = str(self.start_time)
    self.outfile.time_step = str(self.time_step)

    # Add all element properties as variables
    for prop in self.elements.variables:
        if prop in skip_parameters: continue
        var = self.outfile.createVariable(prop, 'f4', ('elementID', 'timeID'))
        for subprop in self.elements.variables[prop].items():
            if subprop[0] != 'dtype':
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
    pass  # To be implemented
