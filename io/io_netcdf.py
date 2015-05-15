import sys
import logging
from netCDF4 import Dataset


def init(self, filename, steps=None):
    
    self.outfile_name = filename
    self.outfile = Dataset(filename, 'w')
    self.outfile.createDimension('elementID', len(self.elements.lon))
    self.outfile.createVariable('elementID', 'i4', ('elementID',))
    self.outfile.createDimension('step', steps)
    self.outfile.createVariable('time', 'f4', ('step',))

    # Add all element properties as variables
    for prop in self.elements.variables:
        if prop == 'ID': continue
        self.outfile.createVariable(prop, 'f4', ('elementID', 'step'))
    # list and number all readers

def write_buffer(self):
    num_steps_to_export = self.steps + 1 - self.steps_exported
    for prop in self.elements.variables:
        if prop == 'ID': continue
        var = self.outfile.variables[prop]
        var[:,self.steps_exported:self.steps_exported+num_steps_to_export] = \
            self.history[prop][:,0:num_steps_to_export]
        
    logging.info('Wrote %s steps to file %s' % (num_steps_to_export,
                                                 self.outfile_name))
    self.history.mask = True  # Reset history array, for new data
    self.steps_exported = self.steps_exported + num_steps_to_export

    # for each element:
    #   - should be possible to select which to write and not 
    # - properties, including ID 
    # - environment variables, including from which reader
    # - which reader for which 

def close(self):
    self.outfile.close()  # close file

def import_file(self, filename, time=None):
    pass
