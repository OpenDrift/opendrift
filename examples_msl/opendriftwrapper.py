#!/usr/bin/env python
#
# Script to be called when Opendrift is run using a self.config file 
# 
# Usage :
# python example_run_self.config.py ./self.config/opendrift.self.config_oceandrift3d.yml
#
# Still in dev....
# 
# 
import os
import numpy as np
import argparse
import yaml
import sys
sys.path.append('C:\github/toolbox_simon/python_tools')
import download_metocean_uds
import download_cmems
from datetime import datetime,timedelta
import opendrift
import logging

class OpenDriftWrapper(object):
    """Wraps the Opendrift model for one implementation.
    Each implementation corresponds to a specific Opendrift
    simulation. Simulations details are are defined in a 
    configuration YAML file.

    Arguments:
    config_yaml_file: self.configuration file with simulation details
    rootdir: the final data root directory
    rundir: The run directory. A.k.a scratch directory

    """
    def __init__(self,
                 config_yaml_file = None,
                 pipe=None,
                 logger=logging,
                 schedule={},
                 **kwargs):

        
        self.pipe = pipe        # scheduler-specific
        self.logger = logger    # scheduler-specific
        self.schedule = schedule# scheduler-specific
        self.config_yaml_file = config_yaml_file
        super(OpenDriftWrapper, self).__init__()

    def read_yaml_config(self):
        """reads a YAML self.config file for OpenDrift, and set self.config attribute of
            OpendriftWrapper object
            see description of YAML self.configuration file
            in XXXX 
        """
        try:
            logging.debug('WRAPPER : Reading self.config file: ' + self.config_yaml_file)
            self.config = yaml.load(open(self.config_yaml_file).read())
        except Exception as ex: # handle arbitrary exception
            logging.error(ex)
            logging.error('WRAPPER : Cannot read ' + self.config_yaml_file)
            sys.exit('Cannot read ' + self.config_yaml_file)

    def release_position_timing(self):
        """ get release position and timing details from config"""
        self.start_time = datetime.strptime(self.config['start_time'],'%d-%m-%Y %H:%M')
        self.end_time = datetime.strptime(self.config['end_time'],'%d-%m-%Y %H:%M')
        if self.start_time == self.end_time:
            self.release_time = self.start_time # one-off release
        else:
            self.release_time = [self.start_time, self.end_time] #constant release over time
        self.lon = self.config['position']['lon']
        self.lat = self.config['position']['lat']
        self.radius = self.config['position']['radius']
        self.z = self.config['position']['z']
        self.nb_parts = self.config['nb_parts']
        # check if we should use cone switch
        if isinstance(self.lon,list) or isinstance(self.lat,list) :
            self.cone = True 
            # this required when lon or lat is a 2-item list 
            # see /examples/example_seed_demonstration.py
        else:
            self.cone = False

        if 'end_position' in self.config:
            self.elon = self.config['end_position']['elon']
            self.elat = self.config['end_position']['elat']
            self.eradius = self.config['end_position']['eradius']
            self.ez = self.config['end_position']['ez']
            # check if cone release
            if self.lon != self.elon or self.lat != self.elat or self.start_time != self.end_time:
                self.lon = [self.lon, self.elon]
                self.lat = [self.lat, self.elat]
                self.radius = [self.radius, self.eradius]
                self.release_time = [start_time, end_time]
                self.cone = True
                logging.debug('WRAPPER : using cone-release')
            else:
                self.cone = False

    def model_selection(self):
        """ import correct class and create (native) OpenDrift object 'opendrift_obj' 
            based on specs in configuration
        """
        if 'extra_model_args' not in self.config.keys():
            extra_model_args = {}
            self.config['extra_model_args'] = {}
        else:
            extra_model_args = self.config['extra_model_args']

        # This could probably be made more generic like the readers' loading part with eval statements
        if self.config['model'] == 'Leeway':
            from opendrift.models.leeway import Leeway
            self.opendrift_obj = Leeway(loglevel=0,**extra_model_args)
            # add extra model arguments if relevant
            #extra_seed_args = {'objectType': ln + 1}
        elif self.config['model'] == 'OceanDrift3D':
            from opendrift.models.oceandrift3D import OceanDrift3D
            # o = OceanDrift3D(loglevel=0,logfile=logfilename,**extra_model_args) # issue with logfile...not sure why
            self.opendrift_obj = OceanDrift3D(loglevel=0,**extra_model_args) #logfile=logfilename
        elif self.config['model'] == 'OpenOil3D':
            from opendrift.models.openoil3D import OpenOil3D
            if 'weathering_model' in self.config['extra_model_args'].keys():
                self.opendrift_obj = OpenOil3D(loglevel=0,**extra_model_args)
            else: # use noaa by default
                self.opendrift_obj = OpenOil3D(weathering_model='noaa',
                              loglevel=0,**extra_model_args)          
        elif self.config['model'] == 'ShipDrift':
            from opendrift.models.shipdrift import ShipDrift
            self.opendrift_obj = ShipDrift(loglevel=0,**extra_model_args)
        elif self.config['model'] == 'SedimentDrift3D':
            from opendrift.models.sedimentdrift3D import SedimentDrift3D
            self.opendrift_obj = SedimentDrift3D(loglevel=0,**extra_model_args)

        logging.debug('WRAPPER : ' + self.config['model'] +' initialized')

    def set_basemap_reader(self):
        """ set Basemap reader to be used by opendrift_obj i.e. shoreline"""
        logging.debug('WRAPPER : adding basemap reader...')
        from opendrift.readers import reader_basemap_landmask
        # Making customised landmask (Basemap)
        reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon= self.config['model_frame']['llcrnrlon'], 
                    llcrnrlat=self.config['model_frame']['llcrnrlat'], 
                    urcrnrlon= self.config['model_frame']['urcrnrlon'], 
                    urcrnrlat=self.config['model_frame']['urcrnrlat'],  
                    resolution=self.config['basemap_resolution'], 
                    projection='merc') # use mercator projection by default
        self.opendrift_obj.add_reader(reader_basemap)

    def set_readers(self):
        """ set readers to be used by opendrift_obj """
        logging.debug('WRAPPER : adding ocean/wave/wind readers...')

        # loop through different readers and initialize
        base_str = 'from opendrift.readers import '
        # The order in which readers are specified in the self.config file matters
        # for "priority" i.e. input highest res/smaller extents first then lowest res/larger extents
        for key in self.config['readers'].keys():
            if self.config['readers'][key] is not None:
                if 'filename' in self.config['readers'][key]: 
                # 'filename' is either input by user in self.config file, 
                # or automatically set after data download
                    logging.debug('WRAPPER : adding reader: ' +  key)
                    for read_cnt,fname in enumerate(self.config['readers'][key]['filename']):
                        # import correct reader class
                        reader_type = self.config['readers'][key]['reader_type'][read_cnt]
                        exec(base_str + reader_type + ' as reader_obj',globals(),locals()) # needed to add globals() to make it work
                        # load reader as reader_obj, then initialize it
                        if 'options' in self.config['readers'][key]: # add specific options as **kwargs
                            opts =  self.config['readers'][key]['options']
                        else: # no specific options
                            opts = {}
                        read_tmp = reader_obj.Reader(filename=fname,**opts)
                        #save reader with new name 
                        exec(key + '_' + str(read_cnt) + '= read_tmp',globals(),locals()) 
                        del read_tmp
                        # add reader to OpenDrift object
                        self.opendrift_obj.add_reader(eval( key + '_' + str(read_cnt)) )

    def set_fallback_values(self):
        """ set fallback to be used by opendrift_obj, as specified in config file """
        # adding fallback values / constants
        if 'fallback_values' in self.config.keys():
            for variable in self.config['fallback_values'].keys():
                self.opendrift_obj.fallback_values[variable] = self.config['fallback_values'][variable]

    def seed_particles(self):
        """ seed particles to opendrift_obj """
        # user-input extra_seed_args is passed as **kwargs to seed_elements() function
        # >requires knowledge of model, but highly flexible
        # Note this is located after the forcing fields loading to allow use of 'seafloor' 
        # or 'seafloor+2' for the vertical z seeding 
        logging.debug('WRAPPER : seeding elements')
        if 'extra_seed_args' not in self.config.keys():
            extra_seed_args = {}
            self.config['extra_seed_args'] = {}
        else:
            extra_seed_args = self.config['extra_seed_args'] # e.g. #extra_seed_args = {'objectType': 26}
        
        self.opendrift_obj.seed_elements(lon = self.lon, 
                                         lat = self.lat, 
                                         z = self.z, 
                                         number = self.nb_parts, 
                                         radius = self.radius,
                                         time = self.release_time, 
                                         cone = self.cone,
                                         **extra_seed_args)
  
    def set_config(self):
        """ applies specific configuration parameters to opendrift_obj
            see examples of configuration here:
            https://github.com/OpenDrift/opendrift/blob/master/opendrift/models/basemodel.py#L123
            Note this overwrites any pre-existing/default parameter value.
         """
        logging.debug('WRAPPER : setting configuration')
        for ic,category in enumerate(self.config['config']):
            for ip,param in enumerate(self.config['config'][category]):
                logging.debug('WRAPPER : CONFIG  %s:%s,%s' % (category,param,self.config['config'][category][param]))
                self.opendrift_obj.set_config('%s:%s' % (category,param), self.config['config'][category][param])

    def set_config(self):
        """ applies specific configuration parameters to opendrift_obj
            see examples of configuration here:
            https://github.com/OpenDrift/opendrift/blob/master/opendrift/models/basemodel.py#L123
            Note this overwrites any pre-existing/default parameter's values.
         """
        logging.debug('WRAPPER : setting configuration')
        for ic,category in enumerate(self.config['config']):
            for ip,param in enumerate(self.config['config'][category]):
                logging.debug('WRAPPER : CONFIG  %s:%s,%s' % (category,param,self.config['config'][category][param]))
                self.opendrift_obj.set_config('%s:%s' % (category,param), self.config['config'][category][param])

    def create_rootdir(self):
        ''' create and/or set output folder if needed
        '''
        if 'rootdir' in self.config.keys():
            if not os.path.exists(self.config['rootdir']):
                os.mkdir(self.config['rootdir'])
        else:
            self.config['rootdir'] = os.getcwd()

    def set_run_specs(self):
        '''
        set all parameters required by opendrift_obj.run() method
        '''
        logging.debug('WRAPPER : setting model run')
        self.time_step = timedelta(seconds = self.config['time_step_sec']) 
        if 'time_step_sec_output' in self.config.keys():
            self.time_step_output = timedelta(seconds=self.config['time_step_sec_output'])
        else:
            self.time_step_output = timedelta(minutes=30)
        # duration of run can be specifed by duration or end_time_run
        if 'duration_hours' in self.config.keys():
            self.duration = timedelta(hours=self.config['duration_hours'])
            self.end_time_run = start_time + duration
        elif 'end_time_run' in self.config.keys():
            self.end_time_run = datetime.strptime(self.config['end_time_run'],'%d-%m-%Y %H:%M')      
        else: # if none input , just use end_time which is end of release
            self.end_time_run = end_time

        if self.config['run_backwards']: # run model backwards 
            self.time_step = -self.time_step
            # here we should have start_time>end_time
            if start_time < end_time_run:
                logging.error('WRAPPER : start_time must be smaller than end_time if run backwards')

        if 'outfile' not in self.config.keys():
            self.outfile = os.path.join(self.config['rootdir'], 'opendrift_' + self.config['imp'] + '_' + start_time.strftime('%d%m%Y_%H%M') + '.nc')
        else:
            self.outfile = os.path.join(self.config['rootdir'],self.config['outfile'])
        
        if 'stop_on_error' not in self.config.keys():
            self.config['stop_on_error'] = False

    def setup_opendrift_from_config(self):
        """ Setup opendrit object based on input self.configuration.
            self.configuration info is contained in the self.self.config dictionary, which
            was loaded by self.read_yaml_config()
        """

        self.create_rootdir()  # set and/or create run directory
        self.release_position_timing() #set release position and timing
        self.model_selection() # select appropriate opendrft model (e.g OpenOil, Leeway etc..
        self.set_basemap_reader() # load shoreline basemap reader
        self.set_readers() # ocean/wave/wind readers from various sources
        self.set_fallback_values() # set constant when for variables with no readers
        self.seed_particles() # particle seeding     
        self.set_config() # set all additional self.configuration (overwriting existing/default values)
        self.set_run_specs() # set last run specifications
        self.start_opendrift_run()

    def start_opendrift_run(self):
        self.opendrift_obj.run(stop_on_error = self.config['stop_on_error'],
                           time_step=self.time_step, 
                           end_time = self.end_time_run, 
                           time_step_output = self.time_step_output,
                           outfile=self.outfile)
          #could use **extra_args like other initialization if needed

    def run(self):
        ''' 
        run() method is required by Scheduler to executes the process. This is the method that is 
        explicitely called by Scheduler, or called in __main__ if wrapper is used as standalone.
        It can return a result as a dictionary which can
        be used to update kwargs for a another action using Scheduler Linked Actions.
        '''
        self.read_yaml_config()
        self.download_reader_data() 
        self.setup_opendrift_from_config()
        self.start_opendrift_run() 

    def download_reader_data(self):
        '''
        download necessary data from UDS, CMEMS or other if required/specified in config file
        '''
            # download data based on self.config file - if necessary
        if any('udshost' in self.config['readers'][block].keys() for block in self.config['readers'].keys()) :
            # some UDS download required
            logging.debug('WRAPPER : Some UDS download required - calling download_metocean_uds.download')
            self.config = download_metocean_uds.download(config_dict = self.config)
        if any('cmems_download' in self.config['readers'][block].keys() for block in self.config['readers'].keys()) :
            # some CMEMS download required
            logging.debug('WRAPPER : Some CMEMS download required - calling download_cmems.download')
            self.config = download_cmems.download(config_dict = self.config)
        
        # if any('shom_data' in self.config['readers'][block].keys() for block in self.config['readers'].keys()) :
            # etc....


    ################################################################################
    # Scheduler specific functions
    # https://github.com/metocean/scheduler/blob/cc81ad31db00278157bb9ddcb3b34d68f8527c26/docs/wrappers.rst
    def set_mpiexec(self):
        '''
        (Optional)

        Called by the scheduler before self.run()

        ncores: int - Total number of cores allocated for all the nodes
        hosts: {str:int} - dictionary with allocated nodes ip's as keys
                           and values as number of cores per node.
        '''
        pass
    def set_ncores(self):
        pass
    def set_cycle(self):
        '''
        (Optional)
        This function will be called by the scheduler after instantiation,
        passing a datetime object, in UTC, representing the current Cycle.
        That's how Scheduler tells to your process at which cycle to run.
        '''
        pass
    def set_progress(self):
        '''
        (Optional)
        This method is not actually part of the wrapper API, just an example
        of how you can implement progress communication. For long running
        Actions, you can implement progress by calling the `send` method from
        the `pipe` object which accepts a tuple as argument with 3 values in it

        (partial, total, message)

        '''
        pass

    ################################################################################

    def post_process(self):
        # post-process/plots
        print o
        if self.config['post_process']['show_plot']:
            o.plot()
        if self.config['post_process']['show_anim']:
            o.animation()
        if self.config['post_process']['show_oil_budget']:
            o.plot_oil_budget()
        # if self.config['post_process']['save_animate']:
            # o.animation()
        # if self.config['post_process']['save_plot']:
            # o.plot()
        # if self.config['post_process']['save_oil_budget']:
             # o.plot_oil_budget()
        # MORE TO ADD HERE
        # density, residence time, plot property etc...
        # at least all built in opendrfit function + more

################################################################################

if __name__ == '__main__':
    # this is executed when opendriftwrapper.py is called as a script, passing the config file as argument
    # for example:
    #     python opendriftwrapper.py ./config/opendrift.config_sediment3d_uds_hc.yml


    parser = argparse.ArgumentParser()
    parser.add_argument('config_file',
                        help='<self.configuration file (YAML-format) with all required infos to run Opendrift>')
    args = parser.parse_args()
    # initialize wrapper object with user-input config file
    wrapper = OpenDriftWrapper(config_yaml_file = args.config_file)
    wrapper.run() # the same methord run() will be called by scheduler
    wrapper.post_process()