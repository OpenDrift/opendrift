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

class OpenDriftWrapper:
    """Wraps the Opendrift model for one implementation.
    Each implementation corresponds to a specific Opendrift
    simulation. Simulations details are are defined in a 
    self.configuration file.

    Arguments:
    imp : implementation name 
    self.config_yaml_file: self.configuration file with simulation details
    rootdir: the final data root directory
    rundir: The run directory. A.k.a scratch directory

    """
    def __init__(self,
                 imp = None,
                 self.config_yaml_file = None,
                 rootdir = None
                 **kwargs):
        super(ERCoreWrapper, self).__init__()
        # should I define all 

    def read_yaml_self.config(self):
        """reads a YAML self.config file for OpenDrift, and set self.config attribute of
            OpendriftWrapper object
            see description of YAML self.configuration file
            in XXXX 
        """
        try:
            logging.debug('WRAPPER : Reading self.config file: ' + self.self.config_yaml_file)
            self.config = yaml.load(open(self.config_yaml_file).read())
            return self.config
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
        if isinstance(lon,list) or isinstance(lat,list) :
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
        """ create correct (native) OpenDrift object 'opendrift_obj' based on specs in configuration
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
                    resolution=self.config['basemap_resolution']['resolution'], 
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
            Note this overwrites any pre-existing/default parameter's values.
         """
        logging.debug('WRAPPER : setting configuration')
        for ic,category in enumerate(self.config['config']):
            for ip,param in enumerate(self.config['config'][category]):
                logging.debug('WRAPPER : CONFIG  %s:%s,%s' % (category,param,self.config['config'][category][param]))
                self.opendrift_obj.set_config('%s:%s' % (category,param), self.config['config'][category][param])

    def run_opendrift_from_config(self):
        """ Setup and run opendrit simulation based on input self.configuration.
            self.configuration info is contained in the self.self.config dictionary, which
            was loaded by self.read_yaml_config()
        """
        #############################################################################################
        # release position and timing
        #############################################################################################
        self.release_position_timing()
        #############################################################################################
        # model selection 
        #############################################################################################
        self.model_selection()
        #############################################################################################
        # readers"
        #############################################################################################
        self.set_basemap_reader() # shoreline basemap
        self.set_readers() # ocean/wave/wind readers from various sources
        self.set_fallback_values
        #############################################################################################
        # particle seeding 
        #############################################################################################
        self.seed_particles()
        #############################################################################################
        # set all additional self.configuration - overwriting existing values)
        #############################################################################################       
        self.set_config()
        #############################################################################################
        # start run
        #############################################################################################
        logging.debug('WRAPPER : setting model run')
        time_step = timedelta(seconds = self.config['time_step_sec']) 
        if 'time_step_sec_output' in self.config.keys():
            time_step_output = timedelta(seconds=self.config['time_step_sec_output'])
        else:
            time_step_output = timedelta(minutes=30)
        # duration of run can be specifed by duration or end_time_run
        if 'duration_hours' in self.config.keys():
            duration = timedelta(hours=self.config['duration_hours'])
            end_time_run = start_time + duration
        elif 'end_time_run' in self.config.keys():
            end_time_run = datetime.strptime(self.config['end_time_run'],'%d-%m-%Y %H:%M')      
        else: # if none input , just use end_time which is end of release
            end_time_run = end_time

        if self.config['run_backwards']: # run model backwards 
            time_step = -time_step
            # here we should have start_time>end_time
            if start_time < end_time_run:
                logging.error('WRAPPER : start_time must be smaller than end_time if run backwards')
        # prepare output folder if needed
        if 'rootdir' in self.config.keys():
            if not os.path.exists(self.config['rootdir']):
                os.mkdir(self.config['rootdir'])
        else:
            self.config['rootdir'] = os.getcwd()

        if 'outfile' not in self.config.keys():
            outfile = os.path.join(self.config['rootdir'], 'opendrift_' + self.config['imp'] + '_' + start_time.strftime('%d%m%Y_%H%M') + '.nc')
        else:
            outfile = os.path.join(self.config['rootdir'],self.config['outfile'])
        
        if 'stop_on_error' not in self.config.keys():
            self.config['stop_on_error'] = False
        
        # run simulation
        o.run(stop_on_error = self.config['stop_on_error'],time_step=time_step,
              end_time = end_time_run, time_step_output = time_step_output,
              outfile=outfile)
              #could use **extra_args like other initialization if needed
        return o

    ########################################
    # Scheduler specific functions
    def set_mpiexec(self):
        pass
    def set_ncores(self):
        pass
    def set_cycle(self):
        pass
    def set_progress(self):
        pass
    def run(self):
        pass
    ########################################

if __name__ == '__main__':

   > allow program to be used as a class or as a script , using a config file as input
        e.g. opendriftwrapper.py config.yml just like example_run_config.py
   > need to add a set of download methods that will read the config of object and download necessary data

    parser = argparse.ArgumentParser()
    parser.add_argument('self.config_file',
                        help='<self.configuration file (YAML-format) with all required infos to run Opendrift>')
    args = parser.parse_args()

    # setup opendrift simulation object using specifications in YAML self.configuration file
    self.config = read_yaml_self.config(args.self.config_file) 

    # download data based on self.config file - if necessary
    if any('udshost' in self.config['readers'][block].keys() for block in self.config['readers'].keys()) :
        # some UDS download required
        logging.debug('WRAPPER : Some UDS download required - calling download_metocean_uds.download')
        self.config = download_metocean_uds.download(self.config_dict = self.config)
    if any('cmems_download' in self.config['readers'][block].keys() for block in self.config['readers'].keys()) :
        # some CMEMS download required
        logging.debug('WRAPPER : Some CMEMS download required - calling download_cmems.download')
        self.config = download_cmems.download(self.config_dict = self.config)
    
    # run simulation
    o = run_opendrift_from_self.config(self.config)
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
