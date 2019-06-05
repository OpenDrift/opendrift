#!/usr/bin/env python
#
# Script to be called when Opendrift is run using a config file 
# 
# Usage :
# python 
#
# Still in dev....
# 
# 
import numpy as np
import argparse
import yaml
import sys
from datetime import datetime,timedelta
import opendrift
import logging

def read_yaml_config(config_yaml_file):
    ''' read a YAML config file for OpenDrift - could be moved to opendrift/__init__.py eventually '''
    try:
        config = yaml.load(open(config_yaml_file).read())
        return config
    except Exception as ex: # handle arbitrary exception
        logging.error(ex)
        logging.error('Cannot read ' + config_yaml_file)
        sys.exit('Cannot read ' + config_yaml_file)

def run_opendrift_from_config(config):
    ''' Setup and run opendrit simulation based on input configuration
        General process is inspired from run_opendrift() from opendrift_gui.py which fetches 
        all infos from GUI, and then start run
    '''
    #############################################################################################
    # release position and timing
    #############################################################################################
    start_time = datetime.strptime(config['start_time'],'%d-%m-%Y %H:%M')
    end_time = datetime.strptime(config['end_time'],'%d-%m-%Y %H:%M')
    lon = config['position']['lon']
    lat = config['position']['lat']
    radius = config['position']['radius']
    z = config['position']['z']
    nb_parts = config['nb_parts']
    cone = False
    if 'end_position' in config:
        elon = config['end_position']['elon']
        elat = config['end_position']['elat']
        eradius = config['end_position']['eradius']
        ez = config['end_position']['ez']
        # check if cone release
        if lon != elon or lat != elat or start_time != end_time:
            lon = [lon, elon]
            lat = [lat, elat]
            radius = [radius, eradius]
            start_time = [start_time, end_time]
            cone = True
            logging.debug('using cone-release')
        else:
            cone = False
    # we could do something on the z here if required e.g. random within range etc..
    # probably more check to do expand options e.g release along line, single point but continuous over time etc.. 
    # 
    #############################################################################################
    # model selection 
    #############################################################################################
    if 'extra_model_args' not in config.keys():
        extra_model_args = {}
    else:
        extra_model_args = config['extra_model_args']
    logfilename = config['imp'] +config['model']+ '.log'
    
    # This could probably be made more generic like the readers' loading part
    if config['model'] == 'Leeway':
        from opendrift.models.leeway import Leeway
        o = Leeway(loglevel=0,logfile = logfilename,**extra_model_args)
        # add extra model arguments if relevant
        #extra_seed_args = {'objectType': ln + 1}

    elif config['model'] == 'OceanDrift3D':
        from opendrift.models.oceandrift3D import OceanDrift3D
        # o = OceanDrift3D(loglevel=0,logfile=logfilename,**extra_model_args) # issue with logfile...not sure why
        o = OceanDrift3D(loglevel=DEBUG,**extra_model_args)
    elif config['model'] == 'OpenOil3D':
        from opendrift.models.openoil3D import OpenOil3D
        if 'weathering_model' in config.keys():
            o = OpenOil3D(weathering_model=config['weathering_model'],
                          loglevel=0,logfile=logfilename,**extra_model_args)
        else:
            o = OpenOil3D(weathering_model='noaa',
                          loglevel=0,logfile=logfilename,**extra_model_args)          

    elif config['model'] == 'ShipDrift':
        from opendrift.models.shipdrift import ShipDrift
        o = ShipDrift(loglevel=0,logfile=logfilename,**extra_model_args)
    logging.debug(config['model'] +' initialized')
    #############################################################################################
    # particle seeding 
    #############################################################################################
    # user-input extra_seed_args is passed as **kwargs to seed_elements() function
    # >requires knowledge of model, but highly flexible 
    logging.debug('seeding elements')
    if 'extra_seed_args' not in config.keys():
        extra_seed_args = {}
    else:
        extra_seed_args = config['extra_seed_args'] # e.g. #extra_seed_args = {'objectType': 26}
    
    o.seed_elements(lon=lon, lat=lat, z=z, number=nb_parts, radius=radius,
                        time=start_time, cone=cone,
                        **extra_seed_args)

    #############################################################################################
    # forcing data "readers"
    #############################################################################################
    # initial code below > reading from a range of sources specificed in data_sources
    #o.add_readers_from_file(o.test_data_folder() + '../../opendrift/scripts/data_sources.txt')
    logging.debug('adding readers')
    from opendrift.readers import reader_basemap_landmask
    # Making customised landmask (Basemap)
    base_map_specs = config['model_frame']
     # add 'resolution'  and 'projection' items  to dictionary
    base_map_specs['resolution'] = config['basemap_resolution'] # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)
    base_map_specs['projection'] = 'merc' 
    reader_basemap = reader_basemap_landmask.Reader(
                llcrnrlon= base_map_specs['llcrnrlon'], llcrnrlat=base_map_specs['llcrnrlat'], 
                urcrnrlon= base_map_specs['urcrnrlon'], urcrnrlat=base_map_specs['urcrnrlat'],  
                resolution=base_map_specs['resolution'], projection='merc') 
    o.add_reader(reader_basemap)
    # use 'merc' projection by default - hardcoded for now
    # loop through different readers and initialize
    base_str = 'from opendrift.readers import '
    for key in config['readers'].keys():
        if config['readers'][key] is not None:
            if 'filename' in config['readers'][key]: # use local file(s)
                for read_cnt,fname in enumerate(config['readers'][key]['filename']):
                    # import correct reader class
                    reader_type = config['readers'][key]['reader_type'][read_cnt]
                    exec(base_str + reader_type + ' as reader_obj',globals(),locals()) # needed to add globals() to make it work
                    # load reader as reader_obj, then initialize it
                    if 'options' in config['readers'][key]: # add specific options as **kwargs
                        opts =  config['readers'][key]['options']
                    else: # no specific options
                        opts = {}

                    read_tmp = reader_obj.Reader(filename=fname,**opts)
                    #save reader with new name 
                    exec(key + '_' + str(read_cnt) + '= read_tmp',globals(),locals()) 
                    del read_tmp
                    # add reader to OpenDrift object
                    o.add_reader(eval( key + '_' + str(read_cnt)) )   
    
            # elif 'uds' in config['readers'][key]:
            # elif 'threds' in config['readers'][key]: 
            # etc...to define depending of needs
    # adding fallback values / constants
    if 'fallback_values' in config.keys():
        for variable in config['fallback_values'].keys():
            o.fallback_values[variable] = config['fallback_values'][variable]
    #############################################################################################
    # set all additional configuration - overwriting if needed
    #############################################################################################       
    logging.debug('setting config')
    for ic,category in enumerate(config['config']):
        for ip,param in enumerate(config['config'][category]):
            o.set_config('%s:%s' % (category,param), config['config'][category][param])
            logging.debug('CONFIG  %s:%s,%s' % (category,param,config['config'][category][param]))

    #############################################################################################
    # start run
    #############################################################################################
    logging.debug('setting model run')
    time_step = timedelta(seconds=config['time_step_sec']) 
    if 'time_step_sec_output' in config.keys():
        time_step_output = timedelta(seconds=config['time_step_sec_output'])
    else:
        time_step_output = timedelta(minutes=30)
    # duration of run can be specifed by duration or end_time_run
    if 'duration_hours' in config.keys():
        duration = timedelta(hours=config['duration_hours'])
        end_time_run = start_time + duration
    elif 'end_time_run' in config.keys():
        end_time_run = datetime.strptime(config['end_time_run'],'%d-%m-%Y %H:%M')      
    else: # if none input , just use end_time which is end of release
        end_time_run = end_time

    if config['run_backwards']: # run model backwards 
        time_step = -time_step
        # here we should have start_time>end_time
        if start_time < end_time_run:
            logging.error('start_time must be smaller than end_time if run backwards')

    if 'outfile' not in config.keys():
        outfile = 'opendrift_' + config['imp'] + '_' + start_time.strftime('%d%m%Y_%H%M')
    else:
        outfile = config['outfile']

    # run simulation
    o.run(stop_on_error = config['stop_on_error'],time_step=time_step,
          end_time = end_time_run, time_step_output = time_step_output,
          outfile=outfile)
          #could use **extra_args like other initialization if needed
    return o


if __name__ == '__main__':
    ###############################
    # read config file input as argument
    ###############################
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file',
                        help='<Configuration file (YAML-format) with all required infos to run Opendrift>')
    args = parser.parse_args()

    # setup opendrift object using specifications in YAML configuration file
    config = read_yaml_config(args.config_file) 
    # run simulation
    o = run_opendrift_from_config(config)
    # post-process/plots
    print o
    if config['post_process']['show_plot']:
        o.plot()
    if config['post_process']['show_animate']:
        o.animation()
    if config['post_process']['show_oil_budget']:
        o.plot_oil_budget()
    # if config['post_process']['save_animate']:
        # o.animation()
    # if config['post_process']['save_plot']:
        # o.plot()
    # if config['post_process']['save_oil_budget']:
         # o.plot_oil_budget()