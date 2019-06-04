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


def read_yaml_config(config_yaml_file):
    ''' read a YAML config file for OpenDrift - could be moved to opendrift/__init__.py eventually '''
    try:
        config = yaml.load(open(config_yaml_file).read())
        return config
    except:
        sys.exit('Cannot read ' + config_yaml_file)

def run_opendrift_from_config(config):
    ''' Setup and run opendrit simulation based on input configuration
        General process is inspired from run_opendrift() from opendrift_gui.py
    '''

    # sys.stdout.write('running OpenDrift')
    # release position and timing
    start_time = datetime.strptime(config['start_time'],'%d-%m-%Y %H:%M')
    end_time = datetime.strptime(config['end_time'],'%d-%m-%Y %H:%M')
    sys.stdout.flush()
    lon = config['position']['lon']
    lat = config['position']['lon']
    radius = config['position']['lon']
    z = config['position']['z']
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
    else:
        cone = False
    # we could do something on the z here if required
    # e.g. random within range etc..
   
    # model selection
    if config['model'] == 'Leeway':
        from opendrift.models.leeway import Leeway
        o = Leeway(loglevel=0)
        # add extra model arguments if relevant
        #extra_seed_args = {'objectType': ln + 1}
    elif config['model'] == 'OceanDrift3D':
        from opendrift.models.oceandrift3D import OceanDrift3D
        o = OceanDrift3D(loglevel=0)
    elif config['model'] == 'OpenOil3D':
        from opendrift.models.openoil3D import OpenOil3D
        o = OpenOil3D(weathering_model='noaa', loglevel=0)
        # use 'noaa' model by default  for now
        # add extra model argument if relevant
        # if extra_model_args:
        #       read which weathering_model is input 'default' or 'noaa'
        #       and initialize model accordingly
    elif config['model'] == 'ShipDrift':
        from opendrift.models.shipdrift import ShipDrift
        o = ShipDrift(loglevel=0)

    import pdb;pdb.set_trace()

    # particle seeding
     >>> to do !
    parse extra_seed_args to initialize seeding correctly
    note these can still be edited further down in config

    extra_seed_args = {}
    for se in self.seed_input:
        if se == 'ocean_only':
            continue  # To be fixed/removed
        val = self.seed_input[se].get()
        try:
            extra_seed_args[se] = np.float(val)
        except:
            extra_seed_args[se] = val
        self.o.set_config('seed:' + se, val)
    if 'object_type' in extra_seed_args:
        extra_seed_args['objectType'] = extra_seed_args['object_type']
        del extra_seed_args['object_type']


    extra_seed_args = {}
    self.o.seed_elements(lon=lon, lat=lat, number=5000, radius=radius,
                    time=start_time, cone=cone,
                    **extra_seed_args)

    self.o.add_readers_from_file(self.o.test_data_folder() +
        '../../opendrift/scripts/data_sources.txt')

    #time_step = o.get_config('general:time_step_minutes')*60
    time_step = 900  # 15 minutes
    time_step_output = timedelta(minutes=30)
    duration = int(self.durationhours.get())*3600/time_step
    if self.directionvar.get() == 'backwards':
        time_step = -time_step
    if self.has_diana is True:
        extra_args = {'outfile': self.outputdir + '/opendrift_' +
            self.model.get() + self.o.start_time.strftime('_%Y%m%d_%H%M.nc')}
    else:
        extra_args = {}

    mapres = self.mapresvar.get()[0]
    self.o.set_config('general:basemap_resolution', mapres)         
    self.simulationname = 'opendrift_' + self.model.get() + \
        self.o.start_time.strftime('_%Y%m%d_%H%M')

    # Starting simulation run
    self.o.run(steps=duration, time_step=time_step,
          time_step_output=time_step_output, **extra_args)
    print(self.o)

    self.results.destroy()
    self.results = tk.Frame(self.seed, bd=2,
                           relief=tk.FLAT, padx=5, pady=0)
    self.results.grid(row=7, column=3, columnspan=1, sticky='ew')
    tk.Button(self.results, text='Show animation',
              command=lambda: self.handle_result(
                'showanimation')).grid(row=1, column=1)
    tk.Button(self.results, text='Save animation',
              command=lambda: self.handle_result(
                'saveanimation')).grid(row=2, column=1)
    tk.Button(self.results, text='Show plot',
              command=lambda: self.handle_result(
                'showplot')).grid(row=3, column=1)
    tk.Button(self.results, text='Save plot',
              command=lambda: self.handle_result(
                'saveplot')).grid(row=4, column=1)
    if self.model.get() == 'OpenOil':
        tk.Button(self.results, text='Save oil budget',
                  command=lambda: self.handle_result(
                    'saveoilbudget')).grid(row=5, column=1)
        tk.Button(self.results, text='Show oil budget',
                  command=lambda: self.handle_result(
                    'showoilbudget')).grid(row=6, column=1)

    if self.has_diana is True:
        diana_filename = self.dianadir + self.simulationname + '.nc'
        self.o.write_netcdf_density_map(diana_filename)
        tk.Button(self.results, text='Show in Diana',
                  command=lambda: os.system('diana &')
                  ).grid(row=8, column=1)


if __name__ == '__main__':
    ###############################
    # read config file input as argument
    ###############################
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file',
                        help='<Configuration file (YAML-format) with all required infos to run Opendrift>')
    args = parser.parse_args()
    # setup opendrift object using specifications in YAML configuration file
    config = read_yaml_config(args.config_file) #
    run_opendrift_from_config(config)
    

    ###############################
    # RUN 
    ###############################
    # Running model
    # here we could re-use  run_opendrift() from opendrift_gui.py. 
    # Instead of using the specs from the GUI, we'd use the specs from the config file
    run_opendrift(config)


    o.plot_oil_budget()

###############################
# PLOTS / ANIMATION
###############################

    # Print and plot results
    print o
    o.animation()
    o.plot()



# def run_opendrift_from_config(config):
#         # o.run(time_step = 1800, 
#         #   end_time = reader_roms_cnz_surface.start_time + timedelta(days = 7.0),
#         #   outfile='opendrift_oil3d_all_inputs.nc',
#         #   time_step_output = 1800)