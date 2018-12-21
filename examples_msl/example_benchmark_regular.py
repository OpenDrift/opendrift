#!/usr/bin/env python

from opendrift.readers import reader_netCDF_MetOcean, reader_basemap_landmask
from opendrift.models.oceandrift3D import OceanDrift3D
import numpy as np
from datetime import datetime, timedelta

o = OceanDrift3D(loglevel=0)  # Set loglevel to 0 for debug information

# Making customised landmask (Basemap)  
reader_basemap = reader_basemap_landmask.Reader(
                llcrnrlon= 168.0199, llcrnrlat=-42.8449, 
                urcrnrlon= 177.4601, urcrnrlat=-37.9051,   
                resolution='i', projection='merc') # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)

reader_struct2D = reader_netCDF_MetOcean.Reader('F:/metocean/R&D_OpenDrift/benchmark_runs/cnz20040101_00z_surf.nc')
# reader_struct3D = reader_netCDF_CF_generic.Reader('F:/metocean/0387_OMV_drillcuttings/opendrift_modelling/flow_fields/cnz20040101_00z_3D.nc')
o.add_reader([reader_basemap,reader_struct2D])


# Seeding some particles
lon = 174.1656; lat = -40.3346; # estuary entrance
nb_parts = 5000 
start_date = None
start_date = reader_struct2D.start_time
# Seed elements at defined position and time
z = -np.random.rand(nb_parts)*2  # Giving elements a random depth - negative down convention
o.seed_elements(lon, lat, z=z, 
	                        radius=0,
	                        number=nb_parts, 
	                        time= [start_date,start_date+timedelta(days=4.0)]) # time=start_date
                          # this will spread the 'number' of particle release,
                          # evenly over the period [tstart_release tend_release]

# Adjusting some configuration
# o.list_config()
# o.list_configspec()
# o.fallback_values['ocean_vertical_diffusivity'] = 0.0001
o.set_config('drift:scheme','runge-kutta4') # 'euler' or 'runge-kutta'
o.set_config('processes:turbulentmixing', False)
o.set_config('processes:verticaladvection' , False) # no vertical current available, so no vertical advection 
o.set_config('turbulentmixing:diffusivitymodel', 'environment') # i.e. specified from model or constant
o.set_config('turbulentmixing:TSprofiles',False)
o.set_config('turbulentmixing:timestep', 90)  # sub-timestep for vertical mixing
o.set_config('drift:current_uncertainty', 0.05)
# o.set_config('drift:wind_uncertainty', 2)
o.set_config('general:coastline_action','stranding') # option('none', 'stranding', 'previous', default='stranding')

# Running model
o.run(stop_on_error = True,time_step=600,
      end_time = start_date + timedelta(days = 5.00), 
       outfile='example_struct.nc',time_step_output = 600.0)
# the start time is defined by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime

# Print and plot results
print(o)
o.plot(linecolor='z')  # Color lines according to depth
# o.animation()
