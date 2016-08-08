#!/usr/bin/env python

import os.path
from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil3D import OpenOil3D


o = OpenOil3D(loglevel=0)  # Set loglevel to 0 for debug information

ncfile = 'oil3Dmixing.nc'
import_file = False  # Set to True to import previous run

if import_file is True:
    o.io_import_file(ncfile)
else:
    reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
    reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

    # Landmask (Basemap)
    reader_basemap = reader_basemap_landmask.Reader(
                        llcrnrlon=3.0, llcrnrlat=61.4,
                        urcrnrlon=6.2, urcrnrlat=62.6,
                        resolution='h', projection='merc')

    o.add_reader([reader_norkyst, reader_basemap, reader_arome])

    # Seeding some particles
    lon = 4.9; lat = 62.1; # Stad

    time = reader_arome.start_time

    # Seed oil elements at defined position and time
    o.seed_elements(lon, lat, z=-0.5, radius=1000, number=2000, time=time,
                    density=880)

    # Adjusting some configuration
    o.config['processes']['turbulentmixing'] = True
    o.config['turbulentmixing']['diffusivitymodel'] = 'windspeed_Sundby1983'
    #o.config['turbulentmixing']['diffusivitymodel'] = 'stepfunction'
    o.config['turbulentmixing']['timestep'] = 2. # seconds

    # Running model (until end of driver data)
    o.run(end_time=reader_arome.start_time + timedelta(hours=12),
          time_step=900, time_step_output=1800, outfile=ncfile)

###########################
# Print and plot results
###########################
print o

o.plot(linecolor='z')
o.animation()
o.plot_vertical_distribution()  # Note the interactive slider at the bottom

