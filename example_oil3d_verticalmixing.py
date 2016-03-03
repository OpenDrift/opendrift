#!/usr/bin/env python

import os.path
from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil3D import OpenOil3D


# Copy the following files (~5 GB) to your local computer for faster run:
# http://super-monitor.met.no/thredds/fileServer/lustreMntB/users/knutfd/public/AROME_MetCoOp_00_DEF.nc
# http://super-monitor.met.no/thredds/fileServer/lustreMntB/users/knutfd/public/subfolder_test/NorKyst-800m_ZDEPTHS_his_00.nc

o = OpenOil3D(loglevel=0)  # Set loglevel to 0 for debug information

ncfile = 'oil3Dmixing.nc'
import_file = False

arome_file = 'test_data/14Jan2016_NorKyst_z_3d/AROME_MetCoOp_00_DEF.nc_20160114_subset'
norkyst_file = 'test_data/14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc'
arome_file = 'gttg'
if not os.path.isfile(arome_file):
    arome_file = 'http://super-monitor.met.no/thredds/dodsC/lustreMntB/users/knutfd/public/AROME_MetCoOp_00_DEF.nc'
    norkyst_file = 'http://super-monitor.met.no/thredds/dodsC/lustreMntB/users/knutfd/public/subfolder_test/NorKyst-800m_ZDEPTHS_his_00.nc'

if import_file is True:
    o.io_import_file(ncfile)
else:
    reader_arome = reader_netCDF_CF_generic.Reader(arome_file)
    reader_norkyst = reader_netCDF_CF_generic.Reader(norkyst_file)

    # Landmask (Basemap)
    reader_basemap = reader_basemap_landmask.Reader(
                        llcrnrlon=3.0, llcrnrlat=61.4,
                        urcrnrlon=6.2, urcrnrlat=62.6,
                        resolution='h', projection='merc')

    o.add_reader([reader_norkyst, reader_basemap, reader_arome])

    # Seeding some particles
    lon = 4.9; lat = 62.1; # Stad

    #time = [reader_arome.start_time,
    #        reader_arome.start_time + timedelta(hours=30)]
    time = reader_arome.start_time

    # Seed oil elements at defined position and time
    o.seed_elements(lon, lat, z=-0.5, radius=0, number=3000, time=time,
                    density=880)

    # Adjusting some configuration
    o.config['processes']['turbulentmixing'] = True
    o.config['turbulentmixing']['diffusivitymodel'] = 'windspeed_Sundby1983'
    #o.config['turbulentmixing']['diffusivitymodel'] = 'stepfunction'
    o.config['turbulentmixing']['timestep'] = 2. # seconds

    # Running model (until end of driver data)
    o.run(end_time=reader_norkyst.start_time + timedelta(hours=48),
          time_step=900, time_step_output=timedelta(hours=1),
          outfile=ncfile)

###########################
# Print and plot results
###########################
print o

o.plot_oil_budget()
o.plot()
#o.plot(linecolor='z')
o.animation()
o.plot_vertical_distribution()

