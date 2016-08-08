#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(
    'test_data/16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(
    'test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=4, llcrnrlat=59.7,
                    urcrnrlon=7, urcrnrlat=61.5,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_arome])

# Seeding some particles
lon = 4.9; lat = 60.0; # Outside Bergen

#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=30)]
time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=3000, number=2000, time=time,
                oiltype='GULLFAKS AB')

# Adjusting some configuration
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = True
o.config['processes']['evaporation'] = True
o.config['processes']['emulsification'] = True

# Running model (until end of driver data)
o.run(steps=66*4, time_step=900, outfile='openoil.nc')

# Print and plot results
print o

o.plot_oil_budget()
o.plot_environment()

o.plot()
o.animation()
