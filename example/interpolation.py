#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(
    'test_data/16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(
    'test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=3, llcrnrlat=59.5,
                    urcrnrlon=7, urcrnrlat=62,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_arome])

# Seeding some particles
lon = 4.2; lat = 60.0; # Outside Bergen
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=50, number=5000, time=time)


# Adjusting some configuration
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = False
o.config['processes']['evaporation'] = False
o.config['processes']['emulsification'] = False
o.config['drift']['current_uncertainty'] = .1
o.config['drift']['wind_uncertainty'] = 2
o.config['drift']['relative_wind'] = False

for r in o.readers:
    o.readers[r].interpolation = 'ndimage'

# Running model
o.run(steps=66*2, time_step=1800)

# Second run, for comparison
o2 = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_basemap, reader_norkyst, reader_arome])
o2.seed_elements(lon, lat, radius=50, number=5000, time=time)
o2.config['processes']['diffusion'] = True
o2.config['processes']['dispersion'] = False
o2.config['processes']['evaporation'] = False
o2.config['processes']['emulsification'] = False
o2.config['drift']['current_uncertainty'] = .1
o2.config['drift']['wind_uncertainty'] = 2
o2.config['drift']['relative_wind'] = False

for r in o.readers:
    o.readers[r].interpolation = 'linearND'

o2.run(steps=66*2, time_step=1800)


# Animate and compare the two runs
o.animation(compare=o2, legend=['ndimage', 'linearND'])
            
