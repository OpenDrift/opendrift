#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.pelagicegg import PelagicEggDrift

o = PelagicEggDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('test_data/arome_subset_16Nov2015.nc')
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('test_data/norkyst800_subset_16Nov2015.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=3, llcrnrlat=59.7,
                    urcrnrlon=7, urcrnrlat=61.5,
                    resolution='h', projection='merc')

o.add_reader([reader_norkyst, reader_basemap, reader_arome])

# Experiment with constant wind instead of arome
#o.add_reader([reader_norkyst, reader_basemap])
#o.fallback_values['x_wind'] = 20
#o.fallback_values['y_wind'] = 20

# Seeding some particles
lon = 4.5; lat = 60.0; # Outside Bergen

#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=30)]
time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, z=-0.5, radius=3000, number=1000, time=time)

# Adjusting some configuration
o.config['drift']['wind_drift_factor'] = .0

o.config['processes']['turbulentmixing'] = True
o.config['turbulentmixing']['diffusivitymodel'] = 'windspeed_Sundby1983'
#o.config['turbulentmixing']['diffusivitymodel'] = 'stepfunction'
o.config['turbulentmixing']['timestep'] = 2. # seconds

o.config['processes']['diffusion'] = False
o.config['drift']['current_uncertainty'] = .1
o.config['drift']['wind_uncertainty'] = 2

# Running model (until end of driver data)
o.run(steps=66*4, time_step=900)

# Print and plot results
print o

#o.plot()
#o.animation()
#o.plot_vertical_distribution()
o.plot_vertical_distribution_slider()

