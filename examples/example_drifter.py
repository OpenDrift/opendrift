#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.oceandrift import OceanDrift

o = OceanDrift()  # Basic drift model suitable for passive tracers or drifters

#######################
# Preparing Readers
#######################

reader_current = reader_netCDF_CF_generic.Reader('../test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
reader_wind = reader_netCDF_CF_generic.Reader('../test_data/16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=4, llcrnrlat=59.8,
                    urcrnrlon=6, urcrnrlat=61,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_current, reader_wind])

#######################
# Seeding elements
#######################

# Elements are moved with the ocean current, in addition to a fraction of 
# the wind speed (wind_drift_factor). This factor depends on the properties
# of the elements. Typical empirical values are:
# - 0.035 (3.5 %) for oil and iSphere driftes
# - 0.01  (1 %) for CODE drifters partly submerged ~0.5 m
# As there are large uncertainties, it makes sence to provide a statistical
# distribution of wind_drift_factors

# Using a constant value for all elements:
#wind_drift_factor = 0.03

# Giving each element a unique (random) wind_drift_factor
wind_drift_factor = np.random.uniform(0, 0.06, 2000)

o.seed_elements(4.7, 59.9, radius=3000, number=2000,
                time=reader_current.start_time,
                wind_drift_factor=wind_drift_factor)

#######################
# Running model
#######################
o.run(time_step=timedelta(minutes=15),
      time_step_output=timedelta(minutes=60))

###########################
# Print and plot results
###########################
print o
o.animation()

# Plot trajectories, colored by the wind_drift_factor of each element
o.plot(linecolor='wind_drift_factor')
