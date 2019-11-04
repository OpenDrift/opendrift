"""
Leeway
===============================

asdf
"""

################################################
# Testing some stuff1
# --------------------
#
# Here's the animation..
#

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway
import matplotlib.pyplot as plt

lw = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(lw.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader(lw.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

################################################
# does this work?

# Adding readers succesively, and specifying which variables they
# shall provide. This way, order of adding readers does not matter,
# except for small rounding differences due to different projection
lw.add_reader(reader_norkyst,
              variables=['x_sea_water_velocity', 'y_sea_water_velocity'])
lw.add_reader(reader_arome,
              variables=['x_wind', 'y_wind'])
lw.fallback_values['x_sea_water_velocity'] = 0
lw.fallback_values['y_sea_water_velocity'] = 0

# Seeding some particles
lon = 4.5; lat = 60.0; # Outside Bergen

# Seed leeway elements at defined position and time
objType = 26  # 26 = Life-raft, no ballast
lw.seed_elements(lon, lat, radius=1000, number=3000,
                 time=reader_arome.start_time, objectType=objType)

lw.set_projection('+proj=merc')

# Running model (until end of driver data)
lw.run(steps=60, time_step=900)
lw.animation(filename='plot_leeway.gif')
plt.close('all')

################################################
# Testing some stuff
# ------------------
#
# Here's the animation..
#
# .. image:: /gallery/plot_leeway.gif

# Print and plot results
print(lw)
#lw.animation(filename='leeway.gif')
lw.plot()
