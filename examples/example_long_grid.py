"""
Grid
=============
"""
#!/usr/bin/env python

import numpy as np

from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Making customised landmask
reader_landmask = reader_global_landmask.Reader(
                    llcrnrlon=3.5, llcrnrlat=59.9,
                    urcrnrlon=5.5, urcrnrlat=61.2)

o.add_reader([reader_landmask, reader_norkyst])

# Seeding some particles
lons = np.linspace(3.5, 5.0, 100)
lats = np.linspace(60, 61, 100)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()

# Seed oil elements at defined position and time
o.seed_elements(lons, lats, radius=0, number=10000,
                time=reader_norkyst.start_time)

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800)

# Print and plot results
print(o)
o.animation(filename='grid.gif')

#%%
# .. image:: https://camo.githubusercontent.com/aaf79a5e243ef601c015c290fc1f935d2ec70340/68747470733a2f2f646c2e64726f70626f7875736572636f6e74656e742e636f6d2f732f337933316c64773235636c6e7362392f6f70656e64726966745f677269642e6769663f646c3d30
#o.plot()
# .. image:: /gallery/animations/grid.gif
