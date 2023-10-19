#!/usr/bin/env python
"""
Model landmask
===============================

Comparing two simulation runs, with landmask from ocean model and GSHHG
"""

from datetime import timedelta

from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

lon = 14.75; lat = 68.1

o = OceanDrift(loglevel=20)

reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

#%%
# First run, with default GSHHG vector landmask
o.add_reader([reader_nordic])
time = reader_nordic.start_time
o.set_config('general:use_auto_landmask', True)
o.seed_elements(lon, lat, radius=3000, number=1000, time=time)
o.run(end_time=reader_nordic.end_time, time_step=1800, time_step_output=3*3600)

#%%
# Second run, with landmask from ocean model
o2 = OceanDrift(loglevel=20)
o2.add_reader([reader_nordic])
lon = 14.75; lat = 68.1
o2.set_config('general:use_auto_landmask', False)
o2.seed_elements(lon, lat, radius=3000, number=1000, time=time)
o2.run(end_time=reader_nordic.end_time, time_step=1800, time_step_output=3*3600)

#%% Prepare cusom colormap/colors for land and ocean
from matplotlib.colors import ListedColormap
import cartopy.feature as cfeature
cmap = ListedColormap([cfeature.COLORS['water'],
                       cfeature.COLORS['land']])

#%%
# .. _model_landmask_only_model:
#
# To only show the landmask from the model, hide the coastline landmask by doing:

o2.plot(background='land_binary_mask', hide_landmask=True, cmap=cmap)


#%%
# Animation illustrating that red particles strand at ocean model land cells, and black particles strand at GSHHG land polygons
o.animation(compare=o2, background='land_binary_mask', cmap=cmap,
            legend=['Default GSHHG landmask', 'Ocean model landmask'])

#%%
# .. image:: /gallery/animations/example_model_landmask_0.gif

o.plot(compare=o2, background='land_binary_mask', cmap=cmap,
       legend=['Default GSHHG landmask', 'Ocean model landmask'])

