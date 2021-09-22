#!/usr/bin/env python
"""
River runoff
===========================
"""

import os
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import opendrift
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_oscillating


outfile = 'runoff.nc'  # Raw simulation output
analysis_file = 'runoff_density.nc'  # Raw simulation output
try:
    os.remove(analysis_file)
except OSError:
    pass

#%%
# First make a simulation with two seedings, marked by *origin_marker*
o = OceanDrift(loglevel=20)
t1 = datetime.now()
t2 = t1 + timedelta(hours=48)
number = 25000
o.seed_elements(time=[t1, t2], lon=9.017931, lat=58.562702, number=number,
                origin_marker_name='River 1')  # River 1
o.seed_elements(time=[t1, t2], lon=8.824815, lat=58.425648, number=number,
                origin_marker_name='River 2')  # River 2
seed_times = o.elements_scheduled_time[0:number]

reader_x = reader_oscillating.Reader('x_sea_water_velocity', period_seconds=3600*24,
                amplitude=1, zero_time=t1)
reader_y = reader_oscillating.Reader('y_sea_water_velocity', period_seconds=3600*72,
                amplitude=.5, zero_time=t2)
o.add_reader([reader_x, reader_y])
o.set_config('drift:horizontal_diffusivity', 300)
o.set_config('general:coastline_action', 'previous')
o.run(duration=timedelta(hours=48),
      time_step=1800, time_step_output=3600, outfile=outfile)

#%%
# Opening the output file lazily with Xarray.
# This will work even if the file is too large to fit in memory, as it
# will read and process data chuck-by-chunk directly from file using Dask.
# Note that the analysis file will be re-used if existing. Thus this file should be deleted after making any changes to the simulation above.
o = opendrift.open_xarray(outfile, analysis_file=analysis_file)

#%%
# We want to extract timeseries of river water at the coordinates of a hypothetical measuring station
# as well as the amount of river water passing through two defined areas/regions
station_lon = 9.4
station_lat = 58.1
box1_lon = [8.4, 8.8]
box1_lat = [57.9, 58.1]
box2_lon = [9.5, 9.9]
box2_lat = [58.3, 58.5]

#%%
# Animation of the spatial density of river runoff water.
# Although there are the same number of elements from each river, the density plots are
# weighted with the actual runoff at time of seeding. This weighting can be done/changed
# afterwards without needing to redo the simulation.
# The calculated density fields will be stored/cached in the analysis file
# for later re-use, as their calculation may be time consuming
# for huge output files.
# Note that other analysis/plotting methods are not yet adapted
# to datasets opened lazily with open_xarray

runoff_river1 = np.abs(np.cos(np.arange(number)*2*np.pi/(number)))  # Impose a temporal variation of runoff
runoff_river2 = 10*runoff_river1  # Let river 2 have 10 times as large runoff as river 1
runoff = np.concatenate((runoff_river1, runoff_river2))

#%%
# Calculate density with given pixel size, weighted by runoff amount per element
o.get_density_xarray(pixelsize_m=1500, weights=runoff)


text = [{'s': o.origin_marker[0], 'x': 8.55, 'y': 58.56, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
        {'s': o.origin_marker[1], 'x': 8.35, 'y': 58.42, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
        {'s': '* Station', 'x': station_lon, 'y': station_lat, 'fontsize': 20, 'color': 'k',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='none', edgecolor='none', alpha=0.4), 'zorder': 1000}]
box = [{'lon': box1_lon, 'lat': box1_lat, 'text': 'Area 1', 'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
       {'lon': box2_lon, 'lat': box2_lat, 'text': 'Area 2', 'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'}]

o.animation(background=o.ads.density.where(o.ads.density>0), bgalpha=1, fast=False,
            show_elements=False, vmin=0, vmax=120, text=text, box=box)

#%%
# .. image:: /gallery/animations/example_river_runoff_0.gif

#%%
# Plotting time series of river runoff, and corresponding water passing through the station and the two defined areas/boxes
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
# Runoff
ax1.plot(seed_times, runoff_river1, label=o.origin_marker[0])
ax1.plot(seed_times, runoff_river2, label=o.origin_marker[1])
ax1.set_ylabel('Runoff  [m3/s]')
ax1.set_title('Runoff')
ax1.margins(x=0)
ax1.legend()
# Area 1
t1, t1_om = o.get_density_timeseries(lon=box1_lon, lat=box1_lat)
t1_om.isel(origin_marker=0).plot(label=o.origin_marker[0], ax=ax2)
t1_om.isel(origin_marker=1).plot(label=o.origin_marker[1], ax=ax2)
t1.plot(label='Total', linestyle='--', ax=ax2)
ax2.legend()
ax2.margins(x=0)
ax2.set_title('Amount of water passing through Area 1')
# Area 2
t2, t2_om = o.get_density_timeseries(lon=box2_lon, lat=box2_lat)
t2_om.isel(origin_marker=0).plot(label=o.origin_marker[0], ax=ax3)
t2_om.isel(origin_marker=1).plot(label=o.origin_marker[1], ax=ax3)
t2.plot(label='Total', linestyle='--', ax=ax3)
ax3.legend()
ax3.margins(x=0)
ax3.set_title('Amount of water passing through Area 2')
# Extracting time series at the location of the station
t, t_om = o.get_density_timeseries(lon=station_lon, lat=station_lat)
t_om.isel(origin_marker=0).plot(label=o.origin_marker[0], ax=ax4)
t_om.isel(origin_marker=1).plot(label=o.origin_marker[1], ax=ax4)
t.plot(label='Total', linestyle='--', ax=ax4)
ax4.legend()
ax4.margins(x=0)
ax4.set_title('Density of water at Station')

plt.show()


# Cleaning up
os.remove(outfile)
os.remove(analysis_file)
