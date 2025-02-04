#!/usr/bin/env python
"""
River runoff
===========================
"""

import os
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import opendrift
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_oscillating


outfile = 'runoff.nc'  # Raw simulation output
histogram_file = 'runoff_histogram.nc'

#%%
# First make a simulation with two seedings, marked by *origin_marker*
o = OceanDrift(loglevel=20)
o.set_config('drift:horizontal_diffusivity', 300)
o.set_config('general:coastline_action', 'previous')
t1 = datetime.now()
t2 = t1 + timedelta(hours=48)
reader_x = reader_oscillating.Reader('x_sea_water_velocity', period_seconds=3600*24,
                amplitude=1, zero_time=t1)
reader_y = reader_oscillating.Reader('y_sea_water_velocity', period_seconds=3600*72,
                amplitude=.5, zero_time=t2)
o.add_reader([reader_x, reader_y])
number = 25000
o.seed_elements(time=[t1, t2], lon=9.017931, lat=58.562702, number=number,
                origin_marker_name='River 1')  # River 1
o.seed_elements(time=[t1, t2], lon=8.824815, lat=58.425648, number=number,
                origin_marker_name='River 2')  # River 2
seed_times = o.elements_scheduled_time[0:number]

o.run(duration=timedelta(hours=48),
      time_step=1800, time_step_output=3600, outfile=outfile)

#%%
# Opening the output file lazily with Xarray.
# This will work even if the file is too large to fit in memory, as it
# will read and process data chuck-by-chunk directly from file using Dask.
o = opendrift.open(outfile)

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
river_water = o.get_histogram(pixelsize_m=1500, weights=runoff, density=False)

rw = river_water.sum(dim='origin_marker')  # For both rivers
#rw = river_water.isel(origin_marker=1)    # For one of the rivers
river_water.name = 'River water [m3/cell]'

text = [{'s': o.origin_marker[0], 'x': 8.55, 'y': 58.56, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
        {'s': o.origin_marker[1], 'x': 8.35, 'y': 58.42, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
        {'s': '* Station', 'x': station_lon, 'y': station_lat, 'fontsize': 20, 'color': 'k',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='none', edgecolor='none', alpha=0.4), 'zorder': 1000}]
box = [{'lon': box1_lon, 'lat': box1_lat, 'text': 'Area 1', 'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'},
       {'lon': box2_lon, 'lat': box2_lat, 'text': 'Area 2', 'fc': 'none', 'alpha': 0.8, 'lw': 1, 'ec': 'k'}]

o.animation(background=rw.where(rw>0), bgalpha=1, fast=False,
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
# Area 1
t1 = river_water.sel(lon_bin=slice(box1_lon[0], box1_lon[1]),
                     lat_bin=slice(box1_lat[0], box1_lat[1]))
t1 = t1.sum(('lon_bin', 'lat_bin'))
t1.isel(origin_marker=0).plot(label=o.origin_marker[0], ax=ax2)
t1.isel(origin_marker=1).plot(label=o.origin_marker[1], ax=ax2)
t1.sum(dim='origin_marker').plot(label='Total', linestyle='--', ax=ax2)
ax2.set_title('Amount of water passing through Area 1')
# Area 2
t2 = river_water.sel(lon_bin=slice(box2_lon[0], box2_lon[1]),
                     lat_bin=slice(box2_lat[0], box2_lat[1]))
t2 = t2.sum(('lon_bin', 'lat_bin'))
t2.isel(origin_marker=0).plot(label=o.origin_marker[0], ax=ax3)
t2.isel(origin_marker=1).plot(label=o.origin_marker[1], ax=ax3)
t2.sum(dim='origin_marker').plot(label='Total', linestyle='--', ax=ax3)
ax3.set_title('Amount of water passing through Area 2')
# Extracting time series at the location of the station
t = river_water.sel(lon_bin=station_lon, lat_bin=station_lat, method='nearest')
t.isel(origin_marker=0).plot(label=o.origin_marker[0], ax=ax4)
t.isel(origin_marker=1).plot(label=o.origin_marker[1], ax=ax4)
t.sum(dim='origin_marker').plot(label='Total', linestyle='--', ax=ax4)
ax4.legend()
ax4.margins(x=0)

for ax in [ax1, ax2, ax3]:
    ax.margins(x=0)
    ax.legend()
    #ax.set_xticks([])
    ax.set_xlabel(None)
ax4.set_title('Density of water at Station')
# TODO disabling due to recent problem with dateformatter
#ax4.xaxis.set_major_formatter(DateFormatter("%d %b %H"))
plt.show()

#%%
# Finally, plot the spatial distribution of mean age of water 
num = o.get_histogram(pixelsize_m=1500, weights=None, density=False)
num.name = 'number'
num.to_netcdf(histogram_file)

mas = o.get_histogram(pixelsize_m=1500, weights=o.result.age_seconds, density=False)
mas = mas/3600  # in hours
mas = mas/num  # per area
mas.name='mean_age'
mas.to_netcdf(histogram_file, 'a')
mas = mas.mean(dim='time').sum(dim='origin_marker')  # Mean time of both rivers
#mas = mas.mean(dim='time').isel(origin_marker=1)  # Mean age of a single river
mas.name='Mean age of water [hours]'

o.plot(background=mas.where(mas>0), fast=True, show_elements=False, show_trajectories=False)


# Cleaning up
os.remove(outfile)
os.remove(histogram_file)
