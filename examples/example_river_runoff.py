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
                origin_marker=0)  # River 1
o.seed_elements(time=[t1, t2], lon=8.824815, lat=58.425648, number=number,
                origin_marker=1)  # River 2
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
# We want to extract timeseries at the coordinates of a hypothetical measuring station
station_lon = 9.4
station_lat = 58.1

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
plt.plot(seed_times, runoff_river1, 'b', label='River 1')
plt.plot(seed_times, runoff_river2, 'r', label='River 2')
plt.ylabel('Runoff  [m3/s]')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d %b %H:%M'))
plt.margins(x=0)
plt.legend()
plt.show()

runoff = np.concatenate((runoff_river1, runoff_river2))
text = [{'s': 'River 1', 'x': 8.55, 'y': 58.56, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
        {'s': 'River 2', 'x': 8.35, 'y': 58.42, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000},
        {'s': '* Station', 'x': station_lon, 'y': station_lat, 'fontsize': 20, 'color': 'k',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='none', edgecolor='none', alpha=0.4), 'zorder': 1000}]
o.animation(density=runoff, density_pixelsize_m=1500, fast=True,
            show_elements=False, vmin=0, vmax=120, text=text)

#%%
# Extracting time series at the location of the station
t, t_om = o.get_density_timeseries(lon=station_lon, lat=station_lat)
t_om.isel(origin_marker=0).plot(label='River 1')
t_om.isel(origin_marker=1).plot(label='River 2')
t.plot(label='Total', linestyle='--')
plt.legend()
plt.title('Density of water at Station') 
plt.show()

#%%
# .. image:: /gallery/animations/example_river_runoff_0.gif

# Cleaning up
os.remove(outfile)
os.remove(analysis_file)
