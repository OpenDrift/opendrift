#!/usr/bin/env python
"""
Vertical advection correction
=============================
"""

#%%
# In OpenDrift, element depth (z) is relative to actual surface, and not mean sea level.
# Vertical velocity from ocean models contain a contribution which is due to
# change in sea surface elevation, however, this should in OpenDrift be subtracted
# due to the choice of defining z relative to actual surface elevation
# This correction is activated with config setting 'drift:vertical_advection_correction'
# The effect is illustrated with the simulations below.


from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift

time = datetime.now()
reader_tidal = reader_oscillating.Reader('sea_surface_height', amplitude=-1,
                                         period=timedelta(hours=6), zero_time=time)
lat=59.8113; lon=10.5517  # Oslo fjord
z = np.arange(0, -60, -5)  # Seeding one particle every 5 meter from surface to 60m depth

# Without stetching
o = OceanDrift(loglevel=0)
o.add_reader(reader_tidal)
o.add_readers_from_list(['https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'])
o.set_config('drift:water_column_stretching', False)
o.set_config('drift:vertical_advection', True)
o.seed_elements(lon=lon, lat=lat, time=time, z=z, number=len(z))
o.run(duration=timedelta(hours=24), time_step=1800)

# With stetching
o2 = OceanDrift(loglevel=0)
o2.add_reader(reader_tidal)
o2.add_readers_from_list(['https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'])
o2.set_config('drift:water_column_stretching', True)
o2.set_config('drift:vertical_advection', True)
o2.seed_elements(lon=lon, lat=lat, time=time, z=z, number=len(z))
o2.run(duration=timedelta(hours=24), time_step=1800)

# With w correction
o3 = OceanDrift(loglevel=0)
o3.add_reader(reader_tidal)
o3.add_readers_from_list(['https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'])
o3.set_config('drift:water_column_stretching', False)
o3.set_config('drift:vertical_advection', True)
o3.set_config('drift:vertical_advection_correction', True)
o3.seed_elements(lon=lon, lat=lat, time=time, z=z, number=len(z))
o3.run(duration=timedelta(hours=24), time_step=1800)


fig = plt.figure(figsize=(10, 10))
gs = fig.add_gridspec(5, 1)
ax1 = fig.add_subplot(gs[0:3, 0])
o.result.z.plot.line(ax=ax1, x='time', add_legend=False, color='k')
o2.result.z.plot.line(ax=ax1, x='time', add_legend=False, color='r')
o3.result.z.plot.line(ax=ax1, x='time', add_legend=False, color='g')
plt.plot([], [], color='k', label='No stretching or correction')
plt.plot([], [], color='r', label='Column stretching')
plt.plot([], [], color='g', label='Vertical velocity correction')
plt.legend()
ax2 = fig.add_subplot(gs[3, 0])
o.result.sea_floor_depth_below_sea_level.plot.line(ax=ax2, x='time', add_legend=False, color='k')
ax3 = fig.add_subplot(gs[4, 0])
o.result.upward_sea_water_velocity.plot.line(ax=ax3, x='time', add_legend=False, color='k')
plt.show()
