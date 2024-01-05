#!/usr/bin/env python
"""
Macondo
==================================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

#%%
# This is a highly simplified 2D simulation, for illustration purpose only

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# For this datasource which does not contain standard_name, we impose a variable mapping
reader_globcurrent = None

try:
    reader_globcurrent = reader_netCDF_CF_generic.Reader(
        'https://tds0.ifremer.fr/thredds/dodsC/GLOBCURRENT-L4-CUREUL_HS-ALT_SUM-V03.0',
        standard_name_mapping={'eastward_eulerian_current_velocity': 'x_sea_water_velocity',
                               'northward_eulerian_current_velocity': 'y_sea_water_velocity'})
except:
    print('Thredds server not available, cannot run example')

if reader_globcurrent is not None:
    reader_oceanwind = reader_netCDF_CF_generic.Reader('https://tds0.ifremer.fr/thredds/dodsC/CERSAT-GLO-CLIM_WIND_L4-OBS_FULL_TIME_SERIE')

    # Add readers
    o.add_reader([reader_globcurrent, reader_oceanwind])

    # Seed some particles
    lon = -88.387161; lat = 28.736669  # Macondo location
    starttime = datetime(2010, 4, 21, 6, 0, 0)  # 4 hours after explosion
    time = [starttime, starttime + timedelta(hours=24*30)]
    o.seed_elements(lon, lat, radius=0, number=15000, time=time, oil_type='LIGHT LOUISIANNA SWEET, BP')

    # Run model
    print(o)
    o.run(duration=timedelta(days=30),
        time_step=timedelta(hours=3),
        time_step_output=timedelta(hours=6))

    # Print and plot results
    print(o)
    o.plot(fast=True)
    o.plot_oil_budget()
    o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_macondo_0.gif
