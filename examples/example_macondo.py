#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

# This is a highly simplified 2D simulation, for illustration purpose only

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

#reader_hycom = reader_netCDF_CF_generic.Reader('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/2010/3hrly')
#print(reader_hycom0)
reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V02.0_FULL_TIME_SERIE')  # Total

reader_oceanwind = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CERSAT-GLO-CLIM_WIND_L4-OBS_FULL_TIME_SERIE')
#print(reader_oceanwind)

# For longer simulations, it is better to pre-generate a basemap:
# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=-94, llcrnrlat=20,
                    urcrnrlon=-80, urcrnrlat=32, resolution='i')

# Add readers
#o.add_reader([reader_basemap, reader_hycom, reader_oceanwind])
o.add_reader([reader_basemap, reader_globcurrent, reader_oceanwind])

# Seed some particles
lon = -88.387161; lat = 28.736669  # Macondo location
starttime = datetime(2010, 4, 21, 6, 0, 0)  # 4 hours after explosion
time = [starttime, starttime + timedelta(hours=24*40)]
o.seed_elements(lon, lat, radius=0, number=5000, time=time)

# Run model
print(o)
o.run(duration=timedelta(days=40),
      time_step=timedelta(hours=3),
      time_step_output=timedelta(days=1))

# Print and plot results
print(o)
o.plot()
o.animation()
