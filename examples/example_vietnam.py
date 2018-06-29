#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Adding readers for global Thredds datasets:
# - Ocean forecast from UK Met Office (FOAM model)
# - Weather forecast from NOAA/NCEP
o.add_readers_from_list([
    'http://data.ncof.co.uk/thredds/dodsC/METOFFICE-GLO-AF-PHYS-HOURLY-CUR',
    'http://oos.soest.hawaii.edu/thredds/dodsC/hioos/model/atm/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])



# Seed some particles
lat=10.0; lon=107.8
#time = datetime(2010, 3, 23, 6, 0, 0)
time = datetime.now()
o.seed_elements(lon, lat, radius=1000, number=1000, time=time)

# Run model
print(o)
o.run(duration=timedelta(days=5),
      time_step=timedelta(minutes=15),
      time_step_output=timedelta(hours=3))

# Print and plot results
print(o)
o.plot()
o.animation()
