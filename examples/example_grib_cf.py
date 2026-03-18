#!/usr/bin/env python
"""
Making GRIB dataset CF-compatible
=================================
"""

from datetime import datetime, timedelta
from opendrift.models.windblow import WindBlow

#%%
# Comparing / using two different wind datasets
#
# The first dataset (UCAR) is a GRIB dataset served by Thredds as a netCDF dataset, but with GRIB attributes
#   CF standard_name are added for identified surface variables (wind) with this method:
#     https://opendrift.github.io/_modules/opendrift/readers.html#add_standard_name_for_surface_grib_variables
#   This should be used with care, making sure that only surface (10m height) variables are used
#
# The second dataset (PACIOOS) is a netCDF4 dataset served by Thredds with CF attributes
wind_datasets = [
    'https://thredds.ucar.edu/thredds/dodsC/grib/NCEP/GFS/Global_0p25deg/Best',
    'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'
    ]

simulations = []

time = datetime.now() - timedelta(days=5)
for wind in wind_datasets:
    o = WindBlow(loglevel=20)
    o.add_readers_from_list(wind)
    o.seed_elements(lon=3, lat=60, time=time, number=10, radius=10000)
    o.run(duration=timedelta(hours=48))
    simulations.append(o)

#%%
# Some differences are expected as datasets may origin from different runs
simulations[0].plot(compare=simulations[1:], lscale='i', legend=['GRIB converted', 'netCDF CF compatible'])
