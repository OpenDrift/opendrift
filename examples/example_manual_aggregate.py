#!/usr/bin/env python
"""
Manual aggregate
==================
"""

from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import open_mfdataset_overlap
from opendrift.models.oceandrift import OceanDrift

#%
# Create manual aggregate from individual URLs, for NorKyst ocean model initialized at 00 hours every day

start_time = datetime.now().date()-timedelta(days=3)
end_time = datetime.now().date()-timedelta(days=2)
ds = open_mfdataset_overlap(
    'https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_his_files/%Y/%m/%d/norkyst800_his_zdepth_%Y%m%dT00Z_m00_AN.nc',
    time_series=pd.date_range(start_time, end_time, freq='1D'))

#%
# Create reader from Xarray dataset
rm = reader_netCDF_CF_generic.Reader(ds, name='NorKyst manual aggregate')
print(rm)

om = OceanDrift()
om.add_reader(rm)
om.seed_elements(lon=4.5, lat=60.0, number=1000, radius=100, time=rm.start_time)
om.run(end_time=rm.end_time)

#%
# Second simulation using ready made aggregate from thredds
ot = OceanDrift()
ot.add_readers_from_list(['https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'])
ot.seed_elements(lon=4.5, lat=60.0, number=1000, radius=100, time=rm.start_time)
ot.run(end_time=rm.end_time)

#%
# Simulation should be identical, but we see that manual aggregate is significantly slower than using thredds aggregate
om.animation(compare=ot,
             legend=[f'NorKyst manual aggregate {om.timing["total time"]}',
                     f'NorKyst thredds aggregate {ot.timing["total time"]}'])

#%%
# .. image:: /gallery/animations/example_manual_aggregate_0.gif
