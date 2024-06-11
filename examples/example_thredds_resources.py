#!/usr/bin/env python
"""
Thredds resources for GUI
=========================
"""

from datetime import datetime
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.reader_netCDF_CF_generic import Reader

o = OceanDrift(loglevel=0)

thredds_resources = open(o.test_data_folder()+'../../opendrift/scripts/data_sources.txt').readlines()
times = {}

#%%
# Open each thredds dataset to check contents and spatial coverage
for t in thredds_resources:
    if t.startswith('http') and not t.startswith('cmems'):
        start = datetime.now()
        print('\n#%%\n%s\n' % t)
        r = Reader(t)
        print(r)
        ts = str(datetime.now() - start)
        times[t] = ts
        print('Time to open reader: ', ts)
        if r.global_coverage():
            lscale = 'coarse'
        else:
            lscale = 'intermediate'
        r.plot(lscale=lscale)

#%%
# Summary of times to open each dataset:
for t, time in times.items():
    print(time, t)
