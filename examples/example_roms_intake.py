#!/usr/bin/env python
"""
ROMS 
Intake reader
==================================
"""
import numpy as np
from opendrift.readers import reader_ROMS_intake
from opendrift.models.oceandrift import OceanDrift
import intake
import xarray as xr
from datetime import timedelta

#from dask.distributed import Client
#client = Client()
#print(client.dashboard_link)

intake_catalog = 'https://mghp.osn.xsede.org/rsignellbucket1/rsignell/testing/cnaps.yml'
dataset = 'CNAPS_Forecast_Archive_64'

#intake_catalog = 'https://usgs-coawst.s3.amazonaws.com/useast-archive/coawst_intake.yml'
#dataset = 'COAWST-USEAST' 

cat = intake.open_catalog(intake_catalog)

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

cnaps = reader_ROMS_intake.Reader(dset=cat[dataset])

o.add_reader(cnaps)


# Seed elements at defined positions, depth and time
o.seed_elements(lon=-75.3, lat=35.4, radius=1000, number=10,
                z=np.linspace(0, -150, 10), time=cnaps.start_time)


# Running model
o.run(time_step=3600, duration=timedelta(hours=24))


# Print and plot results, with lines colored by particle depth
print(o)
o.plot(linecolor='z', fast=True)
#o.animation()

