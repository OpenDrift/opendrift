#!/usr/bin/env python

from datetime import timedelta
import numpy as np
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic

# Drift simulation using 10 member ensemble wind data
# from MEPS model of MET Norway

o = OceanDrift()
r = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_allmembers_extracted_2_5km_latest.nc')
o.add_reader(r)

#o.set_config('general:basemap_resolution', 'c')
o.seed_elements(lat=60, lon=4.5, time=r.start_time,
                radius=500, number=10000)

o.run(duration=timedelta(hours=50), time_step=3600)

# Ensemble members are recycled among the 10000 particles
ensemble_number = np.remainder(o.history['ID'], 10) + 1

o.animation(filename='wind_drift_ensemble.mp4',
            color=ensemble_number, clabel='Ensemble number')
