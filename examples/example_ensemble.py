#!/usr/bin/env python
"""
Ensemble
==================================
"""

from datetime import timedelta
import numpy as np
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic

#%%
# Drift simulation using 30 member ensemble wind data
# from MEPS model of MET Norway

o = OceanDrift(loglevel=20)
o.set_config('drift:vertical_mixing', False)
r = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
o.add_reader(r)

o.seed_elements(lat=60, lon=4.9, time=r.start_time,
                radius=1000, number=10000)

o.run(duration=timedelta(hours=50), time_step=600, time_step_output=3600)

#%%
# Ensemble members are recycled among the 10000 particles
ensemble_number = np.remainder(np.arange(o.num_elements_total()), len(r.realizations)) + 1
o.animation(fast=True, color=ensemble_number, legend=['Member ' + str(i) for i in r.realizations], colorbar=False)

#%%
# We see that elements forced by different wind ensemble members move differently

#%%
# .. image:: /gallery/animations/example_ensemble_0.gif

#%% Finally, plot memory usage during simulation
o.plot_memory_usage()
