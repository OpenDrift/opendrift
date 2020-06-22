#!/usr/bin/env python
"""
Openoil sample output netCDF file
==================================
"""

import os
from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_constant
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20, weathering_model='noaa')

rc = reader_constant.Reader({'x_wind': 5,
                             'y_wind': 4})

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([rc, reader_norkyst, reader_arome])

#%%
# Seed oil elements at defined position and time
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=1)]
o.seed_elements(lon=5.05, lat=59.95, radius=1500, number=100,
                time=time, z=0, m3_per_hour=2, oiltype='TAU 1999')

#%%
# Running model
o.run(steps=4*8, time_step=900, time_step_output=3600, outfile= 'openoil_sample_output.nc')

#%%
# Print and plot results
o.plot_oil_budget()
o.plot()

os.system('ncdump -h openoil_sample_output.nc')
