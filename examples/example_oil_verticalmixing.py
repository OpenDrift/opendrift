#!/usr/bin/env python
"""
Oil vertical mixing
==================================
"""

import os
from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil


o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

ncfile = 'oilmixing.nc'
import_file = False  # Set to True to import previous run

if import_file is True:
    o.io_import_file(ncfile)
else:
    reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
    reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

    o.add_reader([reader_norkyst, reader_arome])

    # Adjusting some configuration
    #o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983')

    # Seed oil elements at defined position and time
    o.seed_elements(lon=4.9, lat=62.1, z=0, radius=1000, number=5000,
                    time=reader_arome.start_time, oil_type='GENERIC DIESEL')

    # Running model
    o.run(end_time=reader_arome.start_time + timedelta(hours=12),
          time_step=900, time_step_output=1800, outfile=ncfile)

#%%
# Print and plot results
print(o)

o.plot(linecolor='z', fast=True)
o.plot_property('z')
o.plot_oil_budget()
o.animation(fast=True)
o.animate_vertical_distribution()

os.remove(ncfile)  # cleaning up

#%%
# .. image:: /gallery/animations/example_oil_verticalmixing_0.gif
# .. image:: /gallery/animations/example_oil_verticalmixing_1.gif

