#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.set_config('general:basemap_resolution', 'i')

# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#o.add_reader([reader_norkyst, reader_arome])
o.add_reader([reader_norkyst])


# Calculating attracting/backwards FTLE/LCS at 20 hours
lcs = o.calculate_ftle(
    time=reader_norkyst.start_time + timedelta(hours=20),
    time_step=timedelta(minutes=30),
    duration=timedelta(hours=5), delta=1000,
    RLCS=False)

# Simulation from beginning and up to 30 hours (time of LCS)
o.reset()
o.seed_elements(lon=4.4, lat=60.2, number=1000, radius=1000,
                time=reader_norkyst.start_time)
o.run(end_time=reader_norkyst.start_time+timedelta(hours=20),
      time_step=timedelta(minutes=30))

o.plot(lcs=lcs, vmin=1e-7, vmax=1e-4, colorbar=True, show_particles=True)
