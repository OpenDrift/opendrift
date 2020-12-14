#!/usr/bin/env python
"""
Analysing huge output files
===========================
"""

from datetime import datetime, timedelta
import opendrift
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_oscillating

#%
# First make a simulation with two seedings, marked by *origin_marker*
o = OceanDrift(loglevel=50)
t1 = datetime.now()
t2 = t1 + timedelta(hours=6)
number = 10000
outfile = 'simulation.nc'  # Raw simulation output
o.seed_elements(time=t1, lon=4, lat=60, number=number,
                origin_marker=0)
o.seed_elements(time=[t1, t2], lon=4.2, lat=60.4, number=number,
                origin_marker=1)

reader_x = reader_oscillating.Reader('x_sea_water_velocity',
                amplitude=1, zero_time=t1)
reader_y = reader_oscillating.Reader('y_sea_water_velocity',
                amplitude=1, zero_time=t2)
o.add_reader([reader_x, reader_y])
o.set_config('drift:horizontal_diffusivity', 10)
o.run(duration=timedelta(hours=24),
      time_step=900, time_step_output=1800, outfile=outfile)

#%
# Opening the output file lazily with Xarray.
# This will work even if the file is too large to fit in memory, as it
# will read and process data chuck-by-chunk directly from file using Dask.
# Note that the analysis file will be re-used if existing. Thus this file should be deleted after making any changes to the simulation above.
o = opendrift.open_xarray(outfile, analysis_file='simulation_density.nc')

#%
# Making two animations, for each of the two seedings / origin_markere.
# The calculated density fields will be stored/cached in the analysis file
# for later re-use, as their calculation may be time consuming
# for huge output files.
# Note that other analysis/plotting methods are not yet adapted
# to datasets opened lazily with open_xarray
for om in [0, 1]:
    o.animation(density=True, density_pixelsize_m=500, fast=False,
                origin_marker=om, show_elements=False, vmin=0, vmax=200)

#%
# First seeding

#%%
# .. image:: /gallery/animations/example_huge_output_0.gif

#%
# Second seeding

#%%
# .. image:: /gallery/animations/example_huge_output_1.gif
