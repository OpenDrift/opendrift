#!/usr/bin/env python
"""
Convolve input
==============

Decreasing the spatial resolution of fields from a reader by convolution.
This may improve accuracy, see: https://doi.org/10.1016/j.rse.2019.01.001
"""

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift


lon = 4.9; lat = 60.0
o = OceanDrift(loglevel=20)

reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
time = reader_norkyst.start_time

o.add_reader([reader_norkyst])
o.seed_elements(lon, lat, radius=1000, number=1000, time=time)
o.run(steps=20)

#%%
# Store final field of x-component of current
original_current = reader_norkyst.var_block_after[list(reader_norkyst.var_block_after.keys())[0]].data_dict['x_sea_water_velocity'].copy()

#%%
# For the second run, the NorKyst currents are convolved with a kernel,
# effectively lowering the spatial resolution.
# <reader>.set_convolution_kernel may also be given as an array (kernel) directly
N = 10  # Convolusion kernel size
reader_norkyst.set_convolution_kernel(N)  # Using convolution kernel for second run
o2 = OceanDrift(loglevel=20)
o2.add_reader([reader_norkyst])
o2.seed_elements(lon, lat, radius=1000, number=1000, time=time)
o2.run(steps=20)

#%%
# Store final field of x-component of (convolved) current
convolved_current = reader_norkyst.var_block_after[
    "['x_sea_water_velocity', 'y_sea_water_velocity']"].data_dict['x_sea_water_velocity']

plt.subplot(2,1,1)
plt.imshow(original_current, interpolation='nearest')
plt.title('Original current field (x-component)')
clim = plt.gci().get_clim()
plt.colorbar()
plt.subplot(2,1,2)
plt.imshow(convolved_current, interpolation='nearest')
plt.clim(clim)  # Make sure plots are comparable
plt.colorbar()
plt.title('Final, convolved current field (x-component)')
plt.show()

#%%
# Print and plot results, with convolved currents as background
print(o)
o.animation(compare=o2, fast=True, legend=[
    'Original currents', 'Current convoled with kernel of size %s' % N],
    background=['x_sea_water_velocity', 'y_sea_water_velocity'])

#%%
# .. image:: /gallery/animations/example_convolve_input_0.gif

o.plot(fast=True)
