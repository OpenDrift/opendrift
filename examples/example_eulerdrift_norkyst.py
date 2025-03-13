"""
Euler simulation / Finite difference of blob with the Norkyst nordic ocean model
================================================================================
"""
import logging
import opendrift
import matplotlib.pyplot as plt
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.models import eulerdrift


reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
# reader_norkyst = reader_netCDF_CF_generic.Reader(OpenDriftSimulation.test_data_folder(None) +  '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
lon0, lat0 = reader_norkyst.xy2lonlat(reader_norkyst.xmin, reader_norkyst.ymin)

s = eulerdrift.ExplSimulation.new(5.21, 60.19, 10., shape = (400, 400))
s.t0 = reader_norkyst.start_time
s.readers.append(eulerdrift.OpendriftReader(reader_norkyst))

print(s.U(0.0))

loc, lac = s.grid.center()
s.source_gaussian_blob(loc, lac, 1., 100, 50.)

#%%
# Plot the initial conditions
s.grid.plot()
plt.show()

#%%
# Integrate
s.integrate(dt = 20., max_steps=300)

#%%
# Plot the result
s.grid.plot()
plt.show()