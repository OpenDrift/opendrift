#!/usr/bin/env python
"""
Mixing down to Mixed Layer Depth
================================
"""

import os
import tempfile
from datetime import datetime, timedelta
import numpy as np
import xarray as xr
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.reader_constant import Reader as ConstantReader
from opendrift.readers.reader_netCDF_CF_generic import Reader

#%%
# Mixed Layer Depth of 20m West of 3 deg E, and 50m to the east
r1 = ConstantReader({'ocean_mixed_layer_thickness': 20})
r2 = ConstantReader({'ocean_mixed_layer_thickness': 50})
r1.xmax = 3
r2.xmin = 3

#%%
# First with Sundby1983 parameterization of diffusivity, based on wind and MLD
o = OceanDrift(loglevel=50)
o.add_reader([r1, r2])
o.set_config('environment:constant:y_wind', 8)  # Some wind for mixing
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983')
# Increasing background diffusivity beyond default (1.2e-5) to avoid artefact due to sharp gradient at MLD
o.set_config('vertical_mixing:background_diffusivity', 0.001)
o.seed_cone(lon=[2, 4], lat=[60, 60], time=datetime.now(), number=5000, z=-1)
o.run(duration=timedelta(hours=48))
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_0.gif

#%%
# Same, but with Large1994 parameterization of diffusivity
o = OceanDrift(loglevel=50)
o.add_reader([r1, r2])
o.set_config('environment:constant:y_wind', 8)  # Some wind for mixing
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
o.set_config('vertical_mixing:background_diffusivity', 0.001)
o.seed_cone(lon=[2, 4], lat=[60, 60], time=datetime.now(), number=5000, z=-1)
o.run(duration=timedelta(hours=48))
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_1.gif

#%%
# Using Large1994, but with 0 diffusivity below MLD
o = OceanDrift(loglevel=50)
o.add_reader([r1, r2])
o.set_config('environment:constant:y_wind', 8)  # Some wind for mixing
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
o.set_config('vertical_mixing:background_diffusivity', 0)
o.seed_cone(lon=[2, 4], lat=[60, 60], time=datetime.now(), number=5000, z=-1)
o.run(duration=timedelta(hours=48))
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_2.gif

#%%
# Diffusivity from an ocean model, given on a few depth levels: two NorKyst800 profiles
# (2026-11-20 00 UTC) with a mixed layer of about 25 m (4.5E 62N) and about 55 m (3.5E 60N),
# each written to a small file and read like ocean model output, and seeded at its own longitude
levels = np.array([0, 1, 2, 3, 5, 7, 10, 15, 25, 50, 65, 75, 100])
K25 = [1e-6, 7e-3, 2.5e-2, 3.3e-2, 4.6e-2, 5.6e-2, 6.3e-2, 4.3e-2, 4.5e-5, 1.1e-6, 1.1e-6, 1.1e-6, 1.1e-6]
K55 = [1e-6, 7.8e-3, 1.9e-2, 2.4e-2, 3.5e-2, 4.5e-2, 6.1e-2, 8.8e-2, 1.2e-1, 6.7e-3, 1.9e-6, 1.1e-6, 1.1e-6]
# The values above are the AKs of NorKyst800 on thredds, rounded; uncomment to read them
# (or profiles from another time or position) from thredds instead:
# r = Reader('https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be')
# env, prof = r.get_variables_interpolated(
#     ['ocean_vertical_diffusivity'], profiles=['ocean_vertical_diffusivity'], profiles_depth=100,
#     time=datetime(2026, 11, 20), lon=np.array([4.5, 3.5]), lat=np.array([62., 60.]), z=np.zeros(2))
# levels = -np.asarray(prof['z'])
# K25, K55 = np.asarray(prof['ocean_vertical_diffusivity']).T
t0 = datetime.now().replace(minute=0, second=0, microsecond=0)


def run_with_profile(K, lons):
    lon, lat = np.arange(0, 6.01, .5), np.arange(59, 61.01, .5)
    ds = xr.Dataset(
        {'AKs': (('time', 'depth', 'lat', 'lon'), np.broadcast_to(np.array(K)[None, :, None, None],
                                                                 (2, len(levels), len(lat), len(lon))).copy(),
                 {'standard_name': 'ocean_vertical_diffusivity', 'units': 'm2 s-1'})},
        coords={'time': [t0, t0 + timedelta(days=3)],
                'depth': ('depth', levels, {'standard_name': 'depth', 'positive': 'down', 'units': 'm'}),
                'lat': ('lat', lat, {'standard_name': 'latitude', 'units': 'degrees_north'}),
                'lon': ('lon', lon, {'standard_name': 'longitude', 'units': 'degrees_east'})})
    diffusivity_file = os.path.join(tempfile.mkdtemp(), 'diffusivity_levels.nc')
    ds.to_netcdf(diffusivity_file)
    o = OceanDrift(loglevel=50)
    o.add_reader(Reader(diffusivity_file))
    o.set_config('drift:vertical_mixing', True)
    o.set_config('vertical_mixing:diffusivitymodel', 'environment')
    o.seed_cone(lon=lons, lat=[60, 60], time=t0, number=2500, z=-1)
    o.run(duration=timedelta(hours=48))
    return o


o55 = run_with_profile(K55, [3.25, 3.75])
o25 = run_with_profile(K25, [4.25, 4.75])
o55.animation_profile(compare=o25, legend=['NorKyst800 3.5E 60N, mixed layer ~55 m',
                                           'NorKyst800 4.5E 62N, mixed layer ~25 m'])

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_3.gif
