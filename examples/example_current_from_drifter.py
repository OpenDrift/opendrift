#!/usr/bin/env python
"""
Current from drifter
====================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_current_from_drifter
from opendrift.models.oceandrift import OceanDrift


o = OceanDrift(loglevel=20)
o.set_config('environment:fallback:land_binary_mask', 0)

#%%
# We make a reader which reconstructs the ocean current from
# observed time series of a drifter
# This is actual data of SLDMB/Code drifter as used in this study:
# Jones, C.E., Dagestad, K.-F., Breivik, O., Holt, B., Rohrs, J., Christensen, K.H., Espeseth, M.M., Brekke, C., Skrunes, S. (2016): Measurement and modeling of oil slick transport. Journal of Geophysical Research - Oceans, Volume 121, Issue 10, October 2016, Pages 7759-7775. DOI: 10.1002/2016JC012113.

drifterlons = [2.407376, 2.405140, 2.403248, 2.401872, 2.400152, 2.398518, 2.397056, 2.395766, 2.394476, 2.393358, 2.392584, 2.391810, 2.390606, 2.389316, 2.388628, 2.388370, 2.387940, 2.387510, 2.387338, 2.387166, 2.387252, 2.387338, 2.387682, 2.387854, 2.388284, 2.388628, 2.389230, 2.390004, 2.390434, 2.390692, 2.391380, 2.391896, 2.392068, 2.392154, 2.392068, 2.391896, 2.391896, 2.391896, 2.391638, 2.391380, 2.391208, 2.391036, 2.390692, 2.390090, 2.389660, 2.389058, 2.388628]
drifterlats = [60.034740, 60.033880, 60.033106, 60.032246, 60.031300, 60.030182, 60.028892, 60.027602, 60.026656, 60.025538, 60.024420, 60.023388, 60.022442, 60.021496, 60.020378, 60.019346, 60.018572, 60.017626, 60.016852, 60.016164, 60.015734, 60.015304, 60.014616, 60.014100, 60.013670, 60.013412, 60.013240, 60.013068, 60.013154, 60.013412, 60.013584, 60.013842, 60.014186, 60.014616, 60.015218, 60.015820, 60.016594, 60.017454, 60.018400, 60.019346, 60.020464, 60.021410, 60.022442, 60.023474, 60.024678, 60.025882, 60.026914]
drifterlats = drifterlats[::-1]
drifterlons = drifterlons[::-1]
driftertimes = [datetime(2015, 6, 10, 5, 50) +
    timedelta(minutes=10)*i for i in range(len(drifterlons))]

r = reader_current_from_drifter.Reader(
        lons=drifterlons, lats=drifterlats, times=driftertimes)
o.add_reader(r)

#%%
# We seed elements within polygon, as could have been extracted
# from remote sensing imagery
lons = [2.39, 2.391, 2.392, 2.393, 2.394, 2.393, 2.392, 2.391, 2.39]
lats = [60.02, 60.02, 60.019, 60.02, 60.021, 60.022, 60.021, 60.021, 60.02]
o.seed_within_polygon(lons=lons, lats=lats,
                      number=2000, time=r.start_time)

#%%
# Finally running simulation
o.run(end_time=r.end_time, time_step=r.time_step)

o.animation(buffer=.01, fast=True, drifter={'time': driftertimes, 'lon': drifterlons, 'lat': drifterlats,
    'label': 'CODE Drifter', 'color': 'b', 'linewidth': 2, 'markersize': 40})

#%%
# .. image:: /gallery/animations/example_current_from_drifter_0.gif

#%%
# Drifter track is shown in red, and simulated trajectories are shown in gray. Oil spill is displaced relative to drifter, but drifter current is assumed to be spatially homogeneous.
o.plot(buffer=.01, fast=True, drifter={
        'lon': drifterlons, 'lat': drifterlats,
        'time': driftertimes, 'linewidth': 2, 'color': 'r'})

