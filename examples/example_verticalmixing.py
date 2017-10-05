#!/usr/bin/env python

from datetime import datetime
import numpy as np

from opendrift.readers import reader_ROMS_native
from opendrift.models.pelagicegg import PelagicEggDrift

lon = 14.0; lat = 68.0;  # Vestfjorden

# Run 1, with full mixing
o1 = PelagicEggDrift(loglevel=0)
o1.set_config('turbulentmixing:max_iterations', 0)  # Default
nordic_native = reader_ROMS_native.Reader(o1.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
o1.add_reader([nordic_native])
o1.seed_elements(lon, lat, radius=1000, number=1000, z=-30,
                 time=nordic_native.start_time)
start = datetime.now()
o1.run(steps=48, time_step=3600)
time1 = datetime.now() - start

# Run 2, with trunkated mixing
o2 = PelagicEggDrift(loglevel=0)
o2.set_config('turbulentmixing:max_iterations', 100)
o2.add_reader([nordic_native])
o2.seed_elements(lon, lat, radius=1000, number=1000, z=-30,
                 time=nordic_native.start_time)
start = datetime.now()
o2.run(steps=48, time_step=3600)
time2 = datetime.now() - start


legend = ['Full vertical mixing (%s)' % time1,
          'Trunkated vertical mixing (%s)' % time2]
o1.animation_profile(compare=o2, legend=legend)
o1.animation(compare=o2, legend=legend)
# Compare the time spent for the two different runs,
# and the vertical mixing in particular
print '#'*50
print legend[0]
print o1.performance()
print legend[1]
print o2.performance()
