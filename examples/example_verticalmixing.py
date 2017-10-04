#!/usr/bin/env python

from datetime import datetime
import numpy as np

from opendrift.readers import reader_ROMS_native
from opendrift.models.pelagicegg import PelagicEggDrift

lon = 14.0; lat = 68.0;  # Vestfjorden

# Run 1, with full mixing
o1 = PelagicEggDrift(loglevel=20)
o1.set_config('turbulentmixing:max_iterations', 0)  # Default
nordic_native = reader_ROMS_native.Reader(o1.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
o1.add_reader([nordic_native])
o1.seed_elements(lon, lat, radius=1000, number=1000, z=-30,
                 time=nordic_native.start_time)
o1.run(steps=50, time_step=3600)

# Run 2, with trunkated mixing
o2 = PelagicEggDrift(loglevel=20)
o2.set_config('turbulentmixing:max_iterations', 100)
o2.add_reader([nordic_native])
o2.seed_elements(lon, lat, radius=1000, number=1000, z=-30,
                 time=nordic_native.start_time)
o2.run(steps=50, time_step=3600)


legend = ['Full vertical mixing', 'Trunkated vertical mixing']
o1.animation_profile(compare=o2, legend=legend)
o1.animation(compare=o2, legend=legend)
