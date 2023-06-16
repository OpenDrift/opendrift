#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2019, Knut-Frode Dagestad, MET Norway
# Copyright 2022, Gaute Hope, MET Norway

import os
import pytest

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

@pytest.mark.veryslow
def test_comparison(tmpdir):
    anifile = os.path.join(tmpdir, 'anim.mp4')
    plotfile = os.path.join(tmpdir, 'plot.png')
    #anifile = None
    #plotfile = None
    o = OceanDrift(loglevel=30)
    rn = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    o.add_reader(rn)
    o.seed_elements(lon=4.8, lat=60.0, number=10, radius=1000,
                    time=rn.start_time)
    o.run(steps=5)

    # Making figures/animations
    o.plot(filename=plotfile, fast=True,
           background=['x_sea_water_velocity', 'y_sea_water_velocity'])
    o.animation(color='lat', filename=anifile, show_trajectories=True, fast=True)
    o.animation(density=True, filename=anifile, show_trajectories=True, fast=True)
    o.animation(filename=anifile, show_trajectories=True, fast=True)
    o.animation(filename=anifile, fast=True,
                background=['x_sea_water_velocity', 'y_sea_water_velocity'])
    o.plot(filename=plotfile, fast=True, linecolor='lat')
    o.plot(filename=plotfile, fast=False)


    # Second run for comparison
    o2 = OceanDrift(loglevel=30)
    o2.add_reader(rn)
    o2.set_config('environment:fallback:x_wind', 15) # Adding wind
    o2.set_config('environment:fallback:y_wind', 0)
    o2.seed_elements(lon=4.8, lat=60.0, number=10, radius=1000,
                    time=rn.start_time)
    o2.run(steps=5)

    o.animation(filename=anifile, compare=o2, fast=True,
                legend=['No wind', '10 m/s wind'])
    o.plot(filename=plotfile, compare=o2, fast=True,
                legend=['No wind', '10 m/s wind'])

    # Check that files have been written
    assert os.path.exists(anifile)
    assert os.path.exists(plotfile)
