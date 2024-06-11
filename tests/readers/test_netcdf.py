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
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import os
from datetime import datetime
import unittest


from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift import OceanDrift

class TestNetCDF(unittest.TestCase):

    def test_MFDataset(self):
        reader_landmask = reader_global_landmask.Reader()

        o = OceanDrift(loglevel=30)
        nordicMF = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic_subset_day*.nc')
        nordicMF_all = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic_subset.nc')
        lon = 14.0
        lat = 67.3

        # nordicMF is missing angle
        nordicMF.Dataset["angle"] = nordicMF_all.Dataset.angle

        #nordic3d = reader_ROMS_native.Reader(o.test_data_folder() +
        #    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        #o.add_reader(nordic3d)  # Slightly different results
        # Subset is made with ncks, and should give identical result
        # e0=0, e1=20, x0=40, x1=70
        # ncks -d eta_rho,0,$e1 -d eta_psi,0,$e1 -d eta_v,0,$e1 -d eta_u,0,$e1 -d xi_rho,$x0,$x1 -d xi_psi,$x0,$x1 -d xi_v,$x0,$x1 -d xi_u,$x0,$x1 Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc Nordic_subset.nc -O --fl_fmt=netcdf4_classic
        # ncks -O -d ocean_time,0 Nordic_subset.nc Nordic_subset_day1.nc
        o.add_reader([reader_landmask, nordicMF_all])
        o.seed_elements(lon, lat, number=100, radius=5000,
                        time=nordicMF_all.start_time)
        o.run(steps=48, time_step=3600)

        print('='*99)
        # Same run, with multi-file dataset
        o2 = OceanDrift(loglevel=30)
        o2.add_reader([reader_landmask, nordicMF])
        o2.seed_elements(lon, lat, number=100, radius=5000,
                        time=nordicMF_all.start_time)
        o2.run(steps=48, time_step=3600)

        #o.plot(filename='o1.png', background='sea_floor_depth_below_sea_level')
        #o2.plot(filename='o2.png', background='sea_floor_depth_below_sea_level')

        assert o.num_elements_active() == 44
        assert o2.num_elements_active() == 44
        assert o.num_elements_deactivated() == 56
        assert o2.num_elements_deactivated() == 56
        self.assertAlmostEqual(o.elements.lon[0], o2.elements.lon[0], 5)


if __name__ == '__main__':
    unittest.main()

