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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import unittest
import os
from datetime import datetime

from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift3D import OceanDrift3D


class Test(unittest.TestCase):
    """Tests for netCDF"""

    def test_lazy(self):
        o = OceanDrift3D()
        w = 'http://thredds.met.no/thredds/dodsC/meps25files/meps_det_pp_2_5km_20191128T06Z.nc'
        o.add_readers_from_list([w], lazy=True)  # Works if lazy=False
        o.fallback_values['x_sea_water_velocity'] = .3
        o.fallback_values['y_sea_water_velocity'] = .3
        o.seed_elements(lon=4.6, lat=60, radius=1000, number=1000,
                        time=datetime(2019, 11, 28, 12))
        print(o)
        o.run(steps=1, outfile='test.nc')  # Fails for Python3
        assert os.path.exists('test.nc')

    def test_MFDataset(self):

        reader_landmask = reader_global_landmask.Reader(
            llcrnrlon=13.5, llcrnrlat=67.1,
            urcrnrlon=14.6, urcrnrlat=67.7)

        o = OceanDrift3D(loglevel=0)
        nordicMF = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic_subset_day*.nc')
        nordicMF_all = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic_subset.nc')
        lon = 14.0
        lat = 67.3
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

        # Same run, with multi-file dataset
        o2 = OceanDrift3D(loglevel=30)
        o2.add_reader([reader_landmask, nordicMF])
        o2.seed_elements(lon, lat, number=100, radius=5000,
                        time=nordicMF_all.start_time)
        o2.run(steps=48, time_step=3600)

        #o.plot(filename='o1.png', background='sea_floor_depth_below_sea_level')
        #o2.plot(filename='o2.png', background='sea_floor_depth_below_sea_level')

        self.assertEqual(o.num_elements_active(), 33)
        self.assertEqual(o2.num_elements_active(), 33)
        self.assertEqual(o.num_elements_deactivated(), 67)
        self.assertEqual(o2.num_elements_deactivated(), 67)
        self.assertEqual(o.elements.lon[0], o2.elements.lon[0])

if __name__ == '__main__':
    unittest.main()
