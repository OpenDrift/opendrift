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

import unittest
from datetime import datetime, timedelta
import numpy as np
import pyproj

from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_constant
from opendrift.models.plastdrift import PlastDrift
from opendrift.models.openoil import OpenOil
from opendrift.models.windblow import WindBlow
from opendrift.models.shipdrift import ShipDrift
from opendrift.models.openberg import OpenBerg
from opendrift.models.larvalfish import LarvalFish

import opendrift
print(opendrift.versions())

class TestModels(unittest.TestCase):
    """Tests for OpenDrift models"""

    def test_seed(self):
        """Test seeding"""
        self.o = OpenOil()
        self.fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
        self.fake_eddy.start_time = datetime(2015, 1, 1)
        self.o.add_reader([self.fake_eddy])
        self.o.seed_elements(lon=4, lat=60, number=100,
                             time=self.fake_eddy.start_time)
        self.assertEqual(len(self.o.elements), 0)
        self.assertEqual(len(self.o.elements_deactivated), 0)
        self.assertEqual(len(self.o.elements_scheduled), 100)

    def test_windblow(self):
        o = WindBlow(loglevel=30)
        reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF_20160202_subset.nc')
        o.add_reader([reader_arome])
        lat = 67.711251; lon = 13.556971  # Lofoten
        o.seed_elements(lon, lat, radius=5000, number=1000,
                        time=reader_arome.start_time)
        o.run(steps=24, time_step=3600)
        self.assertAlmostEqual(o.elements.lon.max(), 15.864, 2)

    def test_shipdrift(self):
        """Sintef case study"""
        s = ShipDrift(loglevel=50)
        s.set_config('drift:current_uncertainty', 0)
        s.set_config('drift:wind_uncertainty', 0)
        c = reader_constant.Reader({
            'sea_surface_wave_significant_height': 5,
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 11,
            'x_wind': 14.14213562,
            'y_wind': -14.14213562,
            'x_sea_water_velocity': 0.05656854249,
            'y_sea_water_velocity': -0.05656854249})
        s.set_config('environment:fallback:land_binary_mask', 0)
        s.add_reader(c)
        s.list_configspec('seed')
        s.set_config('seed:orientation', 'left')
        s.seed_elements(lon=2, lat=60, time=datetime.now(), number=3,
                        length=80, beam=14, height=25, draft=5)
        s.run(time_step=600, duration=timedelta(hours=4))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            s.elements.lon, 2.252, 3))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            s.elements.lat, 59.876, 3))

    def test_shipdrift_defaults(self):
        s = ShipDrift(loglevel=0)
        #s.list_configspec()
        s.set_config('environment:fallback:land_binary_mask', 0)
        c = reader_constant.Reader({
            'sea_surface_wave_significant_height': 5,
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 11,
            'x_wind': 14.14213562,
            'y_wind': -14.14213562,
            'x_sea_water_velocity': 0.05656854249,
            'y_sea_water_velocity': -0.05656854249})
        s.add_reader(c)
        s.set_config('seed:height', 14)
        s.seed_elements(lon=2, lat=60, time=datetime.now(), number=1)
        s.run(duration=timedelta(hours=4))
        #self.assertAlmostEqual(s.elements.lon.max(), 2.1273, 3)  # Without setting config
        self.assertAlmostEqual(s.elements.lon.max(), 2.1990, 3)

    def test_shipdrift_backwards(self):
        """Case above, reversed"""
        s = ShipDrift(loglevel=50)
        c = reader_constant.Reader({
            'sea_surface_wave_significant_height': 5,
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 11,
            'x_wind': 14.14213562,
            'y_wind': -14.14213562,
            'x_sea_water_velocity': 0.05656854249,
            'y_sea_water_velocity': -0.05656854249})
        s.set_config('environment:fallback:land_binary_mask', 0)
        s.add_reader(c)
        s.seed_elements(lon=2.25267706, lat=59.87694775,
                        time=datetime.now(), number=1,
                        length=80, beam=14, height=25, draft=5)
        s.run(time_step=-600, duration=timedelta(hours=4))
        self.assertIsNone(np.testing.assert_array_almost_equal(
                s.elements.lon, 2.0, 3))
        self.assertIsNone(np.testing.assert_array_almost_equal(
                s.elements.lat, 60, 3))

    def test_wind_drift_shear(self):
        """Testing PlastDrift model, with wind-induced current shear"""
        o = PlastDrift(loglevel=30)
        o.set_config('environment:fallback:x_wind',  10)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lat=60, lon=5, time=datetime.now(),
                        number=3,
                        z = np.array([0, -0.05, -.1]))
        o.run(duration=timedelta(hours=10))
        self.assertIsNone(np.testing.assert_array_almost_equal(
                            o.elements.lon,
                          [5.011935,5.01738,5.011235]))
        self.assertAlmostEqual(o.elements.lat[0], o.elements.lat[2], 3)

    def test_openberg(self):
        """Check if weighting array is set correctly
        and if model returns expected positions"""
        o = OpenBerg()
        o.set_config('drift:current_uncertainty', 0)
        o.set_config('drift:wind_uncertainty', 0)

        reader_current = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
                '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

        reader_landmask = reader_global_landmask.Reader()

        o.add_reader([reader_current,reader_landmask])
        o.seed_elements(4.,62.,time=reader_current.start_time)
        o.run(steps=1)

        arr=[0.16072658,0.16466097,0.17384121,0.17325179,0.1715925,0.15592695]

        for indx in range(len(arr)):
            self.assertAlmostEqual(o.uw_weighting[indx],arr[indx],8)

        self.assertAlmostEqual(o.history['lon'].data[0][1],3.991, 3)
        self.assertAlmostEqual(o.history['lat'].data[0][1],62.011, 3)

    def test_oil_in_ice(self):
        """ Testing ice-in-oil transport with
        different values of sea ice concentration as defined by Nordam et al. 2019"""

        c = [0.2, 0.5, 0.8]
        lon = 24; lat = 81

        # Distances calculated with fallback_values and Nordam's equation
        distances = {'0.2': 21.2914, '0.5': 15.1405, '0.8': 7.2}

        geod = pyproj.Geod(ellps='WGS84')

        for i in c:
            o = OpenOil(loglevel=50)
            o.set_config('environment:fallback:x_wind', 0)  # zonal wind
            o.set_config('environment:fallback:y_wind', 4)  # meridional wind

            o.set_config('environment:fallback:x_sea_water_velocity', 0)  # eastward current
            o.set_config('environment:fallback:y_sea_water_velocity', .4)  # meridional current

            o.set_config('environment:fallback:sea_ice_x_velocity', 0)  # zonal ice velocity
            o.set_config('environment:fallback:sea_ice_y_velocity', .2)  # meridional ice velocity

            o.set_config('environment:fallback:sea_surface_wave_stokes_drift_y_velocity', 0.1) # meridional Stokes drif

            o.set_config('processes:dispersion',  False)
            o.set_config('processes:evaporation',  False)
            o.set_config('processes:emulsification',  False)
            o.set_config('drift:stokes_drift',  True)
            o.set_config('processes:update_oilfilm_thickness',  False)
            o.set_config('drift:current_uncertainty',  0)
            o.set_config('drift:wind_uncertainty',  0)

            o.set_config('environment:fallback:sea_ice_area_fraction', i)

            o.seed_elements(lon, lat, radius=1, number=10, time=datetime.now(), wind_drift_factor=0.03)

            o.run(duration=timedelta(hours=10))

            latf = o.history['lat'][0][-1]
            lonf = o.history['lon'][0][-1]

            azimuth1, azimuth2, dist = geod.inv(lon, lat, lonf, latf)

            self.assertAlmostEqual(distances[str(i)], dist/1000, 2)


    def test_larvalfish(self):
        o = LarvalFish()
        # Tests to be added

if __name__ == '__main__':
    unittest.main()
