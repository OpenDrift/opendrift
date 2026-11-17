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
#
# Copyright 2026, Mateusz Matuszak, MET Norway
# Adapted from larvalfish module for Tobis and LADIM. 

import datetime
import numpy as np
from scipy.interpolate import RectBivariateSpline
import logging; logger = logging.getLogger(__name__)
from opendrift.models.larvalfish import LarvalFishElement
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


class TobisElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for Tobis
    """
    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0009}),  # Number from https://eurofish.dk/sand-eels-the-little-fish-with-huge-economic-importance/
        ('stage_fraction', {'dtype': np.float32,  # to track percentage of development time completed
                            'units': '',
                            'default': 0.}),
        ('phase', {'dtype': np.uint8,  # 0 for eggs, 1 for larvae, 2 for fish
                     'units': '',
                     'default': 0}),
        ('hatch_rate', {'dtype': np.float32,
                    'units': '',
                    'default': 0}),
        ('length', {'dtype': np.float32,
                    'units': 'mm',
                    'default': 0})])

class Tobis(OceanDrift):
    """
    Tobis egg and larval drift model based on the OpenDrift framework. 

    Developed at MET Norway.

    Tobis eggs are demersal, meaning that they sink to the sea floor and get stuck there. 
    Larvae can swim actively in the water column.
    """

    ElementType = TobisElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.01, 'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
    }


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(Tobis, self).__init__(*args, **kwargs)

        # IBM configuration options
        self._add_config({
            'IBM:fraction_of_timestep_swimming':
                {'type': 'float', 'default': 0.15,
                 'min': 0.0, 'max': 1.0, 'units': 'fraction',
                 'description': 'Fraction of timestep swimming',
                 'level': CONFIG_LEVEL_ADVANCED},
            })

        # The eggs are considered demersal, meaning that they are stuck to the sea floor. 
        # By default they will therefore stay in place until hatched. 
        # However, egg advection with ocean currents may be turned on if wished. 
        self._add_config({
            'drift:egg_advection':
            {'type': 'bool', 'default':False,
             'description':'Turn on/off egg advection',
             'level':CONFIG_LEVEL_ADVANCED}
        })

        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:vertical_mixing_at_surface', True)
        self._set_config_default('drift:vertical_advection_at_surface', True)

        self.hatch_time = self.get_hatch_time_func()

    def get_hatch_time_func(self):
        """
        Total hatch time based on Smigielski et al. (1984), doi: 10.3354/meps014287

        The underlying data table is taken directly from the paper. Interpolation is
        linear with date and second order in the rate direction.
        """

        days_tab = np.array([
            [61, 51, 39, 25],
            [82, 67, 48, 30],
            [135, 116, 82, 55],
        ])
        temp_tab = np.array([2, 4, 7, 10])
        rate_tab = np.array([0, 0.5, 1])

        spline = RectBivariateSpline(
            x=rate_tab,
            y=temp_tab,
            z=days_tab,
            kx=2,
            ky=1,
        )

        def hatch_time_fn(rate, temp):
            """
            Total hatch time based on Smigielski et al. (1984), doi: 10.3354/meps014287

            The underlying data table is taken directly from the paper. Interpolation is
            linear with date and second order in the rate direction.

            :param rate: 0 = earliest spawners, 0.5 = median spawners, 1 = latest spawners
            :param temp: Ambient temperature, in degrees Celcius
            :return: Total hatch time, in days
            """
            temp = np.minimum(temp_tab[-1], np.maximum(temp_tab[0], temp))
            out = spline(rate.ravel(), temp.ravel(), grid=False)
            return out.reshape(rate.shape)

        return hatch_time_fn
    
    def update_egg_stage_fraction(self):
        """Egg development according to Christensen et al. (2008), doi:10.1139/F08-073

        Using initialization model `e` (variable maturation)"""

    
        eggs = np.where(self.elements.phase==0)[0]
        if len(eggs) > 0:
            dev_days = self.hatch_time(self.elements.hatch_rate[eggs], self.environment.sea_water_temperature[eggs])
            stage_increase = self.time_step.total_seconds()/(dev_days*24*60*60)
            self.elements.stage_fraction[eggs] += stage_increase
            hatching = np.where(self.elements.stage_fraction[eggs]>=1)[0]
            if len(hatching) > 0:
                logger.debug('Hatching %s eggs' % len(hatching))
                self.elements.phase[eggs[hatching]] = 1  # Eggs with stage_fraction >= 1 are hatched (1)
        
        
    def update_larval_stage_fraction(self):
        """Larval development according to Christensen et al. (2008), doi:10.1139/F08-073

        The function modifies (in-place) the variable `stage`. The length of the larvae is
        assumed to be related to the stage as

        L = L0 * (1 - s) + Lm * s

        where L0 = 7.73mm, Lm = 40mm and s = stage - 1"""
        
        larvae = np.where(self.elements.phase==1)[0]

        if len(larvae) > 0:
            #double check simplify ?
            s = self.elements.stage_fraction[larvae] - 1
            Lm = 40
            L0 = 7.73
            L_inf = 218
            
            L = L0 + s * (Lm - L0)
            gam = 0.316
            lamb0 = -1.725
            lamb1 = 0.136
            lamb = np.exp(lamb0 + lamb1 * self.environment.sea_water_temperature[larvae])

            dLdt = lamb * np.power(L / L0, gam) * (1 - L / L_inf)
            L_new = L + dLdt * self.time_step.total_seconds()/86400
            self.elements.stage_fraction[larvae] = 1 + (L_new - L0) / (Lm - L0)
            self.elements.length[larvae] = L_new

            metamorphing = np.where(self.elements.stage_fraction[larvae] >= 2)[0]
            if len(metamorphing) > 0:
                logger.info('Metamorphing %s larvae' % len(metamorphing))
                self.elements.phase[larvae[metamorphing]] = 2  # Larvae with stage_fraction >= 2 are metamorphing (2)
                self.deactivate_metamorphed()

    def larvae_vertical_migration(self):

        larvae = np.where(self.elements.phase==1)[0]
        if len(larvae) == 0:
            return

        # Vertical migration of Larvae
        # Swim function from Peck et al. 2006
        L = self.elements.length[larvae]
        swim_speed = (0.261*(L**(1.552*L**(-0.08))) - 5.289/L) / 1000
        f = self.get_config('IBM:fraction_of_timestep_swimming')
        max_migration_per_timestep = f*swim_speed*self.time_step.total_seconds()

        # Using here UTC hours. Should be changed to local solar time,
        # although a phase shift of some hours should not make much difference
        if self.time.hour < 12:
            direction = -1  # Swimming down when light is increasing
        else:
            direction = 1  # Swimming up when light is decreasing

        self.elements.z[larvae] = np.minimum(0, self.elements.z[larvae] + direction*max_migration_per_timestep)
    
    def freeze_egg(self):
        eggs = np.where(self.elements.phase==0)[0]
        larvae = np.where(self.elements.phase==1)[0]
        if len(eggs) > 0:
            self.elements.moving[eggs] = 0
        if len(larvae) > 0:
            self.elements.moving[larvae] = 1
        
    def deactivate_metamorphed(self):
        metamorphed_indices = [el >= 2 for el in self.elements.phase]            
        self.deactivate_elements(metamorphed_indices, reason='Metamorphosed')

    def update(self):

        self.update_egg_stage_fraction()
        self.update_larval_stage_fraction()
        if self.get_config('drift:egg_advection') is False:
            self.freeze_egg()
        self.advect_ocean_current()
        self.stokes_drift()

        self.vertical_mixing()
        self.larvae_vertical_migration()