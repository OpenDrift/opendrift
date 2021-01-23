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
# Copyright 2021, Trond Kristiansen, Niva
# Added 24 Jan 2021 by Knut-Frode Dagestad, MET Norway

import datetime
import numpy as np
import pandas as pd
import pvlib
from opendrift.elements import LagrangianArray
from opendrift.models.oceandrift import Lagrangian3DArray
from opendrift.models.oceandrift import OceanDrift
import pytz
from scipy.interpolate import interp1d


class LarvalFish(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for larval and juvenile stages of fish
    """

    variables = LagrangianArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': '[]',
                                       'default': 31.25}),  # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 1028.}),
        ('survival', {'dtype': np.float32,
                      'units': 'rate',
                      'default': 1.0}),
        ('development', {'dtype': np.float32,
                         'units': 'rate',
                         'default': 0.0}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),
        ('hatched', {'dtype': np.int64,
                     'units': '',
                     'default': 0}),
        ('age', {'dtype': np.float32,  # Added to track percentage of competency time completed
                 'units': '',
                 'default': 0.}),
        ('length', {'dtype': np.float32,
                    'units': 'mm',
                    'default': 10.0}),
        ('weight', {'dtype': np.float32,
                    'units': 'mg',
                    'default': 0.08}),
        ('light', {'dtype': np.float32,
                   'units': 'ugEm2',
                   'default': 0.}),
        ('light_surface', {'dtype': np.float32,
                           'units': 'ugEm2',
                           'default': 0.}),
        ('growth_rate', {'dtype': np.float32,
                         'units': '%/dt',
                         'default': 0.}),
        ('survival', {'dtype': np.float32,
                      'units': '',
                      'default': 1.})])


class PelagicPlanktonDrift(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

    """

    ElementType = LarvalFish

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_ice_area_fraction': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
    }

    required_profiles_z_range = [0, -200]  # The depth range (in m) which profiles should cover

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}

    def __init__(self, *args, **kwargs):
        # Calling general constructor of parent class
        super(PelagicPlanktonDrift, self).__init__(*args, **kwargs)

        # IBM configuration options
        self._add_config({'IBM:attenuation_coefficient': {'type': 'float', 'default': 0.18,
                                                          'min': 0.0, 'max': 1.0, 'units': 'm-1',
                                                          'description': 'Attenuation',
                                                          'level': self.CONFIG_LEVEL_ADVANCED},
                          'IBM:fraction_of_timestep_swimming': {'type': 'float', 'default': 0.15,
                                                                'min': 0.0, 'max': 1.0, 'units': 'm-1',
                                                                'description': 'Fraction of timestep swimming',
                                                                'level': self.CONFIG_LEVEL_ADVANCED}
                          })

    def update_surface_light(self):
        # Calculate light according to the Inchein clear-sky model using pvlib for the given
        # altitude, latitude, and time of day:
        # https://pvlib-python.readthedocs.io/en/stable/clearsky.html
        dd = (self.start_time + datetime.timedelta(seconds=self.elements.age_seconds[0]))
        current_date = datetime.datetime(year=dd.year,
                                         month=dd.month,
                                         day=dd.day,
                                         hour=dd.hour,
                                         second=dd.second,
                                         tzinfo=pytz.timezone('UTC'))

        # Convert datetime objec to to pandas Timestamp for use in pvlib
        current_date = pd.DatetimeIndex([pd.Timestamp(current_date)])
        att_coeff = self.get_config('IBM:attenuation_coefficient')

        for ind in range(len(self.elements.lat)):

            if self.elements.lat[ind] < 0:
                lat = 360 - self.elements.lat[ind]
            else:
                lat = self.elements.lat[0]

            solpos = pvlib.solarposition.get_solarposition(current_date, self.elements.lat[ind], self.elements.lon[ind])
            apparent_zenith = solpos['apparent_zenith'].to_numpy()
            airmass_rel = pvlib.atmosphere.get_relative_airmass(apparent_zenith, model='kastenyoung1989')
            altitude = 0
            pressure = pvlib.atmosphere.alt2pres(altitude=altitude)
            airmass = pvlib.atmosphere.get_absolute_airmass(airmass_rel, pressure)
            linke_turbidity = pvlib.clearsky.lookup_linke_turbidity(current_date,
                                                                    self.elements.lat[ind],
                                                                    self.elements.lon[ind])
            dni_extra = pvlib.irradiance.get_extra_radiation(current_date)

            ineichen = pvlib.clearsky.ineichen(apparent_zenith, airmass, linke_turbidity, altitude, dni_extra)
            self.elements.light_surface[ind] = ineichen["ghi"].to_numpy()

        # Light at depth of larvae assuming constant attenuation (can be changed to dependent on chlorophyll)
        self.elements.light = self.elements.light_surface * np.exp(att_coeff * self.elements.z)

    def calculate_fish_growth_rate(self, wgt: np.ndarray, temp: np.ndarray) -> np.ndarray:
        # Calculate daily growth rate (SGR in percent %)
        # wgt has to be in milligram. This follows Arild Folkvord 2005
        # Folkvord, A. 2005. “Comparison of Size - at - Age of Larval Atlantic Cod(Gadus Morhua )
        # from Different Populations Based on Size - and Temperature - Dependent Growth
        # Models.” Canadian Journal of Fisheries and Aquatic Sciences.Journal Canadien Des Sciences Halieutiques
        # et Aquatiques 62(5): 1037–52.
        GR = 1.08 + 1.79 * temp - 0.074 * temp * np.log(wgt) \
             - 0.0965 * temp * np.log(wgt) ** 2 \
             + 0.0112 * temp * np.log(wgt) ** 3

        # Growth rate(g) converted to milligram weight (gr_mg) per timestep:
        sec2day = 1 / 86400.0
        g = (np.alog(GR / 100. + 1)) * sec2day * self.time_step.total_seconds()
        return wgt * (np.exp(g) - 1.)

    # Add constant mortality
    def add_daily_constant_mortality(self):
        mort = 0.17
        self.elements.survival = self.elements.survival * np.exp(-mort)

    def update_vertial_position_dynamic_range(self, length,
                                              old_light,
                                              current_light,
                                              current_depth):

        # Update the vertical position of the current larva using length (in mm) to calculate
        # ma hourly move in meters
        mm2m = 1. / 1000.
        # Swim function from Peck et al. 2006
        swim_speed = 0.5 * length * mm2m * self.time_step.total_seconds()
        max_hourly_move = swim_speed * (self.time_step.total_seconds() / 3600.)

        if old_light <= current_light:  # If light increases and stomach is sufficiently full, go down
            depth = min(0.0, current_depth - max_hourly_move)
        else:  # If light decreases or stomach is not sufficiently full, go up
            depth = min(0, current_depth + max_hourly_move)
        return depth

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for Pelagic Egg

        according to
        S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81  # ms-2

        # Pelagic Egg properties that determine buoyancy
        eggsize = self.elements.diameter  # 0.0014 for NEA Cod
        eggsalinity = self.elements.neutral_buoyancy_salinity
        # 31.25 for NEA Cod

        # prepare interpolation of temp, salt
        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper + 1, Tprofiles.shape[0] - 1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                 weight_upper + \
                 Tprofiles[lower, range(Tprofiles.shape[1])] * \
                 (1 - weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                 weight_upper + \
                 Sprofiles[lower, range(Sprofiles.shape[1])] * \
                 (1 - weight_upper)

        # The density difference between a pelagic egg and the ambient water
        # is regulated by their salinity difference through the
        # equation of state for sea water.
        # The Egg has the same temperature as the ambient water and its
        # salinity is regulated by osmosis through the egg shell.
        DENSw = self.sea_water_density(T=T0, S=S0)
        DENSegg = self.sea_water_density(T=T0, S=eggsalinity)
        dr = DENSw - DENSegg  # density difference

        # water viscosity
        my_w = 0.001 * (1.7915 - 0.0538 * T0 + 0.007 * (T0 ** (2.0)) - 0.0023 * S0)
        # ~0.0014 kg m-1 s-1

        # terminal velocity for low Reynolds numbers
        W = (1.0 / my_w) * (1.0 / 18.0) * g * eggsize ** 2 * dr

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W * 1000 * eggsize / my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * T0)  # in cm2/s
        d0 = (eggsize * 100) - 0.4 * \
             (9.0 * my_w ** 2 / (100 * g) * DENSw / dr) ** (1.0 / 3.0)  # cm
        W2 = 19.0 * d0 * (0.001 * dr) ** (2.0 / 3.0) * (my_w * 0.001 * DENSw) ** (-1.0 / 3.0)
        # cm/s
        W2 = W2 / 100.  # back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    def calculate_length(self, wgt, Larval_m):
        # Calculate length based on weight (mg)
        oldLarval_m = Larval_m
        Larval_m = np.exp(2.296 + 0.277 * np.log(wgt) - 0.005128 * np.log(wgt) ** 2)
        if Larval_m < oldLarval_m:
            Larval_m = oldLarval_m
            # print 'Kept the old length', oldLarval_m, Larval_m

        return Larval_m

    def update_fish_larvae(self):

        # LIGHT UPDATE
        # Save the light from previous time-step to use for vertical behavior
        last_light = self.elements.light_surface.copy()
        self.update_surface_light()

        dt = self.time_step.total_seconds()
        temp = self.environment.sea_water_temperature

        # FISH UPDATE
        length = self.elements.length.copy()
        self.elements.weight += self.calculate_fish_growth_rate(dt, self.elements.weight, temp)
        self.elements.length = self.calculate_length(self.elements.weight, self.elements.length)
        self.add_daily_constant_mortality()
        # Update status  according to survival probability
        self.elements.status[np.where(self.elements.survival == 0)] = 0
        self.elements.age += dt

        # VERTICAL POSITION UPDATE
        for ind in range(len(self.elements.lat)):
            self.elements.z[ind] = self.update_vertial_position_dynamic_range(self.elements.length[ind],
                                                                              last_light[ind],
                                                                              self.elements.light[ind],
                                                                              self.elements.z[ind])

    def update(self):

        """
        MAIN:

        Update positions and properties of buoyant particles.
        Here you can turn on or off which physical properties you want to include:
        vertical advection
        horizontal advection
        vertical mixingpelagicegg
        terminal velocity (sinking / buoyancy of eggs/larvae)
        vertical behavior (in update_larval_fish)
        """

        # Update element age
        self.elements.age_seconds += self.time_step.total_seconds()

        # Turn off terminal velocity (sinking) before calculating vertical mixing
        # Turn on mixing through 'processes:turbulentmixing' in run_ibm.py

        if self.stage=="egg":
            self.update_terminal_velocity()
            self.vertical_mixing()
            self.advect_ocean_current()
        else:
            self.elements.terminal_velocity = 0.0
            self.vertical_mixing()

            # Update the growth and vertical position of each individual
            self.update_fish_larvae()
            self.advect_ocean_current()
            self.vertical_advection()
