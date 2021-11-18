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

import numpy as np
import logging; logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
#from opendrift.elements import LagrangianArray


# Defining element properties
class LopheliaLarvae(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for Lophelia larvae
    """

    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0002}),  # For Lophelia
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 1028.}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),
        ('competence', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),
        ('hatched', {'dtype': np.float32,
                     'units': '',
                     'default': 0.})])

class LopheliaLarvaeDrift(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Lophelia larvae assumed to ber neutrally buoyant throughout

        Particles could be e.g. oil droplets, plankton, or sediments

        Under construction.
    """



    ElementType = LopheliaLarvae

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_ice_area_fraction': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
      }

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented)
    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-3000, 0]


    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(LopheliaLarvaeDrift, self).__init__(*args, **kwargs)

        # By default, elements do not strand towards coastline
        self.set_config('general:coastline_action', 'previous')

        # Vertical mixing is enabled by default
        self.set_config('drift:vertical_mixing', True)

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for larvae

        according to
        S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81  # ms-2

        # Properties that determine element buoyancy
        larvaesize = self.elements.diameter  # 0.0002 for Lophelia
        #eggsalinity = self.elements.neutral_buoyancy_salinity

        # Prepare interpolation of temperature and salinity
        # if not (Tprofiles is None and Sprofiles is None):
        #     if z_index is None:
        #         z_i = range(Tprofiles.shape[0])  # Evtl. move out of loop
        #         z_index = interp1d(-self.environment_profiles['z'],
        #                            z_i, bounds_error=False)
        #     zi = z_index(-self.elements.z)
        #     upper = np.maximum(np.floor(zi).astype(np.int), 0)
        #     lower = np.minimum(upper+1, Tprofiles.shape[0]-1)
        #     weight_upper = 1 - (zi - upper)
        #
        # # Do interpolation of temp, salt if profiles were passed into
        # # this function, if not, use reader by calling self.environment
        # if Tprofiles is None:
        #     T0 = self.environment.sea_water_temperature
        # else:
        #     T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
        #         weight_upper + \
        #         Tprofiles[lower, range(Tprofiles.shape[1])] * \
        #         (1-weight_upper)
        # if Sprofiles is None:
        #     S0 = self.environment.sea_water_salinity
        # else:
        #     S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
        #         weight_upper + \
        #         Sprofiles[lower, range(Sprofiles.shape[1])] * \
        #         (1-weight_upper)
        #
        # # The larvae have the same temperature and salinity as the ambient water (can be modified).
        # DENSw = self.sea_water_density(T=T0, S=S0)
        # DENSegg = self.sea_water_density(T=T0, S=S0)
        # dr = DENSw-DENSegg  # density difference
        #
        # W = np.zeros(len(self.elements.ID))

	    # Set vertical velocities depending on stage (age) in larval phase:

        scenario = 'short' # Set scenario, 'short' or 'long' depending on length of pre-competency period
        sday = 24*60*60 # Seconds in one day
        reftemp = 8 # degree C

        if scenario == 'short':
            dup2dwn = 21
        elif scenario == 'long':
    	    dup2dwn = 42

        idd = np.arange(15, dtype=float)
        idd = np.insert(idd, [7, 8], [6.001, 7.001])
        idd = np.append(idd, [dup2dwn, dup2dwn+1])
        idd = idd*sday*reftemp

        w_idd = np.array([0, 0, 0.00005, 0.00005, 0.00005, 0.0002, 0.00025, 0, 0, 0, 0.0003, 0.00035, 0.0004, 0.00045, 0.0005, 0.00055, 0.0006, 0.0006, -0.0006]) # Vertical velocities at specified stages

        W = np.interp(self.elements.competence, idd, w_idd)


        self.elements.terminal_velocity = W

    def update(self):
        """Update positions and properties of particles."""

        # Stokes drift
        self.stokes_drift()

        # Update element age
        self.elements.age_seconds += self.time_step.total_seconds()

        # Update competence
        self.elements.competence += self.time_step.total_seconds() * self.environment.sea_water_temperature

        # Turbulent Mixing
        self.update_terminal_velocity()
        self.vertical_mixing()

        # Horizontal advection
        self.advect_ocean_current()

        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()
