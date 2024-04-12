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
# Copyright 2015, 2020, Knut-Frode Dagestad, MET Norway

"""
This file is a template for an ocean trajectory model
to be implemented within the OpenDrift framework.
To create a new model, make a copy (renamed) of this file, and
modify the class names and details below, according to the comments.
Alternatively, another existing class can be used as a starting point.

For more help, see https://opendrift.github.io/writing_a_new_model.html
For usage of the new model, see example scripts.

First we define an element type to be used by the trajectory model.
An element type (class) may also be imported/reused from
another module - see examples in the folder "elements"

"""

import numpy as np
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift import Lagrangian3DArray

class TemplateElementType(Lagrangian3DArray):
    """Extending LagrangianArray with relevant properties."""

    # Define and name the properties which the elements shall have.
    # These are added to the inherited properties (e.g. from OceanDrift)
    # Property names may be freely chosen, but the "update" function (below)
    # will refer to the names in the specifications of the processes.
    #
    # The "default" value is used for alle elements, unless a custom
    # value is provided at time of seeding:
    # o.seed_elements(......., <property>=<value>)
    #
    # The default value can also be updated by the user
    # before seeding with the config mechanism:
    # o.set_config('seed:<property>', <new default value>)
    #
    # The "units" attribute is used for labeling.

    variables = Lagrangian3DArray.add_variables([
        ('new_property_1', {'dtype': np.float32,
                            'units': '1',
                            'default': 0}),
        ('new_property_2', {'dtype': np.float32,
                            'units': '1',
                            'default': 0})])


class ModelTemplate(OceanDrift):
    """Open source trajectory model based on the OpenDrift framework.

    This is a template model, to be copied and modified according to need.
    This class is a subclass of "OceanDrift", and thus inherits all its
    functionality. More specialised models may also be subclassed.
    """

    # Specify which element/particle type (class) to be used by this model.
    # A model can have/use only one element type
    ElementType = TemplateElementType

    # Specify which environment variables (e.g. wind, waves, currents...)
    # are needed/required by the present model, to be used for updating
    # the element properties (including propagation).
    # The convention of using CF standard_name (https://cfconventions.org)
    # allows reusing readers for many/all models.
    # Vector quantities should be referred to by the x-y components
    # (e.g. x_wind) and not the lon-lat componenets (e.g. east_wind).
    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'x_wind': {'fallback': 0, 'important': False},
        'y_wind': {'fallback': 0, 'important': False},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'important': False, 'profiles': True},
        'land_binary_mask': {'fallback': None},
        }
    # Default (fallback) numerical values will be used for particles 
    # outside of the spatial or tempoaral domain of all the available readers,
    # or if no readers are given for the specific variables.
    # If the fallback value is None, this means that simulation will
    # not run unless a reader has been added to provide this variable.
    # The attribute *important' relates to the concept of Lazy readers:
    # https://opendrift.github.io/tutorial.html#lazy-readers
    # Lazy readers will be initialized successively until there are active
    # readers for all variables, except for variables which are
    # marked as not important, as for wind in the example above.
    # If the attribute 'profiles' is True for a variable (as for
    # ocean_vertical_diffusivity in the example above), a vertical profile
    # of this variable is obtained over the vertical depth interval as 
    # defined by config setting drift:profiles_depth. These profiles are used for the vertical turbulence scheme:
    # https://opendrift.github.io/_modules/opendrift/models/oceandrift.html#OceanDrift.vertical_mixing
    # Profiles are typically required for ocean_vertical_diffusivity, and sometimes for 
    # salinity and temperature in order to calculate the stratification.
    
    def __init__(self, *args, **kwargs):

        # Here any preprocessing and initialisation of a simulation
        # may be performed, e.g. importing external data from file.

        # The constructor of parent class must always be called
        # to perform some necessary common initialisation tasks:
        super(ModelTemplate, self).__init__(*args, **kwargs)
            
        # The above initialisation includes adding common config options.
        # Further config options specific to this subclass can be
        # added with the following method (example options below from OpenOil):
        self._add_config({
            'seed:m3_per_hour': {'type': 'float', 'default': 1,
                'min': 0, 'max': 1e10, 'units': 'm3 per hour',
                'description': 'The amount (volume) of oil released per hour (or total amount if release is instantaneous)',
                'level': self.CONFIG_LEVEL_ESSENTIAL},
            'seed:droplet_diameter_min_subsea': {'type': 'float', 'default': 0.0005,
                'min': 1e-8, 'max': 1, 'units': 'meters',
                'description': 'The minimum dimaeter of oil droplet for a subsea release.',
                'level': self.CONFIG_LEVEL_BASIC},
            'processes:dispersion': {'type': 'bool', 'default': True,
                'description': 'Oil is removed from simulation (dispersed), if entrained as very small droplets.',
                'level': self.CONFIG_LEVEL_BASIC},
            'wave_entrainment:droplet_size_distribution': {'type': 'enum', 'enum': ['Johansen et al. (2015)', 'Li et al. (2017)'], 'default': 'Johansen et al. (2015)',
                'level': self.CONFIG_LEVEL_ADVANCED, 'description':
                'Algorithm to be used for calculating oil droplet size spectrum after entrainment by breaking waves.'},
           })
        
        # The _add_config method takes a dictionary where keys are config settings,
        # which can be adjusted by the user by:
        # >>> o.set_config(key, value)
        # and whose actual value is obtained during the simulation (*update* method below) by
        # value = o.get_config(key)
        # The config keys have a colon (:) separator for logical grouping into categories.
        # The config specification above has the following items:
        # - type: one of ['float', 'int', 'bool', 'enum']
        # - default: the default value
        # - units: the unit, used for labeling (potentially later for automatic scaling)
        # - description: a human readable help text, e.g. to be displayed in the GUI
        # - level:
        #   - CONFIG_LEVEL_ESSENTIAL (1): for config settings which the user must consider.
        #        These config items are e.g. shown on the front/seeding page of a GUI
        #   - CONFIG_LEVEL_BASIC (2): for config settings which most users may consider to change
        #   - CONFIG_LEVEL_ADVANCED (3): for config setting which only advanced users might change
        #
        # - If 'type' is 'float' or 'int' there are two additional items:
        #    - 'min' and 'max' which are the minimum and maximum allowed values of this setting.
        #       The set_config() method will raise an error of provided values are outside of this range.
        # - If 'type' is 'enum' there is an additional item 'enum' which is a list
        #     of the allowed values. These may be strings or numbers.
        #
        # Dafault values of config-settings which are inherited from parent models
        # may be overriden for the current module by:
        self._set_config_default('drift:vertical_mixing', True)

    def update(self):
        """Update positions and properties of particles."""

        # This function (``update``, may not be renamed) defines
        # the actions (processes) to be performed at each time step.
        # The default timestep is one hour, but this may be overriden
        # by the user.
        #
        # The elements and corresponding environment information
        # (wind, waves, current...) is available to this function as
        # the recarrays ``elements`` and ``environment``, respectively.
        # Element properties and environvent variables are named arrays
        # of same length (number of active particles), or scalars.
        #
        # E.g.:
        #    self.elements.z
        #    self.elements.mass
        #    self.environment.x_wind
        #    self.environment.y_wind
        #
        # E.g. to specify evaporation (mass reduction) proportional to wind speed:
        #
        #    self.elements.mass = self.elements.mass - \
        #        np.sqrt(self.environment.x_wind**2 + \
        #                self.environment.y_wind**2)*.01
        #
        # If particles need to be deactived (e.g. if
        # self.environemtn.someproperty is 0), the special function
        # `self.deactivate_elements()` should be used:
        # `self.deactivate_elements(self.environment.someproperty==0)`
        # 
        # Interaction with coast is handled by ``basemodel`` and is
        # chosen with configutation mechanism available to each module.
        #
        # Horizontal position of the particles (lon, lat) are the
        # only properties which should not be modified directly/freely.
        # Due to the need to correct for the curvature of the coordinate
        # systems, and taking into account vector orientations, a dedicated 
        # method ``self.update_positions()`` is
        # available using x- and y-velocity as arguments:
        # E.g. to simply move particles with ambient current:
        #
        # self.update_positions(self.environment.x_sea_water_velocity,
        #                       self.environment.y_sea_water_velocity)``
        #
        # However, instead of using this method directly, most
        # OpenDrift modules use a set of ready-made and reusable
        # advection methods, added below. 

        # Advection with horizontal ocean currents
        self.advect_ocean_current()

        # Advection with vertical ocean currents
        self.vertical_advection()

        # Additional advection due to wind shear in upper 10cm of ocean:
        # ``wind_drift_factor`` multiplied by wind vector
        self.advect_wind()

        # If Stokes drift is not available from a reader,
        # it is parameterized by wind
        self.stokes_drift()

        # Vertical turbulent mixing, taking into account:
        # - element property ``terminal_velocity`` (buoyancy)
        # - environment property ``ocean_vertical_diffusivity``
        # - interaction with surface through method ``surface_wave_mixing``  (should be renamed)
        # - interaction with seafloor through method ``bottom_interaction``
        self.vertical_mixing()

    def prepare_run(self):
        """Code to be run before a simulation loop (``update()``)"""
        # This method must also call the corresponding method of parent class

    def bottom_interaction(self, seafloor_depth):
        """Sub method of vertical_mixing, determines settling"""
