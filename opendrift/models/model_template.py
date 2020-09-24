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
See e.g. module windblow.py for a minimalistic example.

For more help, see https://github.com/OpenDrift/opendrift/wiki
or contact Knut-Frode Dagestad (knutfd@met.no)
For usage of the new model, see example scripts (e.g. example.py)

First we define an element type to be used by the trajectory model.
An element type (class) may also be imported/reused from
another module - see examples in the folder "elements"

"""

import numpy as np
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift import Lagrangian3DArray

class TemplateElementType(Lagrangian3DArray):
    """
    Extending LagrangianArray with relevant properties.

    Define and name the properties which the elements shall have.
    These are added to the inherited properties (e.g. from OceanDrift)
    Property names may be freely chosen, but the "update" function (below)
    will refer to the names in the specifications of the processes.
    """
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
    # Presently, a model can have/use only one element type
    ElementType = TemplateElementType

    # Specify which environment variables (e.g. wind, waves, currents...)
    # are needed/required by the present model, to be used for updating
    # the element properties (including propagation).
    # The convention of using CF standard_name (https://cfconventions.org)
    # allows reusing readers for many/all models.
    # Vector quantities should be referred to by the x-y components
    # (e.g. x_wind) and not the lon-lat componenets (e.g. east_wind).
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'x_wind', 'y_wind', 'land_binary_mask']

    # Default (fallback) numerical values may be specified.
    # These values will be used for particles outside of the
    # spatial or tempoaral domain of all the available readers,
    # or if no readers are given for the specific variables.
    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'x_wind': 0, 'y_wind': 0}

    def __init__(self, *args, **kwargs):
        """

        Here any preprocessing and initialisation of a simulation
        may be performed, e.g. importing external data from file.

        .. code::

            print('preprocessing...')

        Constructor of parent class must always be called
        to perform some necessary common initialisation tasks

        .. code::
            super(ModelTemplate, self).__init__(*args, **kwargs)
        """
        super(ModelTemplate, self).__init__(*args, **kwargs)

    def update(self):
        """Update positions and properties of particles.

        This function (``update``, may not be renamed) defines
        the actions (processes) to be performed at each time step.
        The default timestep is one hour, but this may be overriden
        by the user.

        The elements and corresponding environment information
        (wind, waves, current...) is available to this function as
        the recarrays ``elements`` and ``environment``, respectively.
        Element properties and environvent variables are named arrays
        of same length (number of active particles), or scalars.

        E.g.:

        .. code::

            self.elements.z
            self.elements.mass
            self.environment.x_wind
            self.environment.y_wind

        E.g. to specify evaporation (mass reduction) proportional to wind speed:

        .. code::

            self.elements.mass = self.elements.mass - \
                np.sqrt(self.environment.x_wind**2 + \
                        self.environment.y_wind**2)*.01

        If particles need to be deactived (e.g. if
        self.environemtn.someproperty is 0), the special function
        `self.deactivate_elements()` should be used:
        `self.deactivate_elements(self.environment.someproperty==0)`

        Interaction with coast is handled by ``basemodel`` and is
        chosen with configutation mechanism available to each module.

        Horizontal position of the particles (lon, lat) are the
        only properties which should not be modified directly/freely.
        Due to the need to correct for the curvature of the coordinate
        systems, and taking into account vector orientations, a dedicated 
        method ``self.update_positions()`` is
        available using x- and y-velocity as arguments:
        E.g. to simply move particles with ambient current:
        ``
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)``

        However, instead of using this method directly, most
        OpenDrift modules use a set of ready-made and reusable
        advection methods, added below. 

        """
       
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

    def bottom_interaction(self, seafloor_depth):
        """Sub method of vertical_mixing, determines settling"""
