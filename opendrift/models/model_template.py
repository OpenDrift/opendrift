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

import numpy as np
from opendrift.elements import LagrangianArray
from opendrift.models.basemodel import OpenDriftSimulation


# This file is a template for an ocean trajectory model
# to be implemented within the OpenDrift framework.
# To create a new model, make a copy (renamed) of this file, and
# modify the class names and details below, according to the comments.
# See e.g. module windblow.py for a minimalistic example.
#
# For more help, see https://github.com/OpenDrift/opendrift/wiki
# or contact Knut-Frode Dagestad (knutfd@met.no)
# For usage of the new model, see example scripts (e.g. example.py)


# First we define an element type to be used by the trajectory model.
# An element type (class) may also be imported/reused from
# another module - see examples in the folder "elements"
class TemplateElementType(LagrangianArray):
    """Extending LagrangianArray with relevant properties."""
    # Define and name the properties which the elements shall have.
    # These are added to the four core properties (ID, lon, lat, z)
    # inherited from the basic/generic LagrangianArray class.
    # Property names may be freely chosen, but the "update" function (below)
    # will refer to the names in the specifications of the processes.
    variables = LagrangianArray.add_variables([
        ('mass_oil', {'dtype': np.float32,
                      'unit': 'kg'}),
        ('mass_evaporated', {'dtype': np.float32,
                             'unit': 'kg',
                             'default': 0}),
        ('mass_emulsion', {'dtype': np.float32,
                           'unit': 'kg',
                           'default': 0})])


class ModelTemplate(OpenDriftSimulation):
    """Open source trajectory model based on the OpenDrift framework.

    This is a template model, to be copied and modified according to need.
    This class is a subclass of "OpenDriftSimulation", and thus inherits
    all functionality. More specialised models may also be subclassed.
    """

    # Specify which element/particle type (class) to be used by this model.
    # Presently, a model can have/use only one element type
    ElementType = TemplateElementType

    # Specify which environment variables (e.g. wind, waves, currents...)
    # are needed/required by the present model, to be used for updating
    # the element properties (including propagation).
    # The convention of using CF standard_name (http://cfconventions.org)
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

        # Here any preprocessing may be performed,
        # e.g. importing external data from file.
        print('preprocessing...')

        # Constructor of parent class must always be called
        # to perform some necessary common initialisation tasks
        super(ModelTemplate, self).__init__(*args, **kwargs)

    def update(self):
        """Update positions and properties of particles."""

        # This function ("update", may not be renamed) defines
        # the actions (processes) to be performed at each time step.
        # The default timestep is one hour, but this may be overriden
        # by the user (setting model.time_step, see e.g. example.py)

        # The elements and corresponding environment information
        # (wind, waves, current...) is available to this function as
        # the recarrays "elements" and "environment", respectively.
        # Element properties and environvent variables are namedarrays
        # of same length (number of active particles), or scalars.
        # E.g.:
        # self.elements.z
        # self.elements.mass
        # self.environment.x_wind
        # self.environment.y_wind
        # E.g. to specify evaporation (mass reduction)
        #   proportional to wind speed:
        # self.elements.mass = self.elements.mass - \
        #       np.sqrt(self.environment.x_wind**2 + \
        #               self.environment.y_wind**2)*.01

        # Horizontal position of the particles (lon, lat) are the
        # only properties which should not be modified directly/freely.
        # Due to the need to correct for the curvature of the coordinate
        # systems, a dedicated function (self.update_positions() should be
        # used instead, with x- and y-velocity as arguments:
        # E.g. to simply move particles with ambient current:
        # self.update_positions(self.environment.x_sea_water_velocity,
        #                       self.environment.y_sea_water_velocity)
        # E.g. to move particles with 3 percent of the wind speed:
        # wind_factor = 0.03
        # self.update_positions(self.environment.x_wind*wind_factor,
        #                      self.environment.y_wind*wind_factor)

        # If particles need to be deactived (e.g. when hitting land),
        # the special function self.deactivate_elements() should be used:
        # self.deactivate_elements(self.environment.land_binary_mask==1)
