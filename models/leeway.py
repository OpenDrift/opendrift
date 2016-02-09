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

import os
import logging
from datetime import datetime
from collections import OrderedDict

import numpy as np

from opendrift import OpenDriftSimulation
from elements import LagrangianArray
from readers.reader import pyproj, Reader

RIGHT = 0
LEFT = 1


# Defining the leeway element properties
class LeewayObj(LagrangianArray):
    """Extending LagrangianArray with variables relevant for leeway objects.

    """
    variables = LagrangianArray.add_variables([
        ('objectType', {'dtype': np.int16,
                        'unit': '1',
                        'default': 0}),
        ('orientation', {'dtype': np.int16,
                         'unit': '1',
                         'default': 1}),
        ('jibeProbability', {'dtype': np.float32,
                             'unit': '1/h',
                             'default': 0.04}),
        ('downwindSlope', {'dtype': np.float32,
                           'unit': '%',
                           'default': 1}),
        ('crosswindSlope', {'dtype': np.float32,
                            'unit': '1',
                            'default': 1}),
        ('downwindOffset', {'dtype': np.float32,
                            'unit': 'cm/s',
                            'default': 0}),
        ('crosswindOffset', {'dtype': np.float32,
                             'unit': 'cm/s',
                             'default': 0}),
        ('downwindEps', {'dtype': np.float32,
                         'unit': 'cm/s',
                         'default': 0}),
        ('crosswindEps', {'dtype': np.float32,
                          'unit': 'cm/s',
                          'default': 0})
        ])


class Leeway(OpenDriftSimulation):
    """The Leeway model in the OpenDrift framework.

        Advects a particle (a drifting object) with the ambient current and
        as a function of the wind vector according to the drift properties
        of the object.
    """

    ElementType = LeewayObj

    required_variables = ['x_wind', 'y_wind',
                          'x_sea_water_velocity', 'y_sea_water_velocity',
                          'land_binary_mask']

    fallback_values = {'x_wind': 0.0,
                       'y_wind': 0.0,
                       'x_sea_water_velocity': 0.0,
                       'y_sea_water_velocity': 0.0}

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'missing_data': 'gray', 'stranded': 'red',
                     'evaporated': 'yellow', 'dispersed': 'magenta'}

    # Configuration
    def __init__(self, d=None, *args, **kwargs):

        # Read leeway object properties from file
        if d is None:
            d = os.path.dirname(os.path.realpath(__file__))
            objprop_file = open(d + '/OBJECTPROP.DAT', 'r')
        else:
            objprop_file = open(d, 'r')

        #objprop_file=open('/home/oyvindb/Pro/opendrift/OBJECTPROP.DAT','r')
        objproptxt = objprop_file.readlines()
        objprop_file.close()
        self.leewayprop = OrderedDict({})
        for i in xrange(len(objproptxt)//3+1):
            # Stop at first blank line
            if not objproptxt[i*3].strip():
                break
            elems = objproptxt[i*3].split()[0]

            objKey = objproptxt[i*3].split()[0]
            arr = [float(x) for x in objproptxt[i*3+2].split()]
            props = {'OBJKEY': objKey}
            props['Description'] = objproptxt[i*3+1]
            props['DWSLOPE'] = arr[0]
            props['DWOFFSET'] = arr[1]
            props['DWSTD'] = arr[2]
            props['CWRSLOPE'] = arr[3]
            props['CWROFFSET'] = arr[4]
            props['CWRSTD'] = arr[5]
            props['CWLSLOPE'] = arr[6]
            props['CWLOFFSET'] = arr[7]
            props['CWLSTD'] = arr[8]
            self.leewayprop[i+1] = props

        # Calling general constructor of parent class
        super(Leeway, self).__init__(*args, **kwargs)

    #def seed_leeway(self, lon, lat, radius=0, number=1, time=None,
    def seed_elements(self, lon, lat, radius=0, number=1, time=None,
                      radius1=None, lon1=None, lat1=None, time1=None,
                      objectType=1):
        """Seed particles in a cone-shaped area over a time period."""
        # All particles carry their own objectType (number),
        # but so far we only use one for each sim
        # objtype = np.ones(number)*objectType

        # Probability of jibing (4 % per hour)
        pjibe = 0.04

        if time is None:
            # Use first time of first reader if time is not given for seeding
            try:
                for reader in self.readers.items():
                    if reader[1].start_time is not None:
                        firstReader = reader[1]
                        break
            except:
                raise ValueError('Time must be specified when no '
                                 'reader is added')
            logging.info('Using start time (%s) of reader %s' %
                         (firstReader.start_time, firstReader.name))
            self.start_time = firstReader.start_time
        else:
            self.start_time = time

        if radius1 is None:
            radius1 = radius

        if time1 is None:
            time1 = self.start_time

        if lon1 is None:
            lon1 = lon
            lat1 = lat

        # Drift orientation of particles.  0 is right of downwind,
        # 1 is left of downwind
        # orientation = 1*(np.random.rand(number)>=0.5)
        # Odd numbered particles are left-drifting, even are right of downwind.
        orientation = np.r_[:number] % 2
        ones = np.ones_like(orientation)

        # Downwind leeway properties.
        # Generate normal, N(0,1), random perturbations for leeway coeffs.
        # Negative downwind slope must be avoided as
        # particles should drift downwind.
        # The problem arises because of high error variances (see e.g. PIW-1).
        downwindSlope = ones*self.leewayprop[objectType]['DWSLOPE']
        downwindOffset = ones*self.leewayprop[objectType]['DWOFFSET']
        dwstd = self.leewayprop[objectType]['DWSTD']
        rdw = np.zeros(number)
        epsdw = np.zeros(number)
        for i in xrange(number):
            rdw[i] = np.random.randn(1)
            epsdw[i] = rdw[i]*dwstd
            # Avoid negative downwind slopes
            while downwindSlope[i] + epsdw[i]/20.0 < 0.0:
                rdw[i] = np.random.randn(1)
                epsdw[i] = rdw[i]*dwstd
        downwindEps = epsdw
        # NB
        #downwindEps = np.zeros(number)

        # Crosswind leeway properties
        rcw = np.random.randn(number)
        crosswindSlope = np.zeros(number)
        crosswindOffset = np.zeros(number)
        crosswindEps = np.zeros(number)
        crosswindSlope[orientation == RIGHT] = \
            self.leewayprop[objectType]['CWRSLOPE']
        crosswindSlope[orientation == LEFT] = \
            self.leewayprop[objectType]['CWLSLOPE']
        crosswindOffset[orientation == RIGHT] = \
            self.leewayprop[objectType]['CWROFFSET']
        crosswindOffset[orientation == LEFT] = \
            self.leewayprop[objectType]['CWLOFFSET']
        crosswindEps[orientation == RIGHT] = \
            rcw[orientation == RIGHT] * \
            self.leewayprop[objectType]['CWRSTD']
        crosswindEps[orientation == LEFT] = \
            rcw[orientation == LEFT] * \
            self.leewayprop[objectType]['CWLSTD']

        # NB
        #crosswindEps = np.zeros(number)

        # Jibe probability
        jibeProbability = ones*pjibe

        # Call general seed_elements function of OpenDriftSimulation class
        # with the specific values calculated
        super(Leeway, self).seed_elements(
            lon=lon, lat=lat, radius=radius,
            number=number, time=time, cone=True,
            orientation=orientation, objectType=objectType,
            downwindSlope=downwindSlope,
            crosswindSlope=crosswindSlope,
            downwindOffset=downwindOffset,
            crosswindOffset=crosswindOffset,
            downwindEps=downwindEps, crosswindEps=crosswindEps,
            jibeProbability=jibeProbability)

    def update(self):
        """Update positions and properties of leeway particles."""

        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)
        # CCC update wind direction
        winddir = np.arctan2(self.environment.x_wind, self.environment.y_wind)

        # Move particles with the leeway CCC TODO
        downwind_leeway = ((self.elements.downwindSlope +
                            self.elements.downwindEps/20.0)*windspeed +
                           self.elements.downwindOffset +
                           self.elements.downwindEps/2.0)*.01  # In m/s
        crosswind_leeway = ((self.elements.crosswindSlope +
                            self.elements.crosswindEps/20.0)*windspeed +
                            self.elements.crosswindOffset +
                            self.elements.crosswindEps/2.0)*.01  # In m/s
        sinth = np.sin(winddir)
        costh = np.cos(winddir)
        y_leeway = downwind_leeway*costh+crosswind_leeway*sinth
        x_leeway = -downwind_leeway*sinth+crosswind_leeway*costh
        self.update_positions(-x_leeway, y_leeway)

        # Deactivate elements on land
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

        # Move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Jibe elements randomly according to given probability
        jp_per_timestep = self.elements.jibeProbability * \
            self.time_step.total_seconds() / 3600.0
        jib = jp_per_timestep > np.random.random(self.num_elements_active())
        self.elements.crosswindSlope[jib] = - self.elements.crosswindSlope[jib]
        logging.debug('Jibing %i out of %i elements.' %
                      (sum(jib), self.num_elements_active()))
