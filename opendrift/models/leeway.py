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
"""
Leeway is the search and rescue (SAR) model developed by the US Coast Guard, as originally described in

    Allen, A.A, 2005: Leeway Divergence, USCG R&D Center Technical Report CG-D-05-05. Available through http://www.ntis.gov, reference ADA435435

    Allen A.A. and J.V. Plourde (1999) Review of Leeway; Field Experiments and Implementation, USCG R&D Center Technical Report CG-D-08-99. Available through http://www.ntis.gov, reference ADA366414

and later extended and modified by e.g.

    Breivik, O., A. Allen, C. Maisondieu, J.-C. Roth, and B. Forest, 2012: The leeway of shipping containers at different immersion levels. Ocean Dyn., 62, 741–752, doi:10.1007/s10236-012-0522-z

The Leeway model is based on empirically determined coefficients as tabulated in https://github.com/OpenDrift/opendrift/blob/master/opendrift/models/OBJECTPROP.DAT

The Leeway model is been reprogrammed in Python for OpenDrift by Knut-Frode Dagestad of the Norwegian Meteorological Institute.
"""

from builtins import range
import os
from collections import OrderedDict
import logging

logger = logging.getLogger(__name__)

import numpy as np

from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

RIGHT = 0
LEFT = 1


# Defining the leeway element properties
class LeewayObj(LagrangianArray):
    """Extending LagrangianArray with variables relevant for leeway objects.

    """
    variables = LagrangianArray.add_variables([
        ('object_type', {
            'dtype': np.uint16,
            'units': '1',
            'seed': False,
            'default': 0
        }),
        ('orientation', {
            'dtype': np.uint8,
            'units': '1',
            'description':
            '0/1 is left/right of downwind. Randomly chosen at seed time',
            'seed': False,
            'default': 1
        }),
        ('jibe_probability', {
            'dtype': np.float32,
            'units': '1/h',
            'description':
            'Probability per hour that an object may change orientation (jibing)',
            'default': 0.04
        }),
         ('capsized', {
            'dtype': np.uint8,
            'units': '1',
            'description': '0 is not capsized, changed to 1 after capsizing (irreversible). After capsizing, leeway coeffieiencts are reduced as given by config item capsized:leeway_fraction',
            'seed': True,
            'default': 0
        }),
        ('downwind_slope', {
            'dtype': np.float32,
            'units': '%',
            'seed': False,
            'default': 1
        }),
        ('crosswind_slope', {
            'dtype': np.float32,
            'units': '1',
            'seed': False,
            'default': 1
        }),
        ('downwind_offset', {
            'dtype': np.float32,
            'units': 'cm/s',
            'seed': False,
            'default': 0
        }),
        ('crosswind_offset', {
            'dtype': np.float32,
            'units': 'cm/s',
            'seed': False,
            'default': 0
        }),
        ('downwind_eps', {
            'dtype': np.float32,
            'units': 'cm/s',
            'seed': False,
            'default': 0
        }),
        ('crosswind_eps', {
            'dtype': np.float32,
            'units': 'cm/s',
            'seed': False,
            'default': 0
        }),
        ('current_drift_factor', {
            'dtype':
            np.float32,
            'units':
            '1',
            'description':
            'Elements are moved with this fraction of the '
            'current vector, in addition to currents '
            'and Stokes drift',
            'default':
            1
        })
    ])


class Leeway(OpenDriftSimulation):
    """The Leeway model in the OpenDrift framework.

        Advects a particle (a drifting object) with the ambient current and
        as a function of the wind vector according to the drift properties
        of the object.
    """

    ElementType = LeewayObj

    required_variables = {
        'x_wind': {
            'fallback': None
        },
        'y_wind': {
            'fallback': None
        },
        'x_sea_water_velocity': {
            'fallback': None
        },
        'y_sea_water_velocity': {
            'fallback': None
        },
        'land_binary_mask': {
            'fallback': None
        },
    }

    # Default colors for plotting
    status_colors = {
        'initial': 'green',
        'active': 'blue',
        'missing_data': 'gray',
        'stranded': 'red',
        'evaporated': 'yellow',
        'dispersed': 'magenta'
    }


    # Configuration
    def __init__(self, d=None, *args, **kwargs):

        # Read leeway object properties from file
        if d is None:
            d = os.path.dirname(os.path.realpath(__file__))
            objprop_file = open(d + '/OBJECTPROP.DAT', 'r')
        else:
            objprop_file = open(d, 'r')

        # objprop_file=open('/home/oyvindb/Pro/opendrift/OBJECTPROP.DAT','r')
        objproptxt = objprop_file.readlines()
        objprop_file.close()
        self.leewayprop = OrderedDict({})
        for i in range(len(objproptxt) // 3 + 1):
            # Stop at first blank line
            if not objproptxt[i * 3].strip():
                break
            elems = objproptxt[i * 3].split()[0]

            objKey = objproptxt[i * 3].split()[0].strip()
            arr = [float(x) for x in objproptxt[i * 3 + 2].split()]
            props = {'OBJKEY': objKey}
            props['Description'] = objproptxt[i * 3 + 1].strip()
            props['DWSLOPE'] = arr[0]
            props['DWOFFSET'] = arr[1]
            props['DWSTD'] = arr[2]
            props['CWRSLOPE'] = arr[3]
            props['CWROFFSET'] = arr[4]
            props['CWRSTD'] = arr[5]
            props['CWLSLOPE'] = arr[6]
            props['CWLOFFSET'] = arr[7]
            props['CWLSTD'] = arr[8]
            self.leewayprop[i + 1] = props

        # Config
        descriptions = [
            self.leewayprop[p]['Description'] for p in self.leewayprop
        ]

        # Calling general constructor of parent class
        super(Leeway, self).__init__(*args, **kwargs)

        self._add_config({
            'seed:object_type': {
                'type': 'enum',
                'enum': descriptions,
                'default': descriptions[0],
                'description': 'Leeway object category for this simulation',
                'level': CONFIG_LEVEL_ESSENTIAL
            },
            'seed:jibe_probability': {
                'type': 'float',
                'default': 0.04,
                'min': 0,
                'max': 1,
                'description':
                'Probability per hour for jibing (objects changing orientation)',
                'units': 'probability',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:capsizing': {
                'type': 'bool',
                'default': False,
                'description':
                'If True, elements can be capsized when wind exceeds threshold given by config item capsize:wind_threshold',
                'level': CONFIG_LEVEL_BASIC
            },
             'capsizing:leeway_fraction': {
                'type': 'float',
                'default': 0.4,
                'min': 0,
                'max': 1,
                'description':
                'After capsizing, leeway coefficients are reduced by multiplying by this factor',
                'units': 'fraction',
                'level': CONFIG_LEVEL_BASIC
            },
              'capsizing:wind_threshold': {
                'type': 'float',
                'default': 30,
                'min': 0,
                'max': 50,
                'description':
                'Probability of capsizing per hour is: 0.5 + 0.5tanh((windspeed-wind_threshold)/wind_threshold_sigma)',
                'units': 'm/s',
                'level': CONFIG_LEVEL_BASIC
            },
               'capsizing:wind_threshold_sigma': {
                'type': 'float',
                'default': 5,
                'min': 0,
                'max': 20,
                'description':
                'Sigma parameter in parameterization of capsize probability',
                'units': 'm/s',
                'level': CONFIG_LEVEL_BASIC
            },
        })

        self._set_config_default('general:time_step_minutes', 10)
        self._set_config_default('general:time_step_output_minutes', 60)
        self._set_config_default('drift:max_speed', 5)

    def seed_elements(self, lon, lat, object_type=None, **kwargs):
        """Seed particles in a cone-shaped area over a time period."""
        # All particles carry their own object_type (number),
        # but so far we only use one for each sim
        # objtype = np.ones(number)*object_type

        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        if 'number' in kwargs and kwargs['number'] is not None:
            number = kwargs['number']
        elif len(lon) > 1:
            number = len(lon)
        else:
            number = self.get_config('seed:number')

        if object_type is None:
            object_name = self.get_config('seed:object_type')
            # Get number from name
            found = False
            for object_type in range(1, len(self.leewayprop) + 1):
                if self.leewayprop[object_type]['OBJKEY'] == object_name or (
                        self.leewayprop[object_type]['Description']
                        == object_name):
                    found = True
                    break
            if found is False:
                logger.info(self.list_configspec())
                raise ValueError('Object %s not available' % object_type)

        logger.info('Seeding elements of object type %i: %s (%s)' %
                    (object_type, self.leewayprop[object_type]['OBJKEY'],
                     self.leewayprop[object_type]['Description']))

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
        downwind_slope = ones * self.leewayprop[object_type]['DWSLOPE']
        downwind_offset = ones * self.leewayprop[object_type]['DWOFFSET']
        dwstd = self.leewayprop[object_type]['DWSTD']
        rdw = np.zeros(number)
        epsdw = np.zeros(number)
        for i in range(number):
            rdw[i] = np.random.randn(1)[0]
            epsdw[i] = rdw[i] * dwstd
            # Avoid negative downwind slopes
            while downwind_slope[i] + epsdw[i] / 20.0 < 0.0:
                rdw[i] = np.random.randn(1)
                epsdw[i] = rdw[i] * dwstd
        downwind_eps = epsdw
        # NB
        # downwind_eps = np.zeros(number)

        # Crosswind leeway properties
        rcw = np.random.randn(number)
        crosswind_slope = np.zeros(number)
        crosswind_offset = np.zeros(number)
        crosswind_eps = np.zeros(number)
        crosswind_slope[orientation == RIGHT] = \
            self.leewayprop[object_type]['CWRSLOPE']
        crosswind_slope[orientation == LEFT] = \
            self.leewayprop[object_type]['CWLSLOPE']
        crosswind_offset[orientation == RIGHT] = \
            self.leewayprop[object_type]['CWROFFSET']
        crosswind_offset[orientation == LEFT] = \
            self.leewayprop[object_type]['CWLOFFSET']
        crosswind_eps[orientation == RIGHT] = \
            rcw[orientation == RIGHT] * \
            self.leewayprop[object_type]['CWRSTD']
        crosswind_eps[orientation == LEFT] = \
            rcw[orientation == LEFT] * \
            self.leewayprop[object_type]['CWLSTD']

        # NB
        # crosswind_eps = np.zeros(number)

        # Store seed data for ASCII format output
        if hasattr(self, 'seed_cone_arguments'):
            self.ascii = self.seed_cone_arguments
        else:
            self.ascii = {
                'lon': lon,
                'lat': lat,
                'radius': kwargs['radius'] if 'radius' in kwargs else 0,
                'number': number,
                'time': kwargs['time']
            }

        # Call general seed_elements function of OpenDriftSimulation class
        # with the specific values calculated
        super(Leeway, self).seed_elements(lon,
                                          lat,
                                          orientation=orientation,
                                          object_type=object_type,
                                          downwind_slope=downwind_slope,
                                          crosswind_slope=crosswind_slope,
                                          downwind_offset=downwind_offset,
                                          crosswind_offset=crosswind_offset,
                                          downwind_eps=downwind_eps,
                                          crosswind_eps=crosswind_eps,
                                          **kwargs)

    def list_object_categories(self, substr=None):
        '''Display leeway categories to screen

        Print only objects containing 'substr', if specified'''

        for i, p in enumerate(self.leewayprop):
            description = self.leewayprop[p]['Description']
            objkey = self.leewayprop[p]['OBJKEY']
            if substr is not None:
                if substr.lower() not in description.lower() + objkey.lower():
                    continue
            print('%i %s %s' % (i + 1, objkey, description))

    def plot_capsize_probability(self):
        U = np.linspace(0, 35, 100)
        wind_threshold = self.get_config('capsizing:wind_threshold')
        sigma = self.get_config('capsizing:wind_threshold_sigma')
        p = self.capsize_probability(U, wind_threshold, sigma)
        import matplotlib.pyplot as plt
        plt.plot(U, p)
        plt.title(f'p(u) = 0.5 + 0.5*tanh((u - {wind_threshold} / {sigma})')
        plt.xlabel('Wind speed  [m/s]')
        plt.ylabel('Probability of capsizing per hour')
        plt.show()

    def capsize_probability(self, wind, threshold, sigma):
        return .5 + .5*np.tanh((wind-threshold)/sigma)

    def update(self):
        """Update positions and properties of leeway particles."""

        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)
        # CCC update wind direction
        winddir = np.arctan2(self.environment.x_wind, self.environment.y_wind)

        # Capsizing
        if self.get_config('processes:capsizing') is True:
            wind_threshold = self.get_config('capsizing:wind_threshold')
            wind_threshold_sigma = self.get_config('capsizing:wind_threshold_sigma')
            # For forward run, elements can be capsized, but for backwards run, only capsized elements can be un-capsized
            if self.simulation_direction() == 1:  # forward run
                can_be_capsized = np.where(self.elements.capsized==0)[0]
            else:
                can_be_capsized = np.where(self.elements.capsized==1)[0]
            if len(can_be_capsized) > 0:
                probability = self.capsize_probability(windspeed[can_be_capsized],
                        wind_threshold, wind_threshold_sigma)*np.abs(self.time_step.total_seconds())/3600
                # NB: assuming small timestep
                to_be_capsized = np.where(np.random.rand(len(can_be_capsized)) < probability)[0]
                to_be_capsized = can_be_capsized[to_be_capsized]
                logger.warning(f'Capsizing {len(to_be_capsized)} of {len(can_be_capsized)} elements')
                self.elements.capsized[to_be_capsized] = 1 - self.elements.capsized[to_be_capsized]

        # Move particles with the leeway CCC TODO
        downwind_leeway = (
            (self.elements.downwind_slope + self.elements.downwind_eps / 20.0)
            * windspeed + self.elements.downwind_offset +
            self.elements.downwind_eps / 2.0) * .01  # In m/s
        crosswind_leeway = ((self.elements.crosswind_slope +
                             self.elements.crosswind_eps / 20.0) * windspeed +
                            self.elements.crosswind_offset +
                            self.elements.crosswind_eps / 2.0) * .01  # In m/s
        sinth = np.sin(winddir)
        costh = np.cos(winddir)
        y_leeway = downwind_leeway * costh + crosswind_leeway * sinth
        x_leeway = -downwind_leeway * sinth + crosswind_leeway * costh
        capsize_fraction = self.get_config('capsizing:leeway_fraction')  # Reducing leeway for capsized elements
        x_leeway[self.elements.capsized==1] *= capsize_fraction
        y_leeway[self.elements.capsized==1] *= capsize_fraction
        self.update_positions(-x_leeway, y_leeway)

        # Move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Jibe elements randomly according to given probability
        jibe_rate = -np.log(1 - self.elements.jibe_probability
                            ) / 3600  # Hourly to instantaneous
        jp_per_timestep = 1 - np.exp(
            -jibe_rate * np.abs(self.time_step.total_seconds()))
        jib = jp_per_timestep > np.random.random(self.num_elements_active())
        self.elements.crosswind_slope[
            jib] = -self.elements.crosswind_slope[jib]
        self.elements.orientation[jib] = 1 - self.elements.orientation[jib]
        logger.debug('Jibing %i out of %i elements.' %
                     (np.sum(jib), self.num_elements_active()))

    def export_ascii(self, filename):
        '''Export output to ASCII format of original version'''

        try:
            f = open(filename, 'w')
        except:
            raise ValueError('Could not open file for writing: ' + filename)

        for inp in ['lon', 'lat', 'radius', 'time']:
            if len(np.atleast_1d(self.ascii[inp])) == 1:
                self.ascii[inp] = [self.ascii[inp], self.ascii[inp]]
        f.write('# Drift simulation initiated [UTC]:\n')
        f.write('simDate simTime\n')
        f.write(self.start_time.strftime('%Y-%m-%d\t%H:%M:%S\t#\n'))
        f.write('# Model version:\n'
                'modelVersion\n'
                ' 3.00\n'  # OpenDrift version
                '# Object class id & name:\n'
                'objectClassId  objectClassName\n')
        try:
            objtype = self.elements.object_type[0]
        except:
            objtype = self.elements_deactivated.object_type[0]
        f.write(' %i\t%s\n' % (objtype, self.leewayprop[objtype]['OBJKEY']))
        f.write('# Seeding start time, position & radius:\n'
                'startDate\tstartTime\tstartLon\tstartLat\tstartRad\n')
        f.write('%s\t%s\t%s\t%s\n' %
                (self.ascii['time'][0], self.ascii['lon'][0],
                 self.ascii['lat'][0], self.ascii['radius'][0] / 1000.))
        f.write('# Seeding end time, position & radius:\n'
                'endDate\tendTime\tendLon\tendLat\tendRad\n')
        f.write('%s\t%s\t%s\t%s\n' %
                (self.ascii['time'][1], self.ascii['lon'][1],
                 self.ascii['lat'][1], self.ascii['radius'][1] / 1000.))
        seedDuration = (self.ascii['time'][1] -
                        self.ascii['time'][0]).total_seconds() / 60.
        seedSteps = seedDuration / (self.time_step_output.total_seconds() /
                                    60.)
        seedSteps = np.maximum(1, seedSteps)
        f.write('# Duration of seeding [min] & [timesteps]:\n'
                'seedDuration   seedSteps\n'
                '    %i      %i\n'
                '# Length of timestep [min]:\n'
                'timeStep\n' % (seedDuration, seedSteps))
        f.write('%i\n' % (self.time_step_output.total_seconds() / 60.))
        f.write('# Length of model simulation [min] & [timesteps]:\n'
                'simLength  simSteps\n')
        f.write('%i\t%i\n' % ((self.time - self.start_time).total_seconds() /
                              60., self.steps_output))

        f.write('# Total no of seeded particles:\n'
                'seedTotal\n'
                ' %s\n' % (self.num_elements_activated()))
        f.write('# Particles seeded per timestep:\n'
                'seedRate\n'
                ' %i\n' % (self.num_elements_activated() / seedSteps))

        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()

        beforeseeded = 0
        lons, statuss = self.get_property('lon')
        lats, statuss = self.get_property('lat')
        orientations, statuss = self.get_property('orientation')
        for step in range(self.steps_output):
            lon = lons[step, :]
            lat = lats[step, :]
            orientation = orientations[step, :]
            status = statuss[step, :]
            if sum(status ==
                   0) == 0:  # All elements deactivated: using last position
                lon = lons[step - 1, :]
                lat = lats[step - 1, :]
                orientation = orientations[step - 1, :]
                status = statuss[step - 1, :]
            num_active = np.sum(~status.mask)
            status[status.mask] = 41  # seeded on land
            lon[status.mask] = 0
            lat[status.mask] = 0
            orientation[status.mask] = 0
            status[status == 0] = 11  # active
            status[status == 1] = 41  # stranded
            f.write('\n# Date [UTC]:\nnowDate   nowTime\n')
            f.write((
                self.start_time +
                self.time_step_output * step).strftime('%Y-%m-%d\t%H:%M:%S\n'))
            f.write(
                '# Time passed [min] & [timesteps], now seeded, seeded so far:\ntimePassed  nStep   nowSeeded   nSeeded\n'
            )
            f.write('  %i\t%i\t%i\t%i\n' %
                    ((self.time_step_output * step).total_seconds() / 60,
                     step + 1, num_active - beforeseeded, num_active))
            beforeseeded = num_active
            f.write('# Mean position:\nmeanLon meanLat\n')
            f.write('%f\t%f\n' %
                    (np.mean(lon[status == 11]), np.mean(lat[status == 11])))
            f.write('# Particle data:\n')
            f.write('id  lon lat state   age orientation\n')
            age_minutes = self.time_step_output.total_seconds() * (
                step - index_of_first) / 60
            age_minutes[age_minutes < 0] = 0
            for i in range(num_active):
                f.write('%i\t%.6f\t%.6f\t%i\t%i\t%i\n' %
                        (i + 1, lon[i], lat[i], status[i], age_minutes[i],
                         orientation[i]))

        f.close()

    def _substance_name(self):
        # TODO: find a better algorithm to return name of object category
        if self.history is not None:
            object_type = self.history['object_type'][0, 0]
            if not np.isfinite(object_type):  # For backward simulations
                object_type = self.history['object_type'][-1, 0]
        else:
            object_type = np.atleast_1d(self.elements_scheduled.object_type)[0]
        if np.isfinite(object_type):
            return self.leewayprop[object_type]['OBJKEY']
        else:
            return ''
