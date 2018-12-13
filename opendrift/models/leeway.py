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

from builtins import range
import os
import logging
from collections import OrderedDict

import numpy as np

from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray

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

    #fallback_values = {'x_wind': 0.0,
    #                   'y_wind': 0.0,
    #                   'x_sea_water_velocity': 0.0,
    #                   'y_sea_water_velocity': 0.0}

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'missing_data': 'gray', 'stranded': 'red',
                     'evaporated': 'yellow', 'dispersed': 'magenta'}

    max_speed = 1.5  # Assumed max average speed of any element

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
        for i in range(len(objproptxt)//3+1):
            # Stop at first blank line
            if not objproptxt[i*3].strip():
                break
            elems = objproptxt[i*3].split()[0]

            objKey = objproptxt[i*3].split()[0].strip()
            arr = [float(x) for x in objproptxt[i*3+2].split()]
            props = {'OBJKEY': objKey}
            props['Description'] = objproptxt[i*3+1].strip()
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

        # Config
        descriptions = [self.leewayprop[p]['Description'] for
                        p in self.leewayprop]
        self._add_config('seed:object_type', descriptions, 'Object type',
                         overwrite=True)
        # Default jibe probability is 4% per hour
        self._add_config('seed:jibeProbability',
                         'float(min=0, max=1, default=0.04)',
                         'Jibe probability', overwrite=True)

        # Calling general constructor of parent class
        super(Leeway, self).__init__(*args, **kwargs)

        self._add_config('general:time_step_minutes',
                         'integer(min=1, max=1440, default=10)',
                         'Time step in minutes',
                         overwrite=True)
        self._add_config('general:time_step_output_minutes',
                         'integer(min=1, max=1440, default=60)',
                         'Output time step in minutes',
                         overwrite=True)

    def seed_elements(self, lon, lat, radius=0, number=1, time=None,
                      objectType=None, cone=None, jibeProbability=None):
        """Seed particles in a cone-shaped area over a time period."""
        # All particles carry their own objectType (number),
        # but so far we only use one for each sim
        # objtype = np.ones(number)*objectType
        # Note: cone is not used, simply to provide same interface as others

        if objectType is None:
            object_name = self.get_config('seed:object_type')
            # Get number from name
            found = False
            for objectType in range(1, len(self.leewayprop)+1):
                if self.leewayprop[objectType]['OBJKEY'] == object_name or(
                    self.leewayprop[objectType]['Description'] == object_name):
                    found = True
                    break
            if found is False:
                logging.info(self.list_configspec())
                raise ValueError('Object %s not available' % objectType)

        logging.info('Seeding elements of object type %i: %s (%s)' %
                     (objectType, self.leewayprop[objectType]['OBJKEY'],
                      self.leewayprop[objectType]['Description']))

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
        for i in range(number):
            rdw[i] = np.random.randn(1)
            epsdw[i] = rdw[i]*dwstd
            # Avoid negative downwind slopes
            while downwindSlope[i] + epsdw[i]/20.0 < 0.0:
                rdw[i] = np.random.randn(1)
                epsdw[i] = rdw[i]*dwstd
        downwindEps = epsdw
        # NB
        # downwindEps = np.zeros(number)

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
        # crosswindEps = np.zeros(number)

        # Jibe probability
        if jibeProbability is None:
            jibeProbability = self.get_config('seed:jibeProbability')

        # Store seed data for ASCII format output
        self.ascii = {
             'lon': lon, 'lat': lat, 'radius': radius,
             'number': number, 'time': time}

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

    def list_object_categories(self, substr=None):
        '''Display leeway categories to screen

        Print only objects containing 'substr', if specified'''

        for i, p in enumerate(self.leewayprop):
            description = self.leewayprop[p]['Description']
            objkey = self.leewayprop[p]['OBJKEY']
            if substr is not None:
                if substr.lower() not in description.lower() + objkey.lower():
                    continue
            print('%i %s %s' % (i+1, objkey, description))

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

        # Move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Jibe elements randomly according to given probability
        jp_per_timestep = self.elements.jibeProbability * \
            np.abs(self.time_step.total_seconds()) / 3600.0
        jib = jp_per_timestep > np.random.random(self.num_elements_active())
        self.elements.crosswindSlope[jib] = - self.elements.crosswindSlope[jib]
        self.elements.orientation[jib] = 1 - self.elements.orientation[jib]
        logging.debug('Jibing %i out of %i elements.' %
                      (sum(jib), self.num_elements_active()))

    def export_ascii(self, filename):
        '''Export output to ASCII format of original version'''

        try:
            f = open(filename, 'w')
        except:
            raise ValueError('Could not open file for writing: '
                             + filename)

        start_time = self.start_time
        for inp in ['lon', 'lat', 'radius', 'time']:
            if len(np.atleast_1d(self.ascii[inp])) == 1:
                self.ascii[inp] = [self.ascii[inp], self.ascii[inp]]
        f.write('# Drift simulation initiated [UTC]:\n')
        f.write('simDate simTime\n')
        f.write(self.start_time.strftime('%Y-%m-%d\t%H:%M:%S\t#\n'))
        f.write(
            '# Model version:\n'
            'modelVersion\n'
            ' 3.00\n'  # OpenDrift version
            '# Object class id & name:\n'
            'objectClassId  objectClassName\n')
        try:
            objtype = self.elements.objectType[0]
        except:
            objtype = self.elements_deactivated.objectType[0]
        f.write(' %i\t%s\n' % (objtype,
                self.leewayprop[objtype]['OBJKEY']))
        f.write(
            '# Seeding start time, position & radius:\n'
            'startDate\tstartTime\tstartLon\tstartLat\tstartRad\n')
        f.write('%s\t%s\t%s\t%s\n' % (
            self.ascii['time'][0], self.ascii['lon'][0],
            self.ascii['lat'][0], self.ascii['radius'][0]/1000.))
        f.write(
            '# Seeding end time, position & radius:\n'
            'endDate\tendTime\tendLon\tendLat\tendRad\n')
        f.write('%s\t%s\t%s\t%s\n' % (
            self.ascii['time'][1], self.ascii['lon'][1],
            self.ascii['lat'][1], self.ascii['radius'][1]/1000.))
        seedDuration = (self.ascii['time'][1]-self.ascii['time'][0]).total_seconds()/60.
        seedSteps=seedDuration/(self.time_step_output.total_seconds()/60.)
        seedSteps = np.maximum(1, seedSteps)
        f.write( 
            '# Duration of seeding [min] & [timesteps]:\n'
            'seedDuration   seedSteps\n'
            '    %i      %i\n'
            '# Length of timestep [min]:\n'
            'timeStep\n' % (seedDuration, seedSteps))
        f.write('%i\n' % (self.time_step_output.total_seconds()/60.))
        f.write(
            '# Length of model simulation [min] & [timesteps]:\n'
            'simLength  simSteps\n')
        f.write('%i\t%i\n' % ((self.time - self.start_time).total_seconds()/60., self.steps_output))

        f.write(
            '# Total no of seeded particles:\n'
            'seedTotal\n'
            ' %s\n' % (self.num_elements_activated()))
        f.write(
            '# Particles seeded per timestep:\n'
            'seedRate\n'
            ' %i\n' % (self.num_elements_activated()/seedSteps))

        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()

        beforeseeded = 0
        lons, statuss = self.get_property('lon')
        lats, statuss = self.get_property('lat')
        orientations, statuss = self.get_property('orientation')
        for step in range(self.steps_output):
            lon = lons[step,:]
            lat = lats[step,:]
            orientation = orientations[step,:]
            status = statuss[step,:]
            num_active = np.sum(~status.mask)
            status[status.mask] = 41  # seeded on land
            lon[status.mask] = 0
            lat[status.mask] = 0
            orientation[status.mask] = 0
            status[status==0] = 11 # active
            status[status==1] = 41 # stranded
            ID = np.arange(0, num_active+1)
            f.write('\n# Date [UTC]:\nnowDate   nowTime\n')
            f.write((self.start_time + self.time_step_output*step).strftime('%Y-%m-%d\t%H:%M:%S\n'))
            f.write('# Time passed [min] & [timesteps], now seeded, seeded so far:\ntimePassed  nStep   nowSeeded   nSeeded\n')
            f.write('  %i\t%i\t%i\t%i\n' % (
                (self.time_step_output*step).total_seconds()/60,
                 step+1, num_active-beforeseeded, num_active))
            beforeseeded = num_active
            f.write('# Mean position:\nmeanLon meanLat\n')
            f.write('%f\t%f\n' % (np.mean(lon[status==11]),
                                  np.mean(lat[status==11])))
            f.write('# Particle data:\n')
            f.write('id  lon lat state   age orientation\n')
            age_minutes = self.time_step_output.total_seconds()*(
                           step - index_of_first)/60
            age_minutes[age_minutes<0] = 0
            for i in range(num_active):
                f.write('%i\t%s\t%s\t%i\t%i\t%i\n' % (i+1,
                        lon[i], lat[i], status[i], age_minutes[i],
                        orientation[i]))

        f.close()
