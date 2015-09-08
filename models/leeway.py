import os
import numpy as np
from datetime import datetime

from opendrift import OpenDriftSimulation
from elements import LagrangianArray
from collections import OrderedDict

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
        ('downwindLeeway', {'dtype': np.float32,
                            'unit': '%',
                            'default': 1}),
        ('crosswindLeeway', {'dtype': np.float32,
                            'unit': '%',
                            'default': 1}),
        ('downwindOffset', {'dtype': np.float32,
                          'unit': 'cm/s',
                          'default': 0}),
        ('crosswindOffset', {'dtype': np.float32,
                          'unit': 'cm/s',
                          'default': 0}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
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
        objproptxt=objprop_file.readlines()
        objprop_file.close()
        MAXOBJ=85 # CCC must fix
        self.leewayprop=OrderedDict({})
        for i in xrange(MAXOBJ):
            line = objproptxt[i*3].split()[0]
            key = objproptxt[i*3].split()[0]
            props = {'OBJNUMBER':i}
            props['Description'] = objproptxt[i*3+1]
            arr = [float(x) for x in objproptxt[i*3+2].split()]
            proplist = {}
            proplist['DWSLOPE']  = arr[0]
            proplist['DWOFFSET'] = arr[1]
            proplist['DWSTD']    = arr[2]
            proplist['CWRSLOPE'] = arr[3]
            proplist['CWROFFSET']= arr[4]
            proplist['CWRSTD']   = arr[5]
            proplist['CWLSLOPE'] = arr[6]
            proplist['CWLOFFSET']= arr[7]
            proplist['CWLSTD']   = arr[8]
            self.leewayprop[key] = proplist

        #objno = 0 # CCC default object no
        #objkey = "PIW-1" # CCC default object key

        # Calling general constructor of parent class
        super(Leeway, self).__init__(*args, **kwargs)

    def seed_point_leeway(self, lon, lat, radius, number, time=None, objectType=0):
        """Seed a given number of particles spread around a given position.
        """
        #orientation = np.int_(random.random(number)+0.5)
        # Left is -1, Right of downwind is 1
        orientation = 2*(np.random.rand(number)>0.5)-1
        cwslo=self.leewayprop.values()[objectType]['CWRSLOPE']
        cwoff=self.leewayprop.values()[objectType]['CWRSTD']
        cwstd=self.leewayprop.values()[objectType]['CWROFFSET']
        crosswindLeeway = 0.002*orientation
        downwindLeeway = 0.02

        self.seed_point(lon, lat, radius, number, time,
             orientation=orientation, objectType=objectType,
             downwindLeeway=downwindLeeway, crosswindLeeway=crosswindLeeway)
        #jibeProbability=jibeProbability, 


    def update(self):
        """Update positions and properties of leeway particles."""

        self.elements.age_seconds += self.time_step.total_seconds()

        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)
        # CCC update wind direction
        winddir = np.arctan2(self.environment.x_wind, self.environment.y_wind)
        
        # Move particles with the leeway CCC TODO
        dwl_leeway = self.elements.downwindLeeway*windspeed + \
                     self.elements.downwindOffset
        cwl_leeway = self.elements.orientation * \
                     (self.elements.crosswindLeeway*windspeed + \
                      self.elements.crosswindOffset)
        x_leeway = dwl_leeway*np.cos(winddir)+cwl_leeway*np.sin(winddir)
        y_leeway = -dwl_leeway*np.sin(winddir)+cwl_leeway*np.cos(winddir)
        self.update_positions(x_leeway, y_leeway)

        # Deactivate elements on land
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

        # Move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)
