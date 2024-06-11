from types import LambdaType
from ..basereader import BaseReader
import numpy as np
import matplotlib.pyplot as plt
import pyproj

def none_or_cmp(a, b, cmp):
    if a is None:
        return b
    if b is None:
        return a

    return cmp(a, b)

class Combined(BaseReader):
    """
    A combination of two readers.
    """
    a: BaseReader
    b: BaseReader
    op: LambdaType

    def __init__(self, a, b, op, op_type = "easy", external_params = None):
        '''Combine two readers a and b followinf the operator op. If needed,
        you can ad an op_type that will enable you to use the external parameters
        you want in your op between lines 63 and 73.'''
        self.a = a
        self.b = b
        self.op = op
        self.op_type = op_type
        self.external_params = external_params

        self.variables = list(set(self.a.variables).intersection(self.b.variables))

        self.start_time = none_or_cmp(self.a.start_time, self.b.start_time, max)
        self.end_time = none_or_cmp(self.a.end_time, self.b.end_time, min)

        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90

        self.name = f'Combined({a.name} | {b.name})'

        self.proj4 = '+proj=latlong'
        self.projected = None
        super().__init__()

    def covers_positions(self, lon, lat):
        return np.intersect1d(self.a.covers_positions(lon, lat), self.b.covers_positions(lon, lat))

    def covers_time(self, time):
        return self.a.covers_time(time) and self.b.covers_time(time)

    def get_variables_interpolated(self, variables, profiles=None, profiles_depth=None,time=None,lon=None, lat=None, z=None,rotate_to_proj=None):
        assert set(variables).issubset(self.variables), f"{variables} is not subset of {self.variables}"

        env_a, env_profiles_a = self.a.get_variables_interpolated(variables, time=time,lon=lon, lat=lat, z=z)
        env_b, env_profiles_b = self.b.get_variables_interpolated(variables,time=time,lon=lon, lat=lat, z=z)
        variables = [ var for var in env_a.keys() if var not in ['x', 'y', 'z', 'time'] ]

        #Making disctinction between easy functions or more complex ones that need some external parameters
        env_c = {}
        for var in variables:
            if self.op_type == "easy":
                env_c[var] = self.op(env_a[var], env_b[var])
            elif self.op_type == "combine_gaussian":
                lon_center, lat_center = self.b.lon, self.b.lat
                std = self.external_params
                env_c[var] = self.op(env_a[var], env_b[var], lon, lat, lon_center, lat_center, std)
            else:
                raise ValueError('Op_type not recognised. You should verify the definition of the Reader operator you are using.')

        #Profile
        env_profiles_c = None
        if env_profiles_a is not None:
            env_profiles_c = {}
            variables = [ var for var in env_profiles_a.keys() if var not in ['x', 'y', 'z', 'time'] ]
            for var in variables:
                env_profiles_c[var] = self.op(env_profiles_c[var], env_profiles_c[var])

        return env_c, env_profiles_c

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):

        pm1 = self.a.pixel_size()
        pm2 = self.b.pixel_size()
        if pm1 == None and pm2 == None:
            raise Exception("Neither {} or {} in {} have pixel_size well defined. Plot is not possible.")
        elif pm1 == None:
            pm=pm2
        elif pm2 == None:
            pm=pm1
        else:
            pm = np.lcm(int(pm1), int(pm2))

        delta_x = pm/111000
        delta_y = pm/ 111000 * np.abs(np.cos(np.radians(np.mean(y))))
        x = np.arange(np.min(x), np.max(x), delta_x)
        y = np.arange(np.min(y), np.max(y), delta_y)
        X, Y = np.meshgrid(x, y)
        shape = X.shape
        X = X.flatten()
        Y = Y.flatten()

        variables, _ = self.get_variables_interpolated(requested_variables, lon = X, lat = Y, time = time, z = z)
        for key in variables.keys():
            variables[key] = np.reshape(variables[key], shape)
        variables['x'] = x
        variables['y'] = y

        return variables

