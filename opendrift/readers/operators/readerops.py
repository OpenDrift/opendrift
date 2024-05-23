from types import LambdaType
from ..basereader import BaseReader

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
        you want in your op. '''
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

        super().__init__()

    def covers_positions(self, lon, lat):
        return np.intersect1d(self.a.covers_positions(lon, lat), self.b.covers_positions(lon, lat))

    def covers_time(self, time):
        return self.a.covers_time(time) and self.b.covers_time(time)

    def get_variables_interpolated(self, variables, *args, **kwargs):
        assert set(variables).issubset(self.variables), f"{variables} is not subset of {self.variables}"

        debug=True
        if debug:
            print("args passed are: {}".format(args))
            print("kwargs passed are: {}".format(kwargs))

        env_a, env_profiles_a = self.a.get_variables_interpolated(variables, *args, **kwargs)
        env_b, env_profiles_b = self.b.get_variables_interpolated(variables, *args, **kwargs)

        variables = [
            var for var in env_a.keys() if var not in ['x', 'y', 'z', 'time']
        ]

	#Making disctinction between easy functions or more complex ones that need some external parameters
        for var in variables:
            if self.op_type == "easy":
                env_a[var] = self.op(env_a[var], env_b[var])
            elif self.op_type == "combine_gaussian":
                lon = kwargs['lon']
                lat = kwargs['lat']
                lon_center, lat_center = self.b.lon, self.b.lat
                std = self.external_params
                env_a[var] = self.op(env_a[var], env_b[var], lon, lat, lon_center, lat_center, std)
            else:
                raise ValueError('Op_type not recognised. You should verify the definition of the Reader operator you are using.')

        if env_profiles_a is not None:
            variables = [
                var for var in env_profiles_a.keys() if var not in ['x', 'y', 'z', 'time']
            ]
            for var in variables:
                env_profiles_a[var] = self.op(env_profiles_a[var], env_profiles_b[var])

        return env_a, env_profiles_a

