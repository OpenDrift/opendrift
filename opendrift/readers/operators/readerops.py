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

    def __init__(self, a, b, op):
        self.a = a
        self.b = b
        self.op = op

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

        env_a, env_profiles_a = self.a.get_variables_interpolated(variables, *args, **kwargs)
        env_b, env_profiles_b = self.b.get_variables_interpolated(variables, *args, **kwargs)

        variables = [
            var for var in env_a.keys() if var not in ['x', 'y', 'z', 'time']
        ]

        for var in variables:
            env_a[var] = self.op(env_a[var], env_b[var])

        if env_profiles_a is not None:
            variables = [
                var for var in env_profiles_a.keys() if var not in ['x', 'y', 'z', 'time']
            ]
            for var in variables:
                env_profiles_a[var] = self.op(env_profiles_a[var], env_profiles_b[var])

        return env_a, env_profiles_a

