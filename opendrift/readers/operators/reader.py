from ..basereader import BaseReader

class Combined(BaseReader):
    """
    A combination of two readers.
    """
    a: 'BaseReader'
    b: 'BaseReader'

    def __init__(self, a, b, variables = None):
        self.a = a
        self.b = b

        if variables is not None:
            assert set(self.a.variables).intersection(self.b.variables).issuperset(variables), f"{variables} is not a subset of variables available in both readers {a} and {b}"
            self.variables = variables
        else:
            self.variables = list(set(self.a.variables).intersection(self.b.variables))

        self.start_time = max(self.a.start_time, self.b.start_time)
        self.end_time = min(self.a.end_time, self.b.end_time)

        self.name = f'Combined({a.name} | {b.name})'

        super()

    def covers_positions(self, lon, lat):
        return np.intersect1d(self.a.covers_positions(lon, lat), self.b.covers_positions(lon, lat))

    def covers_time(self, time):
        return self.a.covers_time(time) and self.b.covers_time(time)


class Add(Combined):
    def __init__(self, a, b):
        super(a, b)
