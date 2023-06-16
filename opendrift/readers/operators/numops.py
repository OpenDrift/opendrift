"""
Reader combined with number.
"""

from numbers import Number
from types import LambdaType
from ..basereader import BaseReader

class Combined:
    """
    A reader combined with a number.
    """

    n: Number
    r: BaseReader
    op: LambdaType

    def __init__(self, n, r, op, descop = '|'):
        self.n = n
        self.r = r
        self.op = op

        assert isinstance(n, Number)
        assert isinstance(r, BaseReader)

        self.name = f'NumCombined({self.r} {descop} {self.n})'

    @staticmethod
    def add(n, r):
        return Combined(n, r, lambda x: n + x, '+')

    @staticmethod
    def mul(n, r):
        return Combined(n, r, lambda x: n * x, '*')

    @staticmethod
    def sub(n, r):
        return Combined(n, r, lambda x: x - n, '-')

    @staticmethod
    def div(n, r):
        return Combined(n, r, lambda x: x / n, '/')


    def __getattr__(self, attr):
        """
        Forward all other method calls and attributes to reader.
        """
        return getattr(self.r, attr)

    def get_variables_interpolated(self, variables, *args, **kwargs):
        assert set(variables).issubset(self.r.variables), f"{vars} is not a subset of variables in {self.r}"

        env, env_profiles = self.r.get_variables_interpolated(variables, *args, **kwargs)

        variables = [
            var for var in env.keys() if var not in ['x', 'y', 'z', 'time']
        ]

        for var in variables:
            env[var] = self.op(env[var])

        if env_profiles is not None:
            variables = [
                var for var in env_profiles.keys() if var not in ['x', 'y', 'z', 'time']
            ]
            for var in variables:
                env_profiles[var] = self.op(env_profiles[var])

        return env, env_profiles

