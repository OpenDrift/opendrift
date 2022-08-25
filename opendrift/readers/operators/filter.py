from .ops import Combine, Filter

class FilterVariables(Combine, Filter):
    """
    A reader where only some variables are forwarded.
    """

    v = None

    @property
    def variables(self):
        return self.v

    def __init__(self, r, vars):
        self.r = r

        assert set(vars).issubset(self.r.variables), f"{vars} is not a subset of variables in {self.r}"

        self.v = vars

        self.name = f'Filter({self.r} | {self.v})'

    def __getattr__(self, attr):
        """
        Forward all other method calls and attributes to reader.
        """
        return getattr(self.r, attr)
