import numpy as np

from elements import LagrangianArray


class PassiveTracer(LagrangianArray):
    """Basic implementation of LagrangianArray with no additional properties.

    Contains only the properties of the abstract class LagrangianArray,
    i.e. position (lon, lat, depth) and ID.
    May be used for passive tracer calculations when no properties are needed.
    """

    pass
