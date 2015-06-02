import numpy as np

from elements import LagrangianArray


class Oil(LagrangianArray):
    """Extending LagrangianArray with variables relevant for oil particles."""

    variables = LagrangianArray.add_variables([
        ('massOil', {'dtype': np.float32,
                     'unit': 'kg'}),
        ('viscosity', {'dtype': np.float32,
                            'unit': 'mm2/s (centiStokes)',
                            'default': 100}),
        ('density', {'dtype': np.float32,
                            'unit': 'kg/cm^3',
                            'default': .8}),
        ('massEvaporated', {'dtype': np.float32,
                            'unit': 'kg',
                            'default': 0}),
        ('massEmulsion', {'dtype': np.float32,
                          'unit': 'kg',
                          'default': 0})])
