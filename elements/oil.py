import numpy as np

from elements import LagrangianArray


class Oil(LagrangianArray):

    variables = LagrangianArray.add_variables([
        ('massOil', {'dtype': np.float32,
                     'unit': 'kg'}),
        ('massEvaporated', {'dtype': np.float32,
                            'unit': 'kg',
                            'default': 0}),
        ('massEmulsion', {'dtype': np.float32,
                          'unit': 'kg',
                          'default': 0})])
