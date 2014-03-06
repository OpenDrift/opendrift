import numpy as np

from elements import LagrangianArray

class Oil(LagrangianArray):

    parameters = LagrangianArray.parameters.update({
        'massOil': {'dtype': np.float32,
                    'unit': 'kg'},
        'massEvaporated': {'dtype': np.float32,
                            'default': 0},
        'massEmulsion': {'dtype': np.float32,
                          'default': 0}})


    def update_properties(self):
        # evaporate 10% of oil
        evaporated_mass = self.massOil*0.1
        self.massOil = self.massOil - evaporated_mass
        self.massEvaporated = self.massEvaporated + evaporated_mass
