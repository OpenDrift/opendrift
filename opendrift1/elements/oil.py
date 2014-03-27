import numpy as np

from elements import LagrangianArray


class Oil(LagrangianArray):

    parameters = LagrangianArray.add_parameters([
        ('massOil', {'dtype': np.float32,
                     'unit': 'kg'}),
        ('massEvaporated', {'dtype': np.float32,
                            'unit': 'kg',
                            'default': 0}),
        ('massEmulsion', {'dtype': np.float32,
                          'unit': 'kg',
                          'default': 0})])

    def update_properties(self, environment):
        # evaporate 10% of oil
        evaporated_mass = self.massOil*0.1*environment.windspeed
        self.massOil = self.massOil - evaporated_mass
        self.massEvaporated = self.massEvaporated + evaporated_mass

     def update_position(self, environment):
        self.lon = self.lon + envoronment.wind_east*2
