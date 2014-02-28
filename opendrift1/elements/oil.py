from elements import LagrangianArray

class Oil(LagrangianArray):

    variables = LagrangianArray.variables + [
                        # name, data type[, default value]
                        ('massOil', 'float32'),
                        ('massEvaporated', 'float32', 0),
                        ('massEmulsion', 'float32', 0)]
