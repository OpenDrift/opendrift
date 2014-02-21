from lagrangian_element import LagrangianElement

class OilElement(LagrangianElement):

    variables = LagrangianElement.variables + [
                    ('massOil', 'float32'),
                    ('massEvaporated', 'float32'),
                    ('massEmulsion', 'float32')]
