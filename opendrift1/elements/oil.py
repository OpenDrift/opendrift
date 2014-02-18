from lagrangian_element import LagrangianElement

class OilElement(LagrangianElement):

    variables = LagrangianElement.variables + ['massOil', 'massEvaporated', 'massEmulsion']
