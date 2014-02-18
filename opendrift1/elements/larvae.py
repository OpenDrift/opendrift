from lagrangian_element import LagrangianElement

class LarvaeElement(LagrangianElement):

    variables = LagrangianElement.variables + ['length']

    @classmethod
    def update_properties(cls):
        super(LarvaeElement, cls).update_properties()
        print '...using specialised function for Larvae class'
