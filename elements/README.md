opendrift - Open source framework for ocean trajectory modelling

elements - LagrangianArray
==========================


The module 'elements.py' defines the fundamental class 'LagrangianArray'

The LagrangianArray class holds a dictionary (OrderedDict) defining the variables (properties) of an element.
The variables are defined by:
```
- name
- dtype (typically numpy.float32)
```
and optionally:
```
- unit (following udunits convention)
- standard_name (to support output files compliant with CF convention)
- default (default value, e.g. 0 for depth)
```


The values for these variables are given through the constructor, and are stored as named attributes of the created object:

```
>>> from elements import LagrangianArray

>>> l = LagrangianArray(lon=[5, 6, 7], lat=[60, 60, 60], depth=10)

>>> print l.lat
[ 60.  60.  60.]
```
The input values are arrays, where a single element is simply a special case of LagrangianArray with length 1.

The element properties are updated by the method 'update_properties', which does nothing for a general LagrangianArray


```
l.update_properties()
```

New element types can be created by subclassing the LagrangianArray class and:
- overloading the 'update_properties' class 
- and (possibly) adding more variables through method 'add_variables'

Example for Larvae:

```
class Larvae(LagrangianArray):

    variables = LagrangianArray.add_variables(
        {'length':
            {'dtype': np.float32}})

    def update_properties(self):
        self.length = self.length*1.01  # General larvae grow by 1%

```
New modules may easily be created from a template (template_elements.py).
Multiple inheritance is supported:

```
class CodLarvae(Larvae):

    variables = Larvae.add_variables(
        {'CodLarvaeProperty1':
            {'dtype': np.float32}})


class HalibutLarvae(Larvae):

    variables = Larvae.add_variables(
        {'HalibutLarvaeProperty1':
            {'dtype': np.float32}})

    def update_properties(self):
        self.length = self.length*1.02  # Halibut larvae grow by 2%

```


Example CodLarvae:

```
l = CodLarvae(lon=[5, 6, 7], lat=[60, 60, 60], CodLarvaeProperty1=[1], length=10)
l.update_properties()
l.length
10.1
```

Example HalibutLarvae:

```
l = HalibutLarvae(lon=[5, 6, 7], lat=[60, 60, 60], HalibutLarvaeProperty1=[1], length=10)
l.update_properties()
l.length
10.2
```
