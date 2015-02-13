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


New element types can be created by subclassing the LagrangianArray class and adding more variables through method 'add_variables'

Example for Larvae:

```
class Larvae(LagrangianArray):

    variables = LagrangianArray.add_variables(
        {'length':
            {'dtype': np.float32,
             'unit': 'mm'}})

```
New modules may easily be created from a template (template_elements.py).
Multiple inheritance is supported:

```
class CodLarvae(Larvae):
    """Extending Larvae with variables relevant for cod larvae."""

    variables = Larvae.add_variables(
        {'CodLarvaeProperty1':
            {'dtype': np.float32}})


class HalibutLarvae(Larvae):
    """Extending Larvae with variables relevant for halibut larvae."""

    variables = Larvae.add_variables(
        {'HalibutLarvaeProperty1':
            {'dtype': np.float32}})
```


Example CodLarvae:

```
l = CodLarvae(lon=[5, 6, 7], lat=[60, 60, 60], CodLarvaeProperty1=[1], length=10)
l.length
10
```
