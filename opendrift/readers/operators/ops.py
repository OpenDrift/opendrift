from numbers import Number
import numpy as np

class Combine:
    def __add__(self, other):
        from .reader import Add
        from .numops import Combined as NumCombined
        from ..basereader import BaseReader

        if isinstance(other, Number):
            return NumCombined.add(other, self)
        elif isinstance(other, BaseReader):
            return Add(self, other)
        else:
            return NotImplemented

    def __mul__(self, other):
        from .numops import Combined as NumCombined
        if isinstance(other, Number):
            return NumCombined.mul(other, self)
        else:
            return NotImplemented

    def __div__(self, other):
        from .numops import Combined as NumCombined

        if isinstance(other, Number):
            return NumCombined.div(other, self)
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (-1 * other)

