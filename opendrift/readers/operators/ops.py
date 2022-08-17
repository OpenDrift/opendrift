from abc import abstractmethod
from numbers import Number
from typing import List

class Combine:
    def __add__(self, other):
        from .readerops import Combined as ReaderCombined
        from .numops import Combined as NumCombined
        from ..basereader import BaseReader

        if isinstance(other, Number):
            return NumCombined.add(other, self)
        elif isinstance(other, BaseReader):
            return ReaderCombined(self, other, lambda a, b: a + b)
        else:
            return NotImplemented

    def __mul__(self, other):
        from .numops import Combined as NumCombined
        if isinstance(other, Number):
            return NumCombined.mul(other, self)
        else:
            return NotImplemented

    def __truediv__(self, other):
        from .numops import Combined as NumCombined

        if isinstance(other, Number):
            return NumCombined.div(other, self)
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (-1 * other)

class Filter:
    @property
    @abstractmethod
    def variables(self) -> List[str]:
        pass

    def filter_vars(self, vars):
        """
        Only keep the specified variables.
        """
        from .filter import FilterVariables
        return FilterVariables(self, vars)

    def exclude_vars(self, vars):
        """
        Remove the specified variables.
        """
        from .filter import FilterVariables

        vars = list(set(self.variables) - set(vars))
        return FilterVariables(self, vars)

