from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Optional
import numpy as np


class Spectrum(Enum):
    JONSWAP = 1
    Phillips = 2



class T(ABC):
    """
    Generic Wave Period.
    """
    _t_: float
    js: 'SpecCoeff'
    phillips: 'SpecCoeff'

    @property
    def t(self):
        return self._t_

    @property
    def T(self):
        return self._t_

    def __init__(self, t):
        self._t_ = t

    def tp(self, spectrum: Spectrum) -> 'T':
        """
        Return a interpretation of the wave period in the given spectrum.
        """
        return self.into(Tp, spectrum)

    def into(self, ty: 'T', spectrum) -> 'T':
        if isinstance(self, ty):
            return self

        match spectrum:
            case Spectrum.JONSWAP:
                return ty(self.js.coeff(ty) * self._t_)



@dataclass
class SpecCoeff:
    tp: Optional[float]
    tm_01: float
    tm_10: float
    tm_02: float
    tm_03: float

    def coeff(self, ty) -> float:
        if ty == Tp:
            return self.tp
        elif ty == Tm01:
            return self.tm_01
        elif ty == Tm_10:
            return self.tm_10
        elif ty == Tm02:
            return self.tm_02

class JONSWAP(SpecCoeff): pass
class Phillips(SpecCoeff): pass

class Tp(T):
    js = SpecCoeff(tp=None,
                      tm_01=0.098,
                      tm_10=0.093,
                      tm_02=0.096,
                      tm_03=0.093)


class Tm01(T):
    js = SpecCoeff(tp=1.,
                      tm_01=None,
                      tm_10=0.093,
                      tm_02=0.096,
                      tm_03=0.093)


class Tm_10(T):
    js = SpecCoeff(tp=1.,
                      tm_01=0.098,
                      tm_10=0.093,
                      tm_02=0.096,
                      tm_03=0.093)


class Tm02(T):
    js = SpecCoeff(tp=1.,
                      tm_01=0.098,
                      tm_10=0.093,
                      tm_02=0.096,
                      tm_03=0.093)


class Tm03(T):
    js = SpecCoeff(tp=1.,
                      tm_01=0.098,
                      tm_10=0.093,
                      tm_02=0.096,
                      tm_03=0.093)

