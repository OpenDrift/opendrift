from abc import ABC, abstractmethod, abstractproperty
from enum import Enum

class Spectrum(Enum):
    JONSWAP = 1

class T(ABC):
    """
    Generic Wave Period.
    """

    @abstractmethod
    def tp(self, spectrum: Spectrum) -> float:
        """
        Return a interpretation of the wave period in the given spectrum.
        """
        raise NotImplemented()


class Tp(T):
    _tp_: float

    def __init__(self, tp):
        self._tp_ = tp

    def tp(self, spectrum: Spectrum):
        match spectrum:
            case Spectrum.JONSWAP:
                return self._tp_

        raise ValueError()

class Tm0(T):
    _tm0_: float

    def __init__(self, tm0):
        self._tm0_ = tm0

    def tp(self, spectrum: Spectrum):
        match spectrum:
            case Spectrum.JONSWAP:
                return self._tm0_ * 1.1

        raise ValueError()
