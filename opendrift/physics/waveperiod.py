from abc import ABC, abstractmethod, abstractproperty

class T(ABC):
    """
    Generic Wave Period.
    """

    @abstractproperty
    def jonswap(self):
        """
        Return a JONSWAP interpretation of the wave period.
        """
        pass

class JONSWAP:
    """
    """
    _tp_: float

    def __init__(self, tp):
        self._tp_ = tp

    @classmethod
    def from_tp(cls, tp):
        cls(tp)

    @property
    def tp(self):
        return self.tp

class Tp(T):
    _tp_: float

    def __init__(self, tp):
        self._tp_ = tp

    @property
    def jonswap(self):
        return JONSWAP.from_tp(self._tp_)

class Tm0(T):
    _tm0_: float

    def __init__(self, tm0):
        self._tm0_ = tm0

    @property
    def jonswap(self):
        return JONSWAP.from_tp(1.1 * self._tm0_)
