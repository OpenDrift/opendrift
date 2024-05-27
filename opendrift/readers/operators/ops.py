from abc import abstractmethod
from numbers import Number
from typing import List
import pyproj
import xarray as xr
import numpy as np
class Combine:
    """Combining two readers into a third one. You can use usual operators,
    but also more complex ones such as gaussian combining.
    """

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

    def combine_gaussian(self, measurement_reader, std):
        """Mix two readers with a gaussian, whose std is the one given as an argument.
        The measurment reader have to be of type timeseries, with a lon and lat
        attributes that are taken as the center of the measure.
        """
        from .readerops import Combined as ReaderCombined

        def gaussian_melting(a, b, lon, lat, lon_center, lat_center, std):
            geod = pyproj.Geod(ellps='WGS84')
            if  ( isinstance(lon_center, float) or len(lon_center) == 1 )  and len(lon) != 1 :
                lon_center = lon_center * np.ones(len(lon))
            if ( isinstance(lat_center, float) or len(lat_center) == 1 ) and len(lat) != 1 :
                lat_center = lat_center * np.ones(len(lat))
            _, _, distance = geod.inv(lon, lat, lon_center, lat_center)
            exponential_factor = np.exp( -np.power(distance/std, 2.) / 2)
            res = a * (1 - exponential_factor) + b * exponential_factor
            return res

        return ReaderCombined(a = self, b = measurement_reader, op = gaussian_melting, op_type= "combine_gaussian", external_params = std)



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

