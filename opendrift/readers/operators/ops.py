from abc import abstractmethod
from numbers import Number
from typing import List
import pyproj
import xarray as xr
import numpy as np
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

    def combine_gaussian(self, measurement_reader, std):
        from .readerops import Combined as ReaderCombined

        def gaussian_melting(a, b, lon, lat, lon_center, lat_center, std):
            geod = pyproj.Geod(ellps='WGS84')
            if  ( isinstance(lon_center, float) or len(lon_center) == 1 )  and len(lon) != 1 :
                lon_center = lon_center * np.ones(len(lon))
            if ( isinstance(lat_center, float) or len(lat_center) == 1 ) and len(lat) != 1 :
                lat_center = lat_center * np.ones(len(lat))

            print("lon length: {}".format(len(lon)))
            print("lat length: {}".format(len(lat)))
            print("lon_center length: {}".format(len(lon_center)))
            print("lat_center length: {}".format(len(lat_center)))

            _, _, distance = geod.inv(lon, lat, lon_center, lat_center)
            exponential_factor = np.exp( -np.power(distance/std, 2.) / 2)
            res = a * (1 - exponential_factor) + b * exponential_factor
            print("==================================")
            print("Gaussian calculation")
            print("Lon, lat, lon_center, lat_center = {}".format((lon, lat, lon_center, lat_center)))
            print("Distance used for calculation: {}".format(distance))
            print("a = {}, b = {}, c = {}".format(a, b, res))
            print("=================================")
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

