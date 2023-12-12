# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2021, Gaute Hope, MET Norway
"""
Interface to the ADIOS oil database.
"""

import logging

logger = logging.getLogger(__name__)

import numpy as np
from typing import List

from adios_db.models.oil.oil import Oil as AdiosOil
from adios_db.computation import gnome_oil
from adios_db.computation import physical_properties
from adios_db.computation.estimations import oil_water_surface_tension_from_api


class NotFullOil(Exception):
    pass


def __require_gnome_oil__(f):
    def w(self, *args, **kwargs):

        if self.gnome_oil is None:
            raise NotFullOil()

        return f(self, *args, **kwargs)

    return w


class OpendriftOil:
    id: str
    type: str
    name: str
    API: float
    gnome_suitable: bool
    labels: List[str]
    location: str
    model_completeness: float
    product_type: str
    sample_date: str
    data: dict
    oil: AdiosOil
    gnome_oil: dict

    def __init__(self, o):
        self.data = o

        data = o
        meta = o['metadata']
        self.id = o['oil_id']
        self.name = o['metadata']['name']

        #data = o['data']
        #meta = data['attributes']['metadata']
        #self.id = data['_id']
        #self.name = meta['name']

        logger.debug(f'Parsing Oil: {self.id} / {self.name}')
        #self.oil = AdiosOil.from_py_json(data['attributes'])
        self.oil = AdiosOil.from_py_json(data)

        if not self.oil.metadata.gnome_suitable:
            logger.error(f'{self.id} / {self.name}: is not GNOME suitable')
        else:
            self.gnome_oil = gnome_oil.make_gnome_oil(self.oil)

    def __repr__(self):
        return f"[<adios.Oil> {self.id}] {self.name}"

    def json(self):
        import json
        return json.dumps(self.data)

    def make_full(self) -> 'OpendriftOil':
        return self

    def valid(self):
        """
        Check whether this oil can be used in Opendrift simulations.
        """
        return self.oil.metadata.gnome_suitable

    @__require_gnome_oil__
    def density_at_temp(self, t, unit='K') -> float:
        """
        Return density at temperature (in Kelvin by default).
        """
        return physical_properties.Density(self.oil).at_temp(t, unit)

    @__require_gnome_oil__
    def kvis_at_temp(self, t, unit='K') -> float:
        return physical_properties.KinematicViscosity(self.oil).at_temp(
            t, temp_units=unit)

    @property
    @__require_gnome_oil__
    def mass_fraction(self) -> np.ndarray:
        return np.asarray(self.gnome_oil['mass_fraction'])

    @__require_gnome_oil__
    def oil_water_surface_tension(self) -> float:
        return oil_water_surface_tension_from_api(self.gnome_oil['api'])

    @property
    @__require_gnome_oil__
    def bulltime(self) -> float:
        bulltime = self.gnome_oil['bullwinkle_time']
        if bulltime is None:
            return -999.  # legacy from old oil_library
        else:
            return bulltime

    @property
    @__require_gnome_oil__
    def bullwinkle(self) -> float:
        return self.gnome_oil['bullwinkle_fraction']

    @property
    @__require_gnome_oil__
    def emulsion_water_fraction_max(self) -> float:
        return self.gnome_oil['emulsion_water_fraction_max']

    @__require_gnome_oil__
    def vapor_pressure(self, temp) -> float:
        '''
        Calculate vapor pressure. This method is taken from the old oil_library.

        Args:

            temp: temperature in Kelvin.

        Returns:

            Array of vapor pressures for each component. Pascal.
        '''
        atmos_pressure = 101325.0
        boiling_point = np.asarray(self.gnome_oil['boiling_point'])

        D_Zb = 0.97
        R_cal = 1.987  # calories

        D_S = 8.75 + 1.987 * np.log(boiling_point)
        C_2i = 0.19 * boiling_point - 18

        var = 1. / (boiling_point - C_2i) - 1. / (temp - C_2i)
        ln_Pi_Po = (D_S * (boiling_point - C_2i)**2 /
                    (D_Zb * R_cal * boiling_point) * var)
        Pi = np.exp(ln_Pi_Po) * atmos_pressure

        return Pi

    @property
    @__require_gnome_oil__
    def molecular_weight(self):
        return np.asarray(self.gnome_oil['molecular_weight'])

    @property
    @__require_gnome_oil__
    def k0y(self) -> float:
        return 2.024e-06  # Constant for all oils, should be checked
        #return self.gnome_oil['k0y']
