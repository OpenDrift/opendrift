"""
utilities for doing computation on the physical properties of an
oil record
"""

from operator import itemgetter
import numpy as np

import unit_conversion as uc


class Density:
    """
    class to hold and do calculations on density

    data is stored internally in standard units:
    temperature in Kelvin
    density in kg/m^3
    """
    def __init__(self, oil):
        """
        Initialize a density calculator

        :param oil: an Oil object -- the density data will be extracted

        or

        :param oil: Sequence of density/temperature pairs:
                    ``[(980.0, 288.15), (990.0, 273.15)])``

        If data pairs, units must be kg/m^3 and K
        """
        try:
            data = get_density_data(oil, units='kg/m^3', temp_units="K")
        except AttributeError:
            # not an oil object -- assume it's a table of data in the
            #                      correct form
            data = oil
        if data:
            data = sorted(data, key=itemgetter(1))
            self.densities, self.temps = zip(*data)
        else:
            self.densities = []
            self.temps = []
        self.initialize()

    def initialize(self):
        """
        Initialize the expansion coefficient

        For outside the measured range
        """
        # if there is only one density, use a default
        # Note: no idea where these values came from

        if not np.all(np.diff(self.temps) > 0):
            raise ValueError("temperatures must be discreet")

        if len(self.densities) == 1:
            d = self.densities[0]
            t = self.temps[0]
            if abs(t - 288.16) < 5.0:  # measurement within 5 deg of 15 C
                # API 30 threshold
                self.k_rho_default = -0.0009 if d < 875 else -0.0008
            else:
                self.k_rho_default = -0.00085  # who knows?
        elif len(self.densities) > 1:
            # do a linear fit to the points
            # this should exactly match if there are only two.
            b = self.densities
            A = np.c_[np.ones_like(b), np.array(self.temps)]
            x, _residuals, _rank, _s = np.linalg.lstsq(A, b, rcond=None)
            self.k_rho_default = x[1]
        else:
            raise ValueError("Density needs at least one density value")

    def at_temp(self, temp, unit='K'):
        """
        density(s) at the provided temperature(s)

        :param temp: scalar or sequence of temp in K

        :param unit='K': unit of temperature

        densities will be returned as kg/m^3
        """
        temp = np.asarray(temp)
        scaler = True if temp.shape == () else False
        temp.shape = (-1, )

        if unit != 'K':
            temp = uc.convert(unit, 'K', temp)

        densities = np.interp(temp, self.temps, self.densities,
                              left=-np.inf, right=np.inf)

        left = (densities == -np.inf)
        densities[left] = self.densities[0] + (self.k_rho_default *
                                               (temp[left] - self.temps[0]))

        right = (densities == np.inf)
        densities[right] = self.densities[-1] + (self.k_rho_default *
                                                 (temp[right] - self.temps[-1]))

        return densities if not scaler else densities[0]


class KinematicViscosity:
    """
    class to hold and do calculations on kinematic viscosity

    data is stored internally in standard units:
    temperature in Kelvin
    viscosity in m^2/s
    """
    def __init__(self, oil):
        """
        initialize from an oil object
        """
        try:
            data = get_kinematic_viscosity_data(oil, units='m^2/s',
                                                temp_units="K")
        except AttributeError:
            # not an oil object -- assume it's a table of data in the
            #                      correct form
            data = oil

        if data:
            data = sorted(data, key=itemgetter(1))
            self.kviscs, self.temps = zip(*data)
        else:
            self.kviscs = []
            self.temps = []
        self.initialize()

    def at_temp(self, temp, kvis_units='m^2/s', temp_units="K"):
        """
        Compute the kinematic viscosity of the oil as a function of temperature

        :param temp_k: temperatures to compute at: can be scalar or array
                       of values.  Should be in Kelvin

        viscosity as a function of temp is given by:
        v = A exp(k_v2 / T)

        with constants determined from measured data
        """
        temp = np.asarray(temp)
        temp = uc.convert('temperature', temp_units, 'K', temp)

        kvisc = self._visc_A * np.exp(self._k_v2 / temp)

        kvisc = uc.convert('kinematic viscosity', 'm^2/s', kvis_units, kvisc)
        return kvisc

    def initialize(self):
        '''
        viscosity as a function of temp is given by:

        v = A exp(k_v2 / T)

        The constants, A and k_v2 are determined from the viscosity data:

        If only one data point, a default value for k_vs is used:
           2100 K, based on analysis of data in the ADIOS database as of 2018

        If two data points, the two constants are directly computed

        If three or more, the constants are computed by a least squares fit.
        '''
        # find viscosity measurements with zero weathering

        # this sets:
        self._k_v2 = None  # decay constant for viscosity curve
        self._visc_A = None

        kvis = self.kviscs
        kvis_ref_temps = self.temps

        if len(kvis) == 1:  # use default k_v2
            self._k_v2 = 2100.0
            self._visc_A = kvis[0] * np.exp(-self._k_v2 / kvis_ref_temps[0])
        else:
            # do a least squares fit to the data
            # viscs = np.array(kvis)
            # temps = np.array(kvis_ref_temps)
            b = np.log(kvis)
            A = np.c_[np.ones_like(b), 1.0 / np.array(kvis_ref_temps)]
            x, _residuals, _rank, _s = np.linalg.lstsq(A, b, rcond=None)
            self._k_v2 = x[1]
            self._visc_A = np.exp(x[0])
        return


def get_density_data(oil, units="kg/m^3", temp_units="K"):
    """
    Return a table of density data:

    list of (density, temp) pairs

    :param oil: the oil object to get data from

    :param units="kg/m^3": units you want the density in

    :param temp_units="K": units you want the density in

    """

    densities = [d for d in oil.sub_samples[0].physical_properties.densities
                 if d.density is not None
                 and d.ref_temp is not None]

    # create normalized list of densities
    density_table = []
    for density_point in densities:
        d = density_point.density.converted_to(units).value
        t = density_point.ref_temp.converted_to(temp_units).value
        density_table.append((d, t))
    return density_table


def get_kinematic_viscosity_data(oil, units="m^2/s", temp_units="K"):
    """
    Return a table of kinematic viscosity data:

    list of (viscosity, temp) pairs

    :param oil: the oil object to get data from

    :param units="cSt": units you want the viscosity in

    :param temp_units="K": units you want the viscosity in

    """

    try:
        kvisc = [k for k
                 in oil.sub_samples[0].physical_properties.kinematic_viscosities
                 if k.viscosity is not None
                 and k.ref_temp is not None]
    except IndexError:  # no subsamples at all!
        return []

    if len(kvisc) > 0:
        visc_table = []
        for visc_point in kvisc:
            d = visc_point.viscosity.converted_to(units).value
            t = visc_point.ref_temp.converted_to(temp_units).value
            visc_table.append((d, t))
        return visc_table

    dvisc = oil.sub_samples[0].physical_properties.dynamic_viscosities
    if len(dvisc) > 0:
        dvisc = get_dynamic_viscosity_data(oil, units="Pa s", temp_units="K")
        density = Density(oil)
        visc_table = convert_dvisc_to_kvisc(dvisc, density)
        return visc_table
    else:
        return []


def get_dynamic_viscosity_data(oil, units="Pas", temp_units="K"):
    """
    Return a table of kinematic viscosity data:

    list of (viscosity, temp) pairs

    :param oil: the oil object to get data from

    :param units="cSt": units you want the viscosity in

    :param temp_units="K": units you want the viscosity in

    """

    dvisc = [d for d
             in oil.sub_samples[0].physical_properties.dynamic_viscosities
             if d.viscosity is not None
             and d.ref_temp is not None]

    if len(dvisc) > 0:
        visc_table = []
        for visc_point in dvisc:
            v = visc_point.viscosity.converted_to(units).value
            t = visc_point.ref_temp.converted_to(temp_units).value
            visc_table.append((v, t))
        return visc_table

    kvisc = oil.sub_samples[0].physical_properties.kinematic_viscosities
    if len(kvisc) > 0:
        raise NotImplementedError("can't compute dynamic from kinematic yet")
        kvisc = get_kinematic_viscosity_data(oil)
        # convert here.
    else:
        return []


def convert_dvisc_to_kvisc(dvisc, density):
    """
    convert dynamic viscosity to kinematic viscosity

    :param density: an initialized Density object

    dvisc and densities are tables as returned from:

     - ``get_dynamic_viscosity_data``
     - ``get_density_data``

    units: viscosity: Pas or kg/(m s)
           density: kg/m^3
    """
    kvisc_table = []
    for (dv, temp) in dvisc:
        kv = dv / density.at_temp(temp)
        kvisc_table.append((kv, temp))
    return kvisc_table


# def density_at_temp(densities, temp, units="kg/m^3", temp_units="K"):
#     # sort them to make sure
#     densities = sorted(densities, key=itemgetter(1))
#     dens, temps = zip(*densities)

#     return np.interp(temp, temps, dens)

# def get_kinematic_viscosity_at_temp(temp,
#                                     kvis_units='cSt',
#                                     temp_units='C'):
#     raise NotImplementedError


def get_pour_point(oil):
    """
    Return oil's pour point or None
    """
    phys_props = oil.sub_samples[0].physical_properties

    pour_point = phys_props.pour_point

    return pour_point


def get_distillation_cuts(oil, units="fraction", temp_units="K"):
    """
    Return a table of distillation data:

    list of (cut fraction, temp) pairs

    :param oil: the oil object to get data from

    :param units="fraction": units you want the fraction in

    :param temp_units="K": units you want the temperature in

    """

    distillation_cuts = oil.sub_samples[0].distillation_data.cuts

    # create normalized list of densities
    cuts_table = []
    for cut in distillation_cuts:
        f = cut.fraction.converted_to(units).value
        t = cut.vapor_temp.converted_to(temp_units).value
        cuts_table.append((f, t))
    cuts_table.sort(key=itemgetter(0))
    return cuts_table


def max_water_fraction_emulsion(oil):
    if (oil.metadata.product_type == 'Crude Oil NOS' or
            oil.metadata.product_type == 'Bitumen Blend'):
        return 0.9
    else:
        return 0.0


def emul_water(oil):
    """
    This function computes two terms used in emulsification.
    Ymax is the maximum water fraction of a stable emulsion.
    Smax is the maximum surface area of the water droplets inside
    the emulsion. (from ADIOS2)
    """

    emulsion_max_water = None  # need to drill down in database for this info

    if emulsion_max_water is None:  # max water content not in library
        dens = Density(oil)
        density = dens.at_temp(288.15)
        kvis = KinematicViscosity(oil)
        viscosity = kvis.at_temp(288.15)
        dynamic_viscosity = viscosity * density
        if (dynamic_viscosity > 0.050):
            Ymax = 0.9 - 0.0952 * np.log(dynamic_viscosity / 0.050)
        else:
            Ymax = 0.9
    else:
        Ymax = emulsion_max_water  # stay with the value you have

    # this is done is py_gnome
    # drop_min = 1.0e-6		# min oil droplet size
    # Smax = (6.0 / drop_min) * (Ymax / (1.0 - Ymax))
    return Ymax


def bullwinkle_fraction(oil):

    Ni = 0
    Va = 0
    try:
        bulk_composition = oil.sub_samples[0].bulk_composition
        for element in bulk_composition:
            if element.name == "nickel":
                Ni = element.measurement.value
    except Exception:
        Ni = 0.

    try:
        bulk_composition = oil.sub_samples[0].bulk_composition
        for element in bulk_composition:
            if element.name == "vanadium":
                Va = element.measurement.value
    except Exception:
        Va = 0.

    try:
        f_asph = oil.sub_samples[0].SARA.asphaltenes.value
    except Exception:
        f_asph = 0.

    # need to go through subsamples checking for weathered data
    # if oil.metadata.emuls_constant_max is not None:
    # bullwinkle_fraction = emuls_constant_max

    if (oil.metadata.product_type != "Crude Oil NOS" and
            oil.metadata.product_type != "Bitumen Blend"):
        bullwinkle_fraction = 1.0
    else:
        oil_api = oil.metadata.API
        if (Ni > 0.0 and Va > 0.0 and Ni + Va > 15.0):
            bullwinkle_fraction = 0.0
        elif f_asph > 0:
            bullwinkle_fraction = 0.20219 - 0.168 * np.log10(f_asph)
            bullwinkle_fraction = np.clip(bullwinkle_fraction, 0.0, 0.303)
        elif oil_api < 26.0:
            bullwinkle_fraction = 0.08
        elif oil_api > 50.0:
            bullwinkle_fraction = 0.303
        else:
            bullwinkle_fraction = (-1.038 - 0.78935 * np.log10(1.0 / oil_api))

        bullwinkle_fraction = _adios2_new_bull_calc(bullwinkle_fraction,
                                                    oil_api)

    return bullwinkle_fraction


def _adios2_new_bull_calc(bullwinkle_fraction, oil_api):
    """
    From the Adios2 c++ file OilInitialize.cpp, there is functionality
    inside the function CAdiosData::Bullwinkle() which is annotated
    in the code as 'new bull calc'.

    It uses the following definitions:
    - TG, Documented as the value 'dT/df - evaporation'.
          I can only assume this is the initial fractional rate of
          evaporation.
    - TBP, Documented as the 'ADIOS 1 liquid boiling point
          (bubble pt)'.
    - BullAdios1, which appears to be used to scale-average the
                  initially computed bullwinkle fraction.
      get_density_data
    Regardless, in order to approximate what Adios2 is doing, we
    need this modification of our bullwinkle fraction.
    """

    t_g = 1356.7 - 247.36 * np.log(oil_api)
    t_bp = 532.98 - 3.1295 * oil_api
    bull_adios1 = (483.0 - t_bp) / t_g

    bull_adios1 = np.clip(bull_adios1, 0.0, 0.4)

    return 0.5 * (bullwinkle_fraction + bull_adios1)
