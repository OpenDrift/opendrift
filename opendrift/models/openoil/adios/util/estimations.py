'''
    Oil Properties Estimations.

    These are primitive estimation algorithms to be used primarily for
    estimating oil properties based on known measured values.
'''

import numpy as np


# fixme -- these are unit conversion, and are in the unit_conversion(NUCOS)
#          library.  We should use that lib to be consistent.

def density_from_api(api):
    '''
        Source: Adios2

        API is a measure of how heavy an oil is compared to water.
        So it is a different unit for specific gravity
    '''
    kg_m_3 = 141.5 / (131.5 + api) * 1000.0
    ref_temp_k = 273.15 + 15.0

    return kg_m_3, ref_temp_k


def api_from_density(density):
    '''
        Source: Adios2
    '''
    kg_m_3 = density

    return (141.5 / kg_m_3 * 1000.0) - 131.5


def density_at_temp(ref_density, ref_temp_k, temp_k, k_rho_t=0.0008):
    '''
        Source: Adios2

        If we have an oil density at a reference temperature,
        then we can estimate what its density might be at
        another temperature.

        NOTE: need a reference for the coefficient of expansion
    '''
    return ref_density / (1.0 - k_rho_t * (ref_temp_k - temp_k))


def vol_expansion_coeff(rho_0, t_0, rho_1, t_1):
    '''
        Calculate the volumetric expansion coefficient of a liquid
        based on a set of two densities and their associated temperatures.
    '''
    if t_0 == t_1:
        k_rho_t = 0.0
    else:
        k_rho_t = (rho_0 - rho_1) / (rho_0 * (t_1 - t_0))

    return k_rho_t


def specific_gravity(density):
    '''
        Specific Gravity of Oil with respect to water at 15C (definition used
        for API gravity)
    '''
    return density / 999.103


def dvis_to_kvis(dvis, density):
    '''
        Source: Definition of kinematic viscosity.

        Conversion from dynamic viscosity to kinematic viscosity.
    '''
    return dvis / density


def kvis_at_temp(ref_kvis, ref_temp_k, temp_k, k_v2=2416.0):
    '''
        Source: Adios2

        If we have an oil kinematic viscosity at a reference temperature,
        then we can estimate what its viscosity might be at
        another temperature.

        Note: Bill's most recent viscosity document, and an analysis of the
              multi-KVis oils in our oil library suggest that a value of
              2416.0 (Abu Eishah 1999) would be a good default value for k_v2.
    '''
    return ref_kvis * np.exp(k_v2 / temp_k - k_v2 / ref_temp_k)


def resin_fraction(density, viscosity, f_other=0.0):
    A = _A_coeff(density)
    B = _B_coeff(density, viscosity)

    f_res = 0.033 * A + 0.00087 * B - 0.74
    f_res = np.clip(f_res, 0.0, 1.0 - f_other)

    return f_res


def asphaltene_fraction(density, viscosity, f_other=0.0):
    A = _A_coeff(density)
    B = _B_coeff(density, viscosity)

    f_asph = (0.000014 * A ** 3.0 +
              0.000004 * B ** 2.0 -
              0.18)
    f_asph = np.clip(f_asph, 0.0, 1.0 - f_other)

    return f_asph


def saturates_fraction(density, viscosity, f_other=0.0):
    A = _A_coeff(density)
    B = _B_coeff(density, viscosity)

    f_sat = -2.5 + 76.6 / A + 0.00013 * np.log(B)
    f_sat = np.clip(f_sat, 0.0, 1.0 - f_other)

    return f_sat


def aromatics_fraction(f_res, f_asph, f_sat):
    f_arom = 1.0 - (f_res + f_asph + f_sat)
    f_arom = np.clip(f_arom, 0.0, 1.0)

    return f_arom


def _A_coeff(density):
    '''
        Source: Fingas empirical formulas that are based upon analysis
                of ESTC oil properties database.

        This is an intermediate calculation for a coefficient to be
        used to generate the inert mass fractions of an oil.
    '''
    return 10.0 * np.exp(0.001 * density)


def _B_coeff(density, viscosity):
    '''
        Source: Fingas empirical formulas that are based upon analysis
                of ESTC oil properties database.

        This is an intermediate calculation for a coefficient to be
        used to generate the inert mass fractions of an oil.
    '''
    return 10.0 * np.log(1000.0 * density * viscosity)


def cut_temps_from_api(api, N=5):
    '''
        Source: Adios2 & Jones R. (1997),
                A Simplified Pseudo-component Oil Evaporation Model,
                Proceedings of the 20th Arctic and Marine Oil Spill Program,
                Vancouver, CA,
                Vol. 1, pp. 43-62

        Generate distillation cut temperatures from the oil's API.
    '''
    T_0 = 457.0 - 3.34 * api
    T_G = 1357.0 - 247.7 * np.log(api)

    return np.array([(T_0 + T_G * i / N) for i in range(N)])


def fmasses_from_cuts(f_evap_i):
    '''
        Generate distillation cut fractional masses from the
        cumulative distillation fractions in the cut data.
    '''
    fmass_i = np.array(f_evap_i)
    fmass_i[1:] = np.diff(fmass_i)

    return fmass_i


def fmasses_flat_dist(f_res, f_asph, N=5):
    '''
        Generate a flat distribution of N distillation cut fractional masses.
    '''
    return np.array([(1.0 - f_res - f_asph) / N] * N)


def saturate_mol_wt(boiling_point):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 2.48 and table 2.6

        Note: for this to actually work in every case, we need to limit
              our temperature to:
              - T_i < 1070.0
              - T_i >
              - T_i > 1070.0 - exp(6.98291)  (roughly about == -8.06)
    '''
    T_i = np.clip(np.array(boiling_point),
                  1070.0 - np.exp(6.98291) + 0.00001,
                  1070.0 - 0.00001)
    return (49.677 * (6.98291 - np.log(1070.0 - T_i))) ** (3.0 / 2.0)


def aromatic_mol_wt(boiling_point):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 2.48 and table 2.6

        Note: for this to actually work in every case, we need to limit
              our temperature to:
              - T_i < 1015.0
              - T_i > 1015.0 - exp(6.911)  (roughly about == 11.76)
    '''
    T_i = np.clip(np.array(boiling_point),
                  1015.0 - np.exp(6.911) + 0.00001,
                  1015.0 - 0.00001)
    return (44.504 * (6.911 - np.log(1015.0 - T_i))) ** (3.0 / 2.0)


def resin_mol_wt(_boiling_point):
    '''
        Source: Recommendation from Bill Lehr

        Note: We pass in a boiling point to remain consistent with the other
              molecular weight functions, even though it is not used.
        Note: We return a scalar in all cases.  This should still work with
              numpy array operations, but probably not with regular Python
              sequence types.
              We can fix this if the need arises.
    '''
    return 800.0


def asphaltene_mol_wt(_boiling_point):
    '''
        Source: Recommendation from Bill Lehr

        Note: We pass in a boiling point to remain consistent with the other
              molecular weight functions, even though it is not used.
        Note: We return a scalar in all cases.  This should still work with
              numpy array operations, but probably not with regular Python
              sequence types.
              We can fix this if the need arises.
    '''
    return 1000.0


def trial_densities(boiling_points, watson_factor):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 2.13 and table 9.6

        Generate an initial estimate of volatile oil components based
        on boiling points and the Watson Characterization Factor.
        This is only good for estimating Aromatics & Saturates.
    '''
    boiling_points = np.array(boiling_points)
    return 1000.0 * (1.8 * boiling_points) ** (1.0 / 3.0) / watson_factor


def saturate_densities(boiling_points):
    K_w_sat = 12.0

    return trial_densities(boiling_points, K_w_sat)


def aromatic_densities(boiling_points):
    K_w_arom = 10.0

    return trial_densities(boiling_points, K_w_arom)


def resin_densities(_boiling_points):
    '''
        Note: We pass in a boiling point to remain consistent with the other
              molecular weight functions, even though it is not used.
        Note: We return a scalar in all cases.  This should still work with
              numpy array operations, but probably not with regular Python
              sequence types.
              We can fix this if the need arises.
    '''
    return 1100.0


def asphaltene_densities(_boiling_points):
    '''
        Note: We pass in a boiling point to remain consistent with the other
              molecular weight functions, even though it is not used.
        Note: We return a scalar in all cases.  This should still work with
              numpy array operations, but probably not with regular Python
              sequence types.
              We can fix this if the need arises.
    '''
    return 1100.0


def _hydrocarbon_characterization_param(specific_gravity, temp_k):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 2.115

        Source: Huang, P. K. (1977),
                Characterization and Thermodynamic Correlations for
                Undefined Hydrocarbon Mixtures,
                Ph.D. Dissertation
                Pennsylvania State University,
                University Park, PA,

        This is a characterization parameter, designated as I, that
        was first used by Huang to correlate hydrocarbon properties
    '''
    T_i = np.array(temp_k)
    SG_i = np.array(specific_gravity)

    return 0.3773 * T_i ** (-0.02269) * SG_i ** (0.9182)


def refractive_index(hc_char_param):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 2.114

        This is the refractive index of liquid hydrocarbons at 20C,
        correlated through parameter I
    '''
    i = np.array(hc_char_param)

    return ((1 + 2 * i) / (1 - i)) ** (1.0 / 2.0)


def _hydrocarbon_grouping_param(mol_wt, specific_gravity, temp_k):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 3.50

        This is a hydrocarbon grouping parameter correlated through
        the molecular weight and the refractive index.  Riazi claims
        that it:
        - separates paraffins and aromatics
        - identifies various hydrocarbon types
    '''
    i = _hydrocarbon_characterization_param(specific_gravity, temp_k)
    n = refractive_index(i)

    return mol_wt * (n - 1.475)


def saturate_mass_fraction(fmass_i, mol_wt, specific_gravity, temp_k):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eqs. 3.77 and 3.78

        For a petroleum hydrocarbon fraction in which the molecular
        weight, specific gravity, and boiling point are known, we can
        estimate the fraction of that substance which contains
        saturate compounds.
        When forming his equations, Riazi deals with the PNA composition,
        which essentially means Paraffins, Naphthenes, and Aromatics.
        And for our purposes, the saturates include the paraffins and
        naphthenes.
    '''
    SG_sat_i = specific_gravity
    m = _hydrocarbon_grouping_param(mol_wt, SG_sat_i, temp_k)

    X_P = 3.7387 - 4.0829 * SG_sat_i + 0.014772 * m
    X_N = -1.5027 + 2.10152 * SG_sat_i - 0.02388 * m

    f_sat_i = fmass_i * (X_P + X_N)

    return np.clip(f_sat_i, 0.0, fmass_i)


def oil_water_surface_tension_from_api(api):
    return 0.001 * (39.0 - 0.2571 * api)


def pour_point_from_kvis(ref_kvis, ref_temp_k):
    '''
        Source: Adios2

        If we have an oil kinematic viscosity at a reference temperature,
        then we can estimate what its pour point might be.
    '''
    c_v1 = 5000.0
    ref_kvis = np.array(ref_kvis)
    ref_temp_k = np.array(ref_temp_k)

    T_pp = (c_v1 * ref_temp_k) / (c_v1 - ref_temp_k * np.log(ref_kvis))

    return T_pp


def pour_point_from_sg_mw_kvis(specific_gravity, mol_wt, kvis):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 3.119

        Another way of estimating pour point.
        These inputs may not be available for most imported oil records.
    '''
    SG = specific_gravity

    return (130.47 *
            SG ** 2.970566 *
            mol_wt ** (0.61235 - 0.47357 * SG) *
            kvis * (0.310331 - 0.32834 * SG))


def flash_point_from_bp(temp_k):
    '''
        Source: Reference: Chang A., K. Pashakanti, and Y. Liu (2012),
                           Integrated Process Modeling and Optimization,
                           Wiley Verlag.

    '''
    temp_k = np.array(temp_k)
    return 117.0 + 0.69 * temp_k


def flash_point_from_api(api):
    '''
        Source: Reference: Chang A., K. Pashakanti, and Y. Liu (2012),
                           Integrated Process Modeling and Optimization,
                           Wiley Verlag.

    '''
    api = np.array(api)
    return 457.0 - 3.34 * api


def bullwinkle_fraction_from_asph(f_asph):
    '''
        Source: Adios2
    '''
    return 0.32 - 3.59 * f_asph


def bullwinkle_fraction_from_api(api):
    '''
        Source: Adios2
    '''
    return 0.5762 * np.log10(api) - 0.6353
