"""
Code for estimating various properties that aren't in a oil record.

These can be used to generate complete oil specifications for
modeling, etc.

Note: there is code in here to make various estimations of
    SARA fractionations -- be warned, these are of unceratain
    use, and not very accurate

"""

import numpy as np

# This is  unit conversion -- should be handled by pynucos
# def api_from_density(density):

#     '''
#         Source: Adios2
#     '''
#     kg_m_3 = density

#     return (141.5 / kg_m_3 * 1000.0) - 131.5


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

    return np.array([(T_0 + (T_G * i) / N) for i in range(N)])


def fmasses_flat_dist(f_res=0, f_asph=0, N=5):
    '''
        Generate a flat distribution of N distillation cut fractional masses.
    '''
    return np.array([(1.0 - f_res - f_asph) / N] * N)


def fmasses_from_cuts(f_evap_i):
    '''
        Generate distillation cut fractional masses from the
        cumulative distillation fractions in the cut data.
    '''
    fmass_i = np.array(f_evap_i)
    fmass_i[1:] = np.diff(fmass_i)

    return fmass_i


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
                of ESTS oil properties database.

        This is an intermediate calculation for a coefficient to be
        used to generate the mass fractions of an oil.
    '''
    return 10.0 * np.exp(0.001 * density)


def _B_coeff(density, viscosity):
    '''
        Source: Fingas empirical formulas that are based upon analysis
                of ESTC oil properties database.

        This is an intermediate calculation for a coefficient to be
        used to generate the mass fractions of an oil.
    '''
    return 10.0 * np.log(1000.0 * density * viscosity)


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


def resin_mol_wt(boiling_points):
    '''
        Source: Recommendation from Bill Lehr
    '''
    return np.full_like(boiling_points, 800.0)


def asphaltene_mol_wt(boiling_points):
    '''
        Source: Recommendation from Bill Lehr
    '''
    return np.full_like(boiling_points, 1000.0)


def trial_densities(boiling_points, watson_factor):
    '''
        Source: Dr. M. R. Riazi,
                Characterization and Properties of Petroleum Fractions
                eq. 2.13 and table 9.6

        Generate an initial estimate of volatile oil components based
        on boiling points and the Watson Characterization Factor.
        This is only good for estimating Aromatics & Saturates.
    '''
    rho_i = 1000.0 * (1.8 * boiling_points) ** (1.0 / 3.0) / watson_factor

    return np.clip(np.array(rho_i), 0, 1090)


def saturate_densities(boiling_points):
    K_w_sat = 12.0

    return trial_densities(boiling_points, K_w_sat)


def aromatic_densities(boiling_points):
    K_w_arom = 10.0

    return trial_densities(boiling_points, K_w_arom)


def resin_densities(boiling_points):
    resin_rho =  np.full(len(boiling_points),1100.0)
    return resin_rho


def asphaltene_densities(boiling_points):
    return np.full_like(boiling_points, 1100.0)


def saturate_mass_fraction(fmass_i, temp_k, total_sat=None):
    '''
        Source: Dr. Robert Jones, based on average of 51 Exxon oils
        This assumes we do not known the SARA totals for the oil

    '''
    T_i = temp_k
    k = .0877

    if total_sat is not None:
#         sat_pct_temp = (124.1069 - k * T_i)
#         sat_pct_i = np.clip(sat_pct_temp,0,100)
#         f_sat_i = fmass_i * sat_pct_i / 100.
        k = (124.1069 *fmass_i.sum() - total_sat*100) / (fmass_i*T_i).sum()
        sat_pct_i = 124.1069 - k * T_i
        f_sat_i = fmass_i * sat_pct_i / 100.
#         scale = total_sat / f_sat_i.sum()
#         f_sat_i = scale * f_sat_i
#         sat_pct_i = f_sat_i / fmass_i
    else:
        sat_pct_i = 124.1069 - k * T_i
        f_sat_i = fmass_i * sat_pct_i / 100.

    #sat_pct_i = 124.1069 - k * T_i
    #f_sat_i = fmass_i * sat_pct_i / 100.

    return np.clip(f_sat_i, 0.0, fmass_i)
