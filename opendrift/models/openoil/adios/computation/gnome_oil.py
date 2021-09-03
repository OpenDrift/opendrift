"""
Code for making a "GnomeOil" from an Oil Object

See the PyGNOME code for more about GNOME's requirements

NOTE: This make s JSON compatible Python structure from which to build
a GnomeOil
"""

import copy

import numpy as np

import unit_conversion as uc

import logging
logger = logging.getLogger(__name__)

from . import estimations as est
from .physical_properties import get_density_data, get_kinematic_viscosity_data, get_distillation_cuts
from .physical_properties import bullwinkle_fraction, max_water_fraction_emulsion
from .physical_properties import Density, KinematicViscosity
from .estimations import pour_point_from_kvis, flash_point_from_bp, flash_point_from_api


def get_empty_dict():
    """
    This provides an empty dictionary with everything that is needed
    to generate a GNOME Oil
    """
    return{"name": "",
           # Physical properties
           "api": None,
           "pour_point": None,
           "solubility": None,  # kg/m^3
           # emulsification properties
           "bullwinkle_fraction": None,
           "bullwinkle_time": None,
           "emulsion_water_fraction_max": None,
           "densities": [],
           "density_ref_temps": [],
           "density_weathering": [],
           "kvis": [],  # m/s^2
           "kvis_ref_temps": [],  # K
           "kvis_weathering": [],
           # PCs:
           "mass_fraction": [],
           "boiling_point": [],
           "molecular_weight": [],
           "component_density": [],
           "sara_type": [],
           "flash_point": None,
           "adios_oil_id": None,
           }


def make_gnome_oil(oil):
    """
    Make a dict that a GnomeOil can be built from

    A GnomeOil needs:

              "name,
              "# Physical properties
              "api,
              "pour_point,
              "solubility,  # kg/m^3
              "# emulsification properties
              "bullwinkle_fraction,
              "bullwinkle_time,
              "emulsion_water_fraction_max,
              "densities,
              "density_ref_temps,
              "density_weathering,
              "kvis,
              "kvis_ref_temps,
              "kvis_weathering,
              "# PCs:
              "mass_fraction,
              "boiling_point,
              "molecular_weight,
              "component_density,
              "sara_type,
              "flash_point=None,
              "adios_oil_id=None,

    """
    # metadata:
    go = get_empty_dict()
    # make sure we don't change the original oil object
    # oil = copy.deepcopy(oil)
    go['name'] = oil.metadata.name

    dens = Density(oil)
    ref_density = dens.at_temp(288.7)  # 60F in K
    go['api'] = uc.convert('kg/m^3', 'API', ref_density)
    # for gnome_oil we don't treat api as data, only api from density
    oil.metadata.API = go['api']
    go['adios_oil_id'] = oil.oil_id

    # Physical properties
    phys_props = oil.sub_samples[0].physical_properties

    flash_point = phys_props.flash_point
    if flash_point is None:
        go['flash_point'] = None
    else:
        fp = phys_props.flash_point.measurement.converted_to('K')
        if fp.max_value is not None:
            go['flash_point'] = fp.max_value
        elif fp.value is not None:
            go['flash_point'] = fp.value
        else:  # shouldn't happen, but ...
            go['flash_point'] = None

    pour_point = phys_props.pour_point
    if pour_point is None:
        go['pour_point'] = estimate_pour_point(oil)
    else:
        pp = phys_props.pour_point.measurement.converted_to('K')
        if pp.max_value is not None:
            go['pour_point'] = pp.max_value
        elif pp.value is not None:
            go['pour_point'] = pp.value
        elif pp.min_value is not None:
            go['pour_point'] = pp.min_value
        else:
            go['pour_point'] = None

    # fixme: We need to get the weathered densities, if they are there.
    densities = get_density_data(oil, units="kg/m^3", temp_units="K")

    go['densities'], go['density_ref_temps'] = zip(*densities)
    go['density_weathering'] = [0.0] * len(go['densities'])

    viscosities = get_kinematic_viscosity_data(oil, units="m^2/s", temp_units="K")

    if viscosities:
        go['kvis'], go['kvis_ref_temps'] = zip(*viscosities)
        go['kvis_weathering'] = [0.0] * len(go['kvis'])
    else:
        raise ValueError("Gnome oil needs at least one viscosity value")

    bullwinkle = None
    for sub_sample in oil.sub_samples:
        try:
            frac_weathered = sub_sample.metadata.fraction_weathered.converted_to('fraction').value
            if bullwinkle is None or frac_weathered > bullwinkle:
                bullwinkle = frac_weathered
        except:
            frac_weathered = None

    if bullwinkle is None:
        go['bullwinkle_fraction'] = bullwinkle_fraction(oil)
    else:
        go['bullwinkle_fraction'] = bullwinkle

    #go['bullwinkle_fraction'] = bullwinkle_fraction(oil)
    go['emulsion_water_fraction_max'] = max_water_fraction_emulsion(oil)
    go['solubility'] = 0
    go['k0y'] = 2.024e-06 #do we want this included?

    # pseudocomponents
    cut_temps, frac_evap = normalized_cut_values(oil)

    mass_fraction = component_mass_fractions(oil)
    mask = np.where(mass_fraction==0)
    mol_wt = np.delete(component_mol_wt(cut_temps),mask)
    comp_dens = np.delete(component_densities(cut_temps),mask)
    boiling_pt = np.delete(component_temps(cut_temps),mask)
    sara_type = np.delete(np.array(component_types(cut_temps)),mask)
    mass_frac = np.delete(mass_fraction,mask)

    go['molecular_weight'] = mol_wt.tolist()
    go['component_density'] = comp_dens.tolist()
    go['mass_fraction'] = mass_frac.tolist()
    go['boiling_point'] = boiling_pt.tolist()
    go['sara_type'] = sara_type.tolist()

    return go


def estimate_pour_point(oil):
    """
    estimate pour point from kinematic viscosity

    """

    pour_point = None
    kvis = get_kinematic_viscosity_data(oil)
    c_v1 = 5000.0
    if kvis:
        lowest_kvis = kvis[0]

        #pour_point = pour_point_from_kvis(lowest_kvis[0], lowest_kvis[1])
        pour_point = (c_v1 * lowest_kvis[1]) / (c_v1 - lowest_kvis[1] * np.log(lowest_kvis[0]))

    return pour_point


def estimate_flash_point(oil):
    """
    estimate flash point from api or boiling point

    """

    oil_api = oil.metadata.API

    cuts = get_distillation_cuts(oil)

    # fixme: if we do need this, we should have a better way to get
    #        boiling point -- the first cut is not necessarily the BP!
    #        the IBP is stored in the distillation record, if it is known.
    if len(cuts) > 2:
        lowest_cut = cuts[0]
        flash_point = flash_point_from_bp(lowest_cut[1])
    elif oil_api is not None:
        flash_point = flash_point_from_api(oil_api)
    else:
        flash_point = None

    return flash_point


def component_temps(cut_temps, N=10):
    """
    component temps from boiling point

    """
    component_temps = np.append(cut_temps, list(zip(cut_temps, cut_temps, cut_temps, cut_temps)))
    len_ct = len(cut_temps)
    new_temps = component_temps[len_ct:].copy()

    return np.asarray(new_temps)

def component_types(cut_temps, N=10):
    """
    set component SARA types

    """
    T_i = component_temps(cut_temps, N)

    types_out = ['Saturates', 'Aromatics', 'Resins', 'Asphaltenes'] * int((len(T_i) / 4))

    return types_out


def component_densities(boiling_points):
    """
    estimate component densities from boiling point

    """
    rho_list = list(zip(est.saturate_densities(boiling_points),
                  est.aromatic_densities(boiling_points),
                  est.resin_densities(boiling_points),
                  est.asphaltene_densities(boiling_points)))

    return (np.asarray(rho_list)).flatten()


def component_mol_wt(boiling_points):
    """
    estimate component molecular weight from boiling point

    """
    mw_list = list(zip(est.saturate_mol_wt(boiling_points),
                  est.aromatic_mol_wt(boiling_points),
                  est.resin_mol_wt(boiling_points),
                  est.asphaltene_mol_wt(boiling_points)))

    return (np.asarray(mw_list)).flatten()


def inert_fractions(oil, density=None, viscosity=None):
    """
    resins and asphaltenes from database or estimated if None

    """
    f_res = None
    f_asph = None
    estimated_res = estimated_asph = False

    f_res = oil.sub_samples[0].SARA.resins
    f_asph = oil.sub_samples[0].SARA.asphaltenes

    if f_res is not None:
        f_res = f_res.converted_to('fraction').value
    if f_asph is not None:
        f_asph = f_asph.converted_to('fraction').value

    if f_res is not None and f_asph is not None:
        return f_res, f_asph, estimated_res, estimated_asph
    else:
        dens = Density(oil)
        density = dens.at_temp(288.15)
        kvis = KinematicViscosity(oil)
        viscosity = kvis.at_temp(288.15)

    if f_res is None:
        f_res = est.resin_fraction(density, viscosity)
        estimated_res = True

    if f_asph is None:
        f_asph = est.asphaltene_fraction(density, viscosity, f_res)
        estimated_asph = True

    return f_res, f_asph, estimated_res, estimated_asph

#this will all be removed...
def _linear_curve(x, a, b):
    '''
        Here we describe the form of a linear function for the purpose of
        curve-fitting measured data points.
    '''
    return (a * x + b)


def clamp(x, M, zeta=0.03):
    '''
        We make use of a generalized logistic function or Richard's curve
        to generate a linear function that is clamped at x == M.
        We make use of a zeta value to tune the parameters nu, resulting in a
        smooth transition as we cross the M boundary.
    '''
    return (x -
            (x / (1.0 + np.e ** (-15 * (x - M))) ** (1.0 / (1 + zeta))) +
            (M / (1.0 + np.e ** (-15 * (x - M))) ** (1.0 / (1 - zeta))))


def _inverse_linear_curve(y, a, b, M, zeta=0.12):
    y_c = clamp(y, M, zeta)

    return (y_c - b) / a


# adios version
def normalized_cut_values_adios(oil):
    oil_api = oil.metadata.API
    T0=457.16-3.3447*oil_api
    TG=1356.7-247.36*np.log(oil_api)
    cuts = get_distillation_cuts(oil)
    if len(cuts) == 0:
        # should be a warning if api < 50 or not a crude
        #NumComp=5					 #number of components
        NumComp=10					 #number of components
        BP = np.array([(T0 + (TG *(i+.5)) / NumComp) for i in range(NumComp)])
        fevap = np.cumsum(est.fmasses_flat_dist(N=10))
    else:
        NumCuts = len(cuts)
        NumComp = NumCuts + 1
        BP_i, fevap_i = list(zip(*[(c[1], c[0]) for c in cuts]))
        BP = np.zeros(NumComp)
        fevap = np.zeros(NumComp)
        len_BP = len(BP)
        T1 = BP_i[0]
        T2 = BP_i[1]
        TN = BP_i[NumCuts-1]
        fN = fevap_i[NumCuts-1]
        BP[0] = T1-.5*fevap_i[0]*(T2-T1)/(fevap_i[1]-fevap_i[0])
        fevap[0] = fevap_i[0]
        for i in range(NumCuts-1):
            BP[i+1] = .5*(BP_i[i]+BP_i[i+1])
            fevap[i+1] = fevap_i[i+1]
        BPN =  TN+.5*(1-fN)*(TN-T1)/(fN-fevap_i[0])
        fevap[NumComp-1] = 1
        BP[NumComp-1] = TN+.5*(1-fN)*(TN-T1)/(fN-fevap_i[0])
        if BP[NumComp-1]> 1500.:
            BP[NumComp-1] = 1500
        if BP[NumComp-1]< 1050:
            BP[NumComp-1] = 1050

    return np.asarray(BP), est.fmasses_from_cuts(fevap)

def normalized_cut_values(oil):
    f_res = f_asph = 0 # for now, we are including the resins and asphaltenes
    cuts = get_distillation_cuts(oil)
    oil_api = oil.metadata.API
    iBP = 266
    tBP = 1050

    if len(cuts) == 0:
        # should be a warning if api < 50 or not a crude
        oil_api = oil.metadata.API
        if oil.metadata.product_type != 'Crude Oil NOS':
            logger.error('warning: oil not recommended for use in Gnome')
            # print(WARNINGS['W007'] + "  - oil not recommended for use in Gnome")
        if oil_api < 0:
            raise ValueError("Density is too large for estimations. Oil not suitable for use in Gnome")

        # ADIOS2 method
        BP_i = est.cut_temps_from_api(oil_api)
        fevap_i = np.cumsum(est.fmasses_flat_dist(f_res, f_asph))
        # Robert's new method
        #iBP = 10/9 * (519.3728 - 3.6637 * oil_api) - 1015 / 9
        #tBP = 1015
        #BP_i = [iBP, tBP]
        #fevap_i = [0,1]
    else:
        BP_i, fevap_i = list(zip(*[(c[1], c[0]) for c in cuts]))
        N = len(BP_i)
        if not (fevap_i[1] == fevap_i[0]):
            iBP = BP_i[0] - fevap_i[0] * (BP_i[1] - BP_i[0]) / (fevap_i[1] - fevap_i[0])
        if not (fevap_i[N-1] == fevap_i[0]):
            tBP = BP_i[N-1] + (1 - fevap_i[0]) * (BP_i[N-1] - BP_i[0]) / (fevap_i[N-1] - fevap_i[0])

    iBP = max(266, iBP)
    tBP = min(1050, tBP)

    set_temp = [266,310,353,483,563,650,800,950,1050]
    #set_temp = [266,310,353,423,483,523,543,563,650,800,950,1050]

    N = len(BP_i)

    new_fevap = fevap_i
    new_BP = BP_i
    if not fevap_i[N-1] == 1:
        new_BP = np.append(BP_i,tBP)
        new_fevap = np.append(fevap_i, 1.0)
    if not fevap_i[0] == 0:
        new_BP = np.insert(new_BP,0,iBP)
        new_fevap = np.insert(new_fevap, 0, 0)
    new_evap = np.interp(set_temp, new_BP, new_fevap)

    if new_evap[-1] < 1:
        new_evap[-1] = 1	# put all the extra mass in the last cut

    avg_evap_i = np.asarray(new_evap[1:])
    avg_temp_i = np.asarray([(set_temp[i]+set_temp[i+1])/2 for i in range(0, len(set_temp)-1)])

    avg_temp_i = avg_temp_i[avg_evap_i!=0]
    avg_evap_i = avg_evap_i[avg_evap_i!=0]

    num_ones = 0
    for i in range(len(avg_evap_i)):
        if avg_evap_i[i] == 1:
            num_ones += 1

    if num_ones > 1:
        for i in range(num_ones-1):
            avg_temp_i = np.delete(avg_temp_i,len(avg_evap_i)-1)
            avg_evap_i = np.delete(avg_evap_i,len(avg_evap_i)-1)

    return np.asarray(avg_temp_i), est.fmasses_from_cuts(avg_evap_i)

def normalized_cut_values_james(oil, N=10):
    """
    estimate cut temperatures

    """
    from scipy.optimize import curve_fit	#temporary
    #f_res, f_asph, _estimated_res, _estimated_asph = inert_fractions(oil)
    f_res = f_asph = 0
    cuts = get_distillation_cuts(oil)
    oil_api = oil.metadata.API
    if len(cuts) == 0:
        if oil.metadata.product_type != 'Crude Oil NOS':
            logger.error('warning: oil not recommended for use in Gnome')
            # print(WARNINGS['W007'] + "  - oil not recommended for use in Gnome")
        if oil_api < 0:
            raise ValueError("Density is too large for estimations. Oil not suitable for use in Gnome")
        BP_i = est.cut_temps_from_api(oil_api)
        fevap_i = np.cumsum(est.fmasses_flat_dist(f_res, f_asph))
    else:
        BP_i, fevap_i = list(zip(*[(c[1], c[0]) for c in cuts]))

    popt, _pcov = curve_fit(_linear_curve, BP_i, fevap_i)
    f_cutoff = _linear_curve(732.0, *popt)  # center of asymptote (< 739)
    popt = popt.tolist() + [f_cutoff]

    fevap_i = np.linspace(0.0, 1.0 - f_res - f_asph, (N * 2) + 1)[1:]
    T_i = _inverse_linear_curve(fevap_i, *popt)

    fevap_i = fevap_i.reshape(-1, 2)[:, 1]
    T_i = T_i.reshape(-1, 2)[:, 0]

    above_zero = T_i > 0.0
    T_i = T_i[above_zero]
    fevap_i = fevap_i[above_zero]

    return T_i, est.fmasses_from_cuts(fevap_i)

# need to replace this
def component_mass_fractions(oil):
    """
    estimate pseudocomponent mass fractions

    """
    cut_temps, fmass_i = normalized_cut_values(oil)
    measured_sat = oil.sub_samples[0].SARA.saturates
    sat, arom, res, asph = sara_totals(oil)
    if measured_sat is not None:
        f_sat_i = est.saturate_mass_fraction(fmass_i,cut_temps,sat)	#scale by data value
    else:
        f_sat_i = est.saturate_mass_fraction(fmass_i,cut_temps)
    f_arom_i = fmass_i - f_sat_i

    f_res_i = np.zeros_like(f_sat_i)
    f_asph_i = np.zeros_like(f_sat_i)

    if asph + res > f_arom_i.sum(): # something went wrong
        res = 0
        asph = 0

    f_res_i[len(f_res_i)-1] = res
    f_asph_i[len(f_asph_i)-1] = asph
    f_inert = asph + res

    for i in range(len(f_arom_i)-1,-1,-1):
        if f_inert > 0:
            if f_arom_i[i] > f_inert:
                f_arom_i[i] = f_arom_i[i] - f_inert
                f_inert = 0
            else:
                f_inert = f_inert - f_arom_i[i]
                f_arom_i[i] = 0
        else:
            f_inert = f_inert * 1

    mf_list = list(zip(f_sat_i, f_arom_i, f_res_i, f_asph_i))	# will want to zip all 4 together
    return (np.asarray(mf_list)).flatten()

def sara_totals(oil):
    """
    get SARA from database
    estimate if no data

    """
    SARA = oil.sub_samples[0].SARA
    aromatics = oil.sub_samples[0].SARA.aromatics
    saturates = oil.sub_samples[0].SARA.saturates
    resins = oil.sub_samples[0].SARA.resins
    asphaltenes = oil.sub_samples[0].SARA.asphaltenes

    dens = Density(oil)
    density = dens.at_temp(288.15)
    kvis = KinematicViscosity(oil)
    viscosity = kvis.at_temp(288.15)

    if resins is None:
        resins_total = est.resin_fraction(density, viscosity)
    else:
        resins_total = resins.converted_to('fraction').value
    if asphaltenes is None:
        asphaltenes_total = est.asphaltene_fraction(density, viscosity, resins_total)
    else:
        asphaltenes_total = asphaltenes.converted_to('fraction').value
    if saturates is None:
        saturates_total = est.saturates_fraction(density, viscosity)
    else:
        saturates_total = saturates.converted_to('fraction').value
    if aromatics is None:
        aromatics_total = est.aromatics_fraction(resins_total,asphaltenes_total,saturates_total)
    else:
        aromatics_total = aromatics.converted_to('fraction').value

    #print("SARA totals = ", saturates_total,aromatics_total,resins_total, asphaltenes_total)

    return saturates_total, aromatics_total, resins_total, asphaltenes_total
