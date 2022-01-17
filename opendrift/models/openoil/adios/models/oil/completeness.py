'''
    Calculate the completeness of the data contained in an oil record.

    The top-level interface is a single function:
        completeness(pyjson_of_an_oil_record)

    (pyjson is a python dict that is a match for the JSON)

    It performs calculations that were designed by RobertJ, and returns a value
    with a scale of 0->100%.
'''
import logging

from .oil import Oil

# from ..common.measurement import MassFraction, Temperature

logger = logging.getLogger(__name__)


# is this function needed???
# answer: SRP, it's the right thing to do
def set_completeness(oil):
    oil.metadata.model_completeness = completeness(oil)


def completeness(oil):
    '''
        Calculate the completeness of the data contained in an oil record.

        :param oil: The oil record to be validated, in json-compatible python
                    data structure.
    '''
    res = 0
    for check_func in CHECKS:
        res += check_func(oil)

    return round(res * 10.0)


def check_emulsion_water_content(oil):
    '''
        One emulsion water content in any subsample. Score = 2.5
    '''
    sub_samples = oil.sub_samples

    for sample in sub_samples:
        emuls = sample.environmental_behavior.emulsions
        for e in emuls:
            if (is_measurement_good(e.water_content)):
                return 2.5  # only need one valid emulsion water content

    return 0.0


def check_density(oil):
    '''
        Fresh oil: One density or API. Score = 1
    '''

    if oil.metadata.API is not None:
        return 1.0

    if len(oil.sub_samples) > 0:
        ss = oil.sub_samples[0]

        densities = ss.physical_properties.densities

        for d in densities:
            if (is_measurement_good(d.density) and
                    is_measurement_good(d.ref_temp)):
                return 1.0

    return 0.0


def check_second_density(oil):
    '''
    Fresh oil: Second density separated by temperature.

    Score = deltaT/40 but not greater than 0.5

    maxDeltaT: The difference between the lowest and highest
    measurement in the set.
    '''
    if len(oil.sub_samples) > 0:
        ss = oil.sub_samples[0]
        densities = ss.physical_properties.densities

        temps = [d.ref_temp.converted_to('C').value
                 for d in densities
                 if d.ref_temp is not None]

        if len(temps) >= 2:
            t1, *_, t2 = sorted([t for t in temps if t is not None])
            delta_t = t2 - t1

            if delta_t > 0.0:
                return min(delta_t / 40.0, 0.5)

    return 0.0


def check_viscosity(oil):
    '''
    Fresh oil: One viscosity. Score = 0.5
    '''
    if len(oil.sub_samples) > 0:
        ss = oil.sub_samples[0]
        kvis = ss.physical_properties.kinematic_viscosities
        dvis = ss.physical_properties.dynamic_viscosities

        for v in kvis:
            if (is_measurement_good(v.viscosity) and
                    is_measurement_good(v.ref_temp)):
                return 0.5

        for v in dvis:
            if (is_measurement_good(v.viscosity) and
                    is_measurement_good(v.ref_temp)):
                return 0.5

    return 0.0


def check_second_viscosity(oil):
    '''
    Fresh oil: Second viscosity at a different temperature.

    Score = maxDeltaT/40, but not greater than 0.5

    maxDeltaT: The difference between the lowest and highest
    measurement in the set.
    '''
    if len(oil.sub_samples) > 0:
        ss = oil.sub_samples[0]
        kvis = ss.physical_properties.kinematic_viscosities
        dvis = ss.physical_properties.dynamic_viscosities

        temps = []
        for v_i in (kvis, dvis):
            temps.extend([v.ref_temp.converted_to('C').value
                          for v in v_i
                          if v.ref_temp is not None])

        if len(temps) >= 2:
            t1, *_, t2 = sorted([t for t in temps if t is not None])
            delta_t = t2 - t1

            if delta_t > 0.0:
                return min(delta_t / 40.0, 0.5)

    return 0.0


def check_distillation(oil):
    '''
    Fresh oil: Two Distillation cuts separated by mass or volume fraction.

    Score = 3 * maxDeltaFraction

    maxDeltaFraction: The difference between the lowest and
    highest measurement in the set
    '''
    if len(oil.sub_samples) > 0:
        ss = oil.sub_samples[0]
        cuts = ss.distillation_data.cuts

        fractions = [c.fraction.converted_to('fraction').value for c in cuts]

        if len(fractions) >= 2:
            f1, *_, f2 = sorted([f for f in fractions if f is not None])

            return 3.0 * (f2 - f1)

    return 0.0


def check_weathered_density(oil):
    '''
        One Evaporated oil: Density. Score = 1
    '''
    if len(oil.sub_samples) > 1:
        ss = oil.sub_samples[1]
        densities = ss.physical_properties.densities

        if (len(densities) > 0 and
                is_measurement_good(densities[0].density) and
                is_measurement_good(densities[0].ref_temp)):
            return 1.0

    return 0.0


def check_weathered_viscosity(oil):
    '''
        One Evaporated oil: Viscosity. Score = 1
    '''
    if len(oil.sub_samples) > 1:
        ss = oil.sub_samples[1]
        kvis = ss.physical_properties.kinematic_viscosities
        dvis = ss.physical_properties.dynamic_viscosities

        if (len(kvis) > 0 and
                is_measurement_good(kvis[0].viscosity) and
                is_measurement_good(kvis[0].ref_temp)):
            return 1.0

        if (len(dvis) > 0 and
                is_measurement_good(dvis[0].viscosity) and
                is_measurement_good(dvis[0].ref_temp)):
            return 1.0

    return 0.0


def is_measurement_good(measurement):
    return not any([(getattr(measurement, a, None) is None)
                    for a in ('value', 'unit', 'unit_type')])


# build a list of all the check function:
CHECKS = [val for name, val in vars().items() if name.startswith("check_")]
