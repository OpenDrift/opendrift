'''
Add labels to a record -- for those that we can use some automatic criteria for

For a given product type:
    some labels can be determined from API and/or Viscosity ranges.

The criteria follows the ASTM (and others) standards, where we can
'''
from math import inf

from ..product_type import types_to_labels
from ....computation.physical_properties import KinematicViscosity

# # Here are all the Product Types:
# ('Crude Oil NOS',
#  'Tight Oil',
#  'Condensate',
#  'Bitumen Blend',
#  'Refined Product NOS',
#  'Fuel Oil NOS',
#  'Distillate Fuel Oil',
#  'Residual Fuel Oil',
#  'Refinery Intermediate',
#  'Solvent',
#  'Bio-fuel Oil',
#  'Bio-Petro Fuel Oil',
#  'Natural Plant Oil',
#  'Lube Oil',
#  'Dielectric Oil',
#  'Other')

# # These are the current labels that aren't mapped yet:

# 'MDO', 'Vacuum Gas Oil'

# these are the labels with no criteria for density or viscosity
# e.g, if it's a "Crude Oil NOS", it's a 'Crude Oil'
synonyms_for_product_types = {
    'Crude Oil', 'Shale Oil', 'Fracking Oil', 'Fuel Oil', 'Residual Fuel', 'Distillate Fuel',
    'Refined Product', 'Condensate', 'Transformer Oil'
}
# If it's an exact match, then it's definitely a synonym
for pt, labels in types_to_labels.labels.items():
    for label in labels:
        if label == pt:
            # print(f"adding: {label}")
            synonyms_for_product_types.add(label)

# these are labels that are synonymous to other labels
synonyms_for_labels = {
    'Heavy Fuel Oil': ['HFO', 'No. 6 Fuel Oil', 'Bunker C'],
    'Kerosene': ['Jet Fuel'],
    'No. 2 Fuel Oil': ['Diesel', 'Home Heating Oil'],
}

no_criteria = {
    "api_min": -inf,
    "api_max": inf,
    "kvis_min": -inf,
    "kvis_max": inf,
    'kvis_temp': 15
}
label_map = {label: no_criteria for label in synonyms_for_product_types}

# this maps the labels according to API and kinematic viscosity (cSt at given temp in C) ranges.
label_map.update({
    # 'example_label': {"api_min": -inf, "api_max": inf, "kvis_min": -inf, "kvis_max": inf, 'kvis_temp': None},
    'Condensate': {
        "api_min": 50,
        "api_max": inf,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 15
    },
    'Light Crude': {
        "api_min": 35,
        "api_max": 50,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 15
    },
    'Medium Crude': {
        "api_min": 20,
        "api_max": 35,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 15
    },
    'Heavy Crude': {
        "api_min": -inf,
        "api_max": 20,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 15
    },
    'Group V': {
        "api_min": -inf,
        "api_max": 10.0,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 15
    },
    'Heavy Fuel Oil': {
        "api_min": -inf,
        "api_max": 15.0,
        "kvis_min": 200,
        "kvis_max": inf,
        'kvis_temp': 50
    },

    # pretty much made this up ... non-newtonian
    'Bitumen': {
        "api_min": -inf,
        "api_max": 10,
        "kvis_min": 1000,
        "kvis_max": inf,
        'kvis_temp': 40
    },
    # I went through all the bitumen records that we have, and most of them do not have
    # kinematic viscosity measurements. However, they all have dynamic viscosity
    # measurements which were > 10^6 cP at 0C and > 10^5 cP at 15C.
    # I propose changing the criteria to reflect these measurements.

    #***** The following need more verification
    # Refined light products
    'No. 2 Fuel Oil': {
        "api_min": 30,
        "api_max": 39,
        "kvis_min": 2.5,
        "kvis_max": 4,
        'kvis_temp': 38
    },
    'Kerosene': {
        "api_min": 47.6,
        "api_max": 67.8,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 38
    },
    'Aviation Gas': {
        "api_min": 47.6,
        "api_max": 70.8,
        "kvis_min": -inf,
        "kvis_max": inf,
        'kvis_temp': 38
    },
    'Gasoline': {
        "api_min": 59.7,
        "api_max": 76.6,
        "kvis_min": -inf,
        "kvis_max": 2.5,
        'kvis_temp': 38
    },

    # Intermediate Fuel Oils
    'MDO': {
        "api_min": 30,
        "api_max": 42,
        "kvis_min": -inf,
        "kvis_max": 11,
        'kvis_temp': 40
    },
    'IFO': {
        "api_min": 15,
        "api_max": 30,
        "kvis_min": 4,
        "kvis_max": 200,
        'kvis_temp': 38
    },
    #IFO needs an additional limitation to not be tied to biodiesel
})

for label, synonyms in synonyms_for_labels.items():
    label_map.update({syn: label_map[label] for syn in synonyms})


def get_suggested_labels(oil):
    """
    get the labels suggested for this oil

    :param oil: the oil object to get the labels for
    :type oil: Oil object
    """
    labels = set()
    pt = oil.metadata.product_type
    # everything gets its product type as a label as well
    # unless it has no product type
    # if pt:
    #     labels.add(pt)
    if pt == "Other":  # we don't want any labels auto added for Other
        return sorted(labels)
    try:
        for label in types_to_labels.left[oil.metadata.product_type]:
            if is_label(oil, label):
                labels.add(label)
    except KeyError:
        pass

    # sorting so they'll be in consistent order
    return sorted(labels)


def add_labels_to_oil(oil):
    """
    add labels to the passed in oil.

    this adds labels, but does not remove any existing ones

    :param oil: the oil object to add labels to
    :type oil: Oil object
    """
    labels = set(oil.metadata.labels)
    for label in types_to_labels.left['oil.metadata.product_type']:
        if is_label(oil, label):
            labels.add(label)
    oil.metadata.labels = sorted(labels)


def is_label(oil, label):
    try:
        data = label_map[label]
    except KeyError:
        return False

    api = oil.metadata.API

    # check API:
    if ((data['api_min'] != -inf) or (data['api_max'] != inf)):
        if api is None:
            is_label = False
        else:
            is_label = True if data['api_min'] <= api < data['api_max'] else False
    else:
        is_label = True

    if is_label and ((data['kvis_min'] != -inf) or
                     (data['kvis_max'] != inf)):  # check viscosity limits
        try:
            KV = KinematicViscosity(oil)
            kvis = KV.at_temp(temp=data['kvis_temp'], kvis_units='cSt', temp_units='C')
            is_label = True if data['kvis_min'] <= kvis < data['kvis_max'] else False
        except (ZeroDivisionError,
                ValueError):  # if it can't do this, we don't apply the label
            is_label = True

    return is_label


# def link_oil_to_labels(oil):
#     '''
#         Here, we have a single oil and we would like to link it to one or more
#         labels based on its properties.
#     '''
#     try:
#         labels = oil['metadata']['labels']
#     except TypeError:
#         labels = oil.metadata.labels
#     except KeyError:
#         labels = []

#     sample = OilEstimation(oil).get_sample()

#     if sample is None:
#         logger.warn('oil: {} has no fresh sample.  Skipping Categorization.'
#                     .format(oil['oil_id']))
#         return

#     if is_crude_light(oil):
#         labels.extend(['Crude Oil', 'Light Crude'])

#     if is_crude_medium(oil):
#         labels.extend(['Crude Oil', 'Medium Crude'])

#     if is_crude_heavy(oil):
#         labels.extend(['Crude Oil', 'Heavy Crude'])

#     if is_refined_fuel_oil_1(oil, sample):
#         labels.extend(['Refined Product',
#                        'Light',
#                        'No. 1 Fuel Oil',
#                        'Gasoline',
#                        'Kerosene'])

#     if is_refined_fuel_oil_2(oil, sample):
#         labels.extend(['Refined Product',
#                        'No. 2 Fuel Oil',
#                        'Diesel',
#                        'Heating Oil'])

#     if is_refined_ifo(oil, sample):
#         labels.extend(['Refined Product',
#                        'Intermediate',
#                        'Fuel Oil'])

#     if is_refined_fuel_oil_6(oil, sample):
#         labels.extend(['Refined Product',
#                        'Heavy',
#                        'Fuel Oil 6',
#                        'HFO',
#                        'Bunker',
#                        'Group V'])

#     if is_generic(oil):
#         labels.extend(['Other', 'Generic'])

#     if len(labels) == 0:
#         labels.extend(['Other'])

#     try:
#         oil['metadata']['labels'] = labels
#     except TypeError:
#         oil.metadata.labels = labels

# def is_crude(oil):
#     try:
#         return ('product_type' in oil['metadata'] and
#                 oil['metadata']['product_type'] is not None and
#                 oil['metadata']['product_type'].lower() == 'crude')
#     except KeyError:
#         return (hasattr(oil.metadata, 'product_type') and
#                 oil.metadata.product_type is not None and
#                 oil.metadata.product_type.lower() == 'crude')

# def is_refined(oil):
#     try:
#         return ('product_type' in oil['metadata'] and
#                 oil['metadata']['product_type'] is not None and
#                 oil['metadata']['product_type'].lower() == 'refined')
#     except KeyError:
#         return (hasattr(oil.metadata, 'product_type') and
#                 oil.metadata.product_type is not None and
#                 oil.metadata.product_type.lower() == 'refined')

# def api_min(oil, oil_api):
#     api = oil['metadata'].get('API', None)

#     return api is not None and api > oil_api

# def api_max(oil, oil_api):
#     api = oil['metadata'].get('API', None)

#     return api is not None and api < oil_api

# def is_crude_light(oil):
#     return is_crude(oil) and api_min(oil, 31.1)

# def is_crude_medium(oil):
#     return is_crude(oil) and api_max(oil, 31.1) and api_min(oil, 22.3)

# def is_crude_heavy(oil):
#     return is_crude(oil) and api_max(oil, 22.3)

# def is_refined_light_products(oil):
#     '''
#        Category Name:
#        - Light Products
#        Parent:
#        - Refined
#        Sample Oils:
#        - Cooper Basin Light Naphtha
#        - kerosene
#        - JP-4
#        - avgas
#        Density Criteria:
#        - API >= 35
#        Kinematic Viscosity Criteria:
#        - v > 0.0 cSt @ 38 degrees Celcius
#     '''
#     raise NotImplementedError

# def is_refined_fuel_oil_1(oil, sample):
#     '''
#        Category Name:
#        - Fuel oil #1/gasoline/kerosene
#        Sample Oils:
#        - gasoline
#        - kerosene
#        - JP-4
#        - avgas
#        Density Criteria:
#        - API >= 35
#        Kinematic Viscosity Criteria:
#        - v <= 2.5 cSt @ 38 degrees Celcius
#     '''
#     return (is_refined(oil) and
#             api_min(oil, 35.0) and
#             is_within_viscosity_range(sample, kvis_max=2.5))

# def is_refined_fuel_oil_2(oil, sample):
#     '''
#        Category Name:
#        - Fuel oil #2/Diesel/Heating Oil
#        Sample Oils:
#        - Diesel
#        - Heating Oil
#        - No. 2 Distillate
#        Density Criteria:
#        - 30 <= API < 39
#        Kinematic Viscosity Criteria:
#        - 2.5 < v <= 4.0 cSt @ 38 degrees Celcius
#     '''
#     return (is_refined(oil) and
#             api_min(oil, 30.0) and
#             api_max(oil, 39.0) and
#             is_within_viscosity_range(sample, kvis_min=2.5, kvis_max=4.0))

# def is_refined_ifo(oil, sample):
#     '''
#        Category Name:
#        - Intermediate Fuel Oil
#        Sample Oils:
#        - IFO 180
#        - Fuel Oil #4
#        - Marine Diesel
#        Density Criteria:
#        - 15 <= API < 30
#        Kinematic Viscosity Criteria:
#        - 4.0 < v < 200.0 cSt @ 38 degrees Celcius
#     '''
#     return (is_refined(oil) and
#             api_min(oil, 15.0) and
#             api_max(oil, 30.0) and
#             is_within_viscosity_range(sample, kvis_min=4.0, kvis_max=200.0))

# def is_refined_fuel_oil_6(oil, sample):
#     '''
#        Category Name:
#        - Fuel Oil #6/Bunker/Heavy Fuel Oil/Group V
#        Sample Oils:
#        - Bunker C
#        - Residual Oil
#        Density Criteria:
#        - API < 15
#        Kinematic Viscosity Criteria:
#        - 200.0 <= v cSt @ 50 degrees Celcius
#     '''
#     return (is_refined(oil) and
#             api_max(oil, 15.0) and
#             is_within_viscosity_range(sample, kvis_min=200.0))

# def is_generic(oil):
#     '''
#         Category Name:
#         - Other->Generic
#         Criteria:
#         - Any oils that have been generically generated.  These are found
#           in the OilLibTest data file.  Basically these oils have a name
#           that is prefixed with '*GENERIC'.
#     '''
#     try:
#         ret = oil['metadata'].get('name', None)
#     except AttributeError:
#         ret = oil.metadata.name

#     if ret is not None:
#         return ret.startswith('*GENERIC')
#     else:
#         return None

# def is_within_viscosity_range(oil_sample, kvis_min=None, kvis_max=None):
#     category_temp = 273.15 + 38

#     viscosity = oil_sample.kvis_at_temp(category_temp)

#     if viscosity is None:
#         return False

#     if kvis_min is not None and kvis_max is not None:
#         return (viscosity > kvis_min) and (viscosity <= kvis_max)
#     elif kvis_min is not None:
#         return viscosity > kvis_min
#     elif kvis_max is not None:
#         return viscosity <= kvis_max
#     else:
#         return True

# Data from: https://www.sigma2c.com/16api_gravity.html

#                Specific Gravity  API
# Naphtha light       0.66-0.70
# Naphtha medium      0.70-0.75
# Naphtha heavy       0.75-0.80
# Crude oil           0.80-0.97
# Aviation gasoline   0.70-0.78....47.6--70.6
# Kerosene            0.71-0.79    47.6--67.8
# Gasoline            0.68-0.74    59.7--76.6
# Gas oil             0.78-0.86    33.0--49.9
# Diesel oil          0.82-0.90    25.7--41.1
# Lubricating oil     0.82-0.92
# Fuel oil            0.92-0.99
# Asphaslitc bitumen  1.00-1.10
