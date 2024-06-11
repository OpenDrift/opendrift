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
# Copyright 2020, Manuel Aghito, MET Norway

"""
ChemicalDrift is an OpenDrift module for drift and fate of chemicals.
The module is under development within the scope of the Horizon2020 project EMERGE
Manuel Aghito. Norwegian Meteorological Institute. 2021.
The initial version is based on Radionuclides module by Magne Simonsen
"""

import numpy as np
import logging; logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED
import pyproj
from datetime import datetime

# Defining the Chemical element properties
class Chemical(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for chemicals
    """

    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.}),
        #('neutral_buoyancy_salinity', {'dtype': np.float32,
        #                               'units': '[]',
        #                               'default': 31.25}),  # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 2650.}),  # Mineral particles
        ('specie', {'dtype': np.int32,
                    'units': '',
                    'default': 0}),
        ('mass', {'dtype': np.float32,
                      'units': 'ug',
                      'seed': True,
                      'default': 1e3}),
        ('mass_degraded', {'dtype': np.float32,
                             'units': 'ug',
                             'seed': True,
                             'default': 0}),
        ('mass_degraded_water', {'dtype': np.float32,
                             'units': 'ug',
                             'seed': True,
                             'default': 0}),
        ('mass_degraded_sediment', {'dtype': np.float32,
                             'units': 'ug',
                             'seed': True,
                             'default': 0}),
        ('mass_volatilized', {'dtype': np.float32,
                             'units': 'ug',
                             'seed': True,
                             'default': 0})

        ])


class ChemicalDrift(OceanDrift):
    """Chemical particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Chemical functionality include interactions with solid matter
        (particles and sediments) through transformation processes, implemented
        with stochastic approach for dynamic partitioning.

        Under construction.
    """

    ElementType = Chemical

    required_variables = {
        'x_sea_water_velocity': {'fallback': None},
        'y_sea_water_velocity': {'fallback': None},
        'sea_surface_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'ocean_vertical_diffusivity': {'fallback': 0.0001, 'profiles': True},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'upward_sea_water_velocity': {'fallback': 0},
        #'conc3': {'fallback': 1.e-3},
        'spm': {'fallback': 1},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'active_sediment_layer_thickness': {'fallback': 0.03}, # TODO - currently not used, redundant with 'chemical:sediment:mixing_depth'
        'doc': {'fallback': 0.0},
        # Variables for dissociation
        'sea_water_ph_reported_on_total_scale':{'fallback': 8.1, 'profiles': True}, # water_pH from CMENS with standard name #
        'pH_sediment':{'fallback': 6.9, 'profiles': False}, # supplied by the user, with pH_sediment as standard name #
        }


    def specie_num2name(self,num):
        return self.name_species[num]

    def specie_name2num(self,name):
        num = self.name_species.index(name)
        return num

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(ChemicalDrift, self).__init__(*args, **kwargs)

        # TODO: descriptions and units must be added in config setting below
        self._add_config({
            'chemical:transfer_setup': {'type': 'enum',
                'enum': ['Sandnesfj_Al','metals', '137Cs_rev', 'custom', 'organics'], 'default': 'custom',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:dynamic_partitioning': {'type': 'bool', 'default': True,
                'level': CONFIG_LEVEL_BASIC, 'description': 'Toggle dynamic partitioning'},
            'chemical:slowly_fraction': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:irreversible_fraction': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:dissolved_diameter': {'type': 'float', 'default': 0,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:particle_diameter': {'type': 'float', 'default': 5e-6,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:particle_concentration_half_depth': {'type': 'float', 'default': 20,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:doc_concentration_half_depth': {'type': 'float', 'default': 1000, # TODO: check better
                'min': 0, 'max': 1000, 'units': 'm',                                     # Vertical conc drops more slowly slower than for SPM
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},                # example: 10.3389/fmars.2017.00436. lower limit around 40 umol/L
            'chemical:particle_diameter_uncertainty': {'type': 'float', 'default': 1e-7,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'seed:LMM_fraction': {'type': 'float','default': .1,
                'min': 0, 'max': 1, 'units': '',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'seed:particle_fraction': {'type': 'float','default': 0.9,
                'min': 0, 'max': 1, 'units': '',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            # Species
            'chemical:species:LMM': {'type': 'bool', 'default': True,
                'level': CONFIG_LEVEL_BASIC, 'description': 'Toggle LMM species'},
            'chemical:species:LMMcation': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:LMManion': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:Colloid': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:Humic_colloid': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:species:Polymer': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:species:Particle_reversible': {'type': 'bool', 'default': True,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:Particle_slowly_reversible': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:Particle_irreversible': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:species:Sediment_reversible': {'type': 'bool', 'default': True,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:Sediment_slowly_reversible': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'chemical:species:Sediment_irreversible': {'type': 'bool', 'default': False,
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            # Transformations
            'chemical:transformations:Kd': {'type': 'float', 'default': 2.0,
                'min': 0, 'max': 1e9, 'units': 'm3/kg',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:S0': {'type': 'float', 'default': 0.0,
                'min': 0, 'max': 100, 'units': 'PSU',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'parameter controlling salinity dependency of Kd for metals'},
            'chemical:transformations:Dc': {'type': 'float', 'default': 1.16e-5,                # Simonsen 2019
                'min': 0, 'max': 1e6, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:slow_coeff': {'type': 'float', 'default': 0, #1.2e-7,         # Simonsen 2019
                'min': 0, 'max': 1e6, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:volatilization': {'type': 'bool', 'default': False,
                'description': 'Chemical is evaporated.',
                'level': CONFIG_LEVEL_BASIC},
            'chemical:transformations:degradation': {'type': 'bool', 'default': False,
                'description': 'Chemical mass is degraded.',
                'level': CONFIG_LEVEL_BASIC},
            'chemical:transformations:degradation_mode': {'type': 'enum',
                'enum': ['OverallRateConstants'], 'default': 'OverallRateConstants',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            # sorption/desorption
            'chemical:transformations:dissociation': {'type': 'enum',
                'enum': ['nondiss','acid', 'base', 'amphoter'], 'default': 'nondiss',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:LogKOW': {'type': 'float', 'default': 3.361,          # Naphthalene
                'min': -3, 'max': 10, 'units': 'Log L/Kg',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:TrefKOW': {'type': 'float', 'default': 25.,           # Naphthalene
                'min': -3, 'max': 30, 'units': 'C',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:DeltaH_KOC_Sed': {'type': 'float', 'default': -21036., # Naphthalene
                'min': -100000., 'max': 100000., 'units': 'J/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:DeltaH_KOC_DOM': {'type': 'float', 'default': -25900., # Naphthalene
                'min': -100000., 'max': 100000., 'units': 'J/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:Setchenow': {'type': 'float', 'default': 0.2503,      # Naphthalene
                'min': 0, 'max': 1, 'units': 'L/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:pKa_acid': {'type': 'float', 'default': -1,
                'min': -1, 'max': 14, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:pKa_base': {'type': 'float', 'default': -1,
                'min': -1, 'max': 14, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:KOC_DOM': {'type': 'float', 'default': -1,
                'min': -1, 'max': 10000000000, 'units': 'L/KgOC',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:KOC_sed': {'type': 'float', 'default': -1,
                'min': -1, 'max': 10000000000, 'units': 'L/KgOC',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:fOC_SPM': {'type': 'float', 'default': 0.05,
                'min': 0.01, 'max': 0.1, 'units': 'gOC/g',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:fOC_sed': {'type': 'float', 'default': 0.05,
                'min': 0.01, 'max': 0.1, 'units': 'gOC/g',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:transformations:aggregation_rate': {'type': 'float', 'default': 0,
                'min': 0, 'max': 1, 'units': 's-1',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            # Degradation in water column
            'chemical:transformations:t12_W_tot': {'type': 'float', 'default': 224.08,      # Naphthalene
                'min': 1, 'max': None, 'units': 'hours',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'half life in water, total'},
            'chemical:transformations:Tref_kWt': {'type': 'float', 'default': 25.,          # Naphthalene
                'min': -3, 'max': 30, 'units': 'C',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:DeltaH_kWt': {'type': 'float', 'default': 50000.,     # generic
                'min': -100000., 'max': 100000., 'units': 'J/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            # Degradation in sediment layer
            'chemical:transformations:t12_S_tot': {'type': 'float', 'default': 5012.4,      # Naphthalene
                'min': 1, 'max': None, 'units': 'hours',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'half life in sediments, total'},
            'chemical:transformations:Tref_kSt': {'type': 'float', 'default': 25.,          # Naphthalene
                'min': -3, 'max': 30, 'units': 'C',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            'chemical:transformations:DeltaH_kSt': {'type': 'float', 'default': 50000.,     # generic
                'min': -100000., 'max': 100000., 'units': 'J/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            # Volatilization
            'chemical:transformations:MolWt': {'type': 'float', 'default': 128.1705,         # Naphthalene
                'min': 50, 'max': 1000, 'units': 'amu',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'molecular weight'},
            'chemical:transformations:Henry': {'type': 'float', 'default': 4.551e-4,        # Napththalene
                'min': None, 'max': None, 'units': 'atm m3 mol-1',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Henry constant'},
            # vapour pressure
            'chemical:transformations:Vpress': {'type': 'float', 'default': 11.2,           # Naphthalene
                'min': None, 'max': None, 'units': 'Pa',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Vapour pressure'},
            'chemical:transformations:Tref_Vpress': {'type': 'float', 'default': 25.,        # Naphthalene
                'min': None, 'max': None, 'units': 'C',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Vapour pressure ref temp'},
            'chemical:transformations:DeltaH_Vpress': {'type': 'float', 'default': 55925.,   # Naphthalene
                'min': -100000., 'max': 115000., 'units': 'J/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Enthalpy of volatilization'},
            # solubility
            'chemical:transformations:Solub': {'type': 'float', 'default': 31.4,            # Naphthalene
                'min': None, 'max': None, 'units': 'g/m3',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Solubility'},
            'chemical:transformations:Tref_Solub': {'type': 'float', 'default': 25.,         # Naphthalene
                'min': None, 'max': None, 'units': 'C',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Solubility ref temp'},
            'chemical:transformations:DeltaH_Solub': {'type': 'float', 'default': 25300.,    # Naphthalene
                'min': -100000., 'max': 100000., 'units': 'J/mol',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Enthalpy of solubilization'},
            # Sedimentation/Resuspension
            'chemical:sediment:mixing_depth': {'type': 'float', 'default': 0.03,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:density': {'type': 'float', 'default': 2600,
                'min': 0, 'max': 10000, 'units': 'kg/m3',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:effective_fraction': {'type': 'float', 'default': 0.9,
                'min': 0, 'max': 1, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:corr_factor': {'type': 'float', 'default': 0.1,
                'min': 0, 'max': 10, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:porosity': {'type': 'float', 'default': 0.6,
                'min': 0, 'max': 1, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:layer_thickness': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:desorption_depth': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:desorption_depth_uncert': {'type': 'float', 'default': .5,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:resuspension_depth': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:resuspension_depth_uncert': {'type': 'float', 'default': .5,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:resuspension_critvel': {'type': 'float', 'default': .01,
                'min': 0, 'max': 1, 'units': 'm/s',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:burial_rate': {'type': 'float', 'default': .00003,   # MacKay
                'min': 0, 'max': 10, 'units': 'm/year',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'chemical:sediment:buried_leaking_rate': {'type': 'float', 'default': 0,
                'min': 0, 'max': 10, 'units': 's-1',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            #
            'chemical:compound': {'type': 'enum',
                'enum': ['Naphthalene','Phenanthrene','Fluoranthene',
                         'Benzo-a-anthracene','Benzo-a-pyrene','Dibenzo-ah-anthracene',
                         'C1-Naphthalene','Acenaphthene','Acenaphthylene','Fluorene',
                         'Dibenzothiophene','C2-Naphthalene','Anthracene','C3-Naphthalene','C1-Dibenzothiophene',
                         'Pyrene','C1-Phenanthrene','C2-Dibenzothiophene',
                         'C2-Phenanthrene','Benzo-b-fluoranthene','Chrysene',
                         'C3-Dibenzothiophene','C3-Phenanthrene',
                         'Benzo-k-fluoranthene','Benzo-ghi-perylene','Indeno-123cd-pyrene',
                         'Copper','Cadmium','Chromium','Lead','Vanadium','Zinc','Nickel',None],
                'default': None,
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': ''},
            })


    def prepare_run(self):

        logger.info( 'Number of species: {}'.format(self.nspecies) )
        for i,sp in enumerate(self.name_species):
            logger.info( '{:>3} {}'.format( i, sp ) )


        logger.info( 'transfer setup: %s' % self.get_config('chemical:transfer_setup'))

        logger.info('nspecies: %s' % self.nspecies)
        logger.info('Transfer rates:\n %s' % self.transfer_rates)

        self.SPM_vertical_levels_given = False
        for key, value in self.env.readers.items():
            if 'spm' in value.variables:
                if (hasattr(value,'sigma') or hasattr(value,'z') ):
                    self.SPM_vertical_levels_given = True

        self.DOC_vertical_levels_given = False
        for key, value in self.env.readers.items():
            if 'doc' in value.variables:
                if (hasattr(value,'sigma') or hasattr(value,'z') ):
                    self.DOC_vertical_levels_given = True

        super(ChemicalDrift, self).prepare_run()

    def init_species(self):
        # Initialize specie types
        if self.get_config('chemical:transfer_setup')=='metals':
            self.set_config('chemical:species:LMM',True)
            self.set_config('chemical:species:Particle_reversible', True)
            self.set_config('chemical:species:Particle_slowly_reversible', True)
            self.set_config('chemical:species:Sediment_reversible', True)
            self.set_config('chemical:species:Sediment_slowly_reversible', True)
        elif self.get_config('chemical:transfer_setup')=='137Cs_rev':
            self.set_config('chemical:species:LMM',True)
            self.set_config('chemical:species:Particle_reversible', True)
            self.set_config('chemical:species:Sediment_reversible', True)
        elif self.get_config('chemical:transfer_setup')=='Sandnesfj_Al':
            self.set_config('chemical:species:LMM', False)
            self.set_config('chemical:species:LMMcation', True)
            self.set_config('chemical:species:LMManion', True)
            self.set_config('chemical:species:Humic_colloid', True)
            self.set_config('chemical:species:Polymer', True)
            self.set_config('chemical:species:Particle_reversible', True)
            self.set_config('chemical:species:Sediment_reversible', True)
        elif self.get_config('chemical:transfer_setup')=='organics':
            self.set_config('chemical:species:LMM',True)
            self.set_config('chemical:species:Particle_reversible', True)
            self.set_config('chemical:species:Particle_slowly_reversible', False)
            self.set_config('chemical:species:Sediment_reversible', True)
            self.set_config('chemical:species:Sediment_slowly_reversible', True)
            self.set_config('chemical:species:Humic_colloid', True)
        elif self.get_config('chemical:transfer_setup')=='custom':
            # Do nothing, species must be set manually
            pass
        else:
            logger.error('No valid transfer_setup {}'.format(self.get_config('chemical:transfer_setup')))


        self.name_species=[]
        if self.get_config('chemical:species:LMM'):
            self.name_species.append('LMM')
        if self.get_config('chemical:species:LMMcation'):
            self.name_species.append('LMMcation')
        if self.get_config('chemical:species:LMManion'):
            self.name_species.append('LMManion')
        if self.get_config('chemical:species:Colloid'):
            self.name_species.append('Colloid')
        if self.get_config('chemical:species:Humic_colloid'):
            self.name_species.append('Humic colloid')
        if self.get_config('chemical:species:Polymer'):
            self.name_species.append('Polymer')
        if self.get_config('chemical:species:Particle_reversible'):
            self.name_species.append('Particle reversible')
        if self.get_config('chemical:species:Particle_slowly_reversible'):
            self.name_species.append('Particle slowly reversible')
        if self.get_config('chemical:species:Particle_irreversible'):
            self.name_species.append('Particle irreversible')
        if self.get_config('chemical:species:Sediment_reversible'):
            self.name_species.append('Sediment reversible')
        if self.get_config('chemical:species:Sediment_slowly_reversible'):
            self.name_species.append('Sediment slowly reversible')
        if self.get_config('chemical:species:Sediment_irreversible'):
            self.name_species.append('Sediment irreversible')


        if self.get_config('chemical:species:Sediment_slowly_reversible') and \
                    self.get_config('chemical:species:Particle_slowly_reversible'):
            self.set_config('chemical:slowly_fraction', True)
        if self.get_config('chemical:species:Sediment_irreversible') and \
                    self.get_config('chemical:species:Particle_irreversible'):
            self.set_config('chemical:irreversible_fraction', True)


        self.nspecies      = len(self.name_species)
#         logger.info( 'Number of species: {}'.format(self.nspecies) )
#         for i,sp in enumerate(self.name_species):
#             logger.info( '{:>3} {}'.format( i, sp ) )



    def seed_elements(self, *args, **kwargs):

        if hasattr(self,'name_species') == False:
            self.init_species()
            self.init_transfer_rates()



        if 'number' in kwargs:
            num_elements = kwargs['number']
        else:
            num_elements = self.get_config('seed:number')


        if 'specie' in kwargs:
            print('num_elements', num_elements)
            try:
                print('len specie:',len(kwargs['specie']))
            except:
                print('specie:',kwargs['specie'])

            init_specie = np.ones(num_elements,dtype=int)
            init_specie[:] = kwargs['specie']


        else:

            # Set initial partitioning
            if 'particle_fraction' in kwargs:
                particle_frac = kwargs['particle_fraction']
            else:
                particle_frac = self.get_config('seed:particle_fraction')

            if 'LMM_fraction' in kwargs:
                lmm_frac = kwargs['LMM_fraction']
            else:
                lmm_frac = self.get_config('seed:LMM_fraction')

            if not lmm_frac + particle_frac == 1.:
                logger.error('Fraction does not sum up to 1: %s' % str(lmm_frac+particle_frac) )
                logger.error('LMM fraction: %s ' % str(lmm_frac))
                logger.error( 'Particle fraction %s '% str(particle_frac) )
                raise ValueError('Illegal specie fraction combination : ' + str(lmm_frac) + ' '+ str(particle_frac) )

            init_specie = np.ones(num_elements, int)

            dissolved=np.random.rand(num_elements)<lmm_frac
            if self.get_config('chemical:transfer_setup')=='Sandnesfj_Al':
                init_specie[dissolved]=self.num_lmmcation
            else:
                init_specie[dissolved]=self.num_lmm
            init_specie[~dissolved]=self.num_prev


            kwargs['specie'] = init_specie


        logger.debug('Initial partitioning:')
        for i,sp in enumerate(self.name_species):
            logger.debug( '{:>9} {:>3} {:24} '.format(  np.sum(init_specie==i), i, sp ) )

        # Set initial particle size
        if 'diameter' in kwargs:
            diameter = kwargs['diameter']
        else:
            diameter = self.get_config('chemical:particle_diameter')

        std = self.get_config('chemical:particle_diameter_uncertainty')

        init_diam = np.zeros(num_elements,float)
        init_diam[init_specie==self.num_prev] = diameter + np.random.normal(0, std, sum(init_specie==self.num_prev))
        kwargs['diameter'] = init_diam



        super(ChemicalDrift, self).seed_elements(*args, **kwargs)

    def tempcorr(self,mode,DeltaH,T_C,Tref_C):
        ''' Temperature correction using Arrhenius or Q10 method
        '''
        if mode == 'Arrhenius':
            R = 8.3145 # J/(mol*K)
            T_K = T_C + 273.15
            Tref_K = Tref_C + 273.15
            corr = np.e**(-(DeltaH/R)*(1/T_K - 1/Tref_K))
        elif mode =='Q10':
            corr = 2**((T_C - Tref_C)/10)
        return corr

    def salinitycorr(self,Setschenow,Temperature,Salinity):
        ''' Salinity correction
        '''
        # Setschenow constant for the given chemical (L/mol)
        # Salinity   (PSU =g/Kg)
        # Temperature (Celsius)

        MWsalt = 68.35 # average mass of sea water salt (g/mol) Schwarzenbach Gschwend Imboden Environmental Organic Chemistry

        Dens_sw = self.sea_water_density(T=Temperature, S=Salinity)*1e-3 # (Kg/L)

        # ConcSalt= (Salinitypsu/MWsalt)∙Dens_sw
        #         = (     g/Kg    /    g/mol  )∙  Kg/L
        #         = mol/Kg ∙ Kg/L = mol/L

        ConcSalt = (Salinity/MWsalt)*Dens_sw

        # Log(Kd_fin)=(Setschenow ∙ ConcSalt)+Log(Kd_T)
        # Kd_fin = 10^(Setschenow ∙ ConcSalt) * Kd_T

        corr = 10**(Setschenow*ConcSalt)

        return corr

### Functions to update partitioning coefficients

    def calc_KOC_sedcorr(self, KOC_sed_initial, KOC_sed_n, pKa_acid, pKa_base, KOW, pH_sed, diss):
        ''' Calculate correction of KOC due to pH of sediments
        '''
        # Estimate KOC for dissociated forms from KOW
        KOC_sed_diss_acid = (10**(0.11*np.log10(KOW)+1.54)) # KOC for dissociated acid species (L/kg_OC), from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202.
        KOC_sed_diss_base = 10**(pKa_acid**(0.65*((KOW/(KOW+1))**0.14))) # KOC for ionized form of base species (L/kg_OC) # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202
        # Updated values of KOC to calculate correction factor
        KOC_sed_updated = np.empty_like(pH_sed)
        KOC_sedcorr= np.empty_like(pH_sed)

        for i in (range(len(pH_sed))):
            if diss=='acid':

                Phi_n_sed    = 1/(1 + 10**(pH_sed[i]-pKa_acid))
                Phi_diss_sed = 1-Phi_n_sed
                KOC_sed_updated[i] = (KOC_sed_n * Phi_n_sed) + (Phi_diss_sed * KOC_sed_diss_acid)
                KOC_sedcorr[i] = KOC_sed_updated[i]/KOC_sed_initial

            elif diss=='base':

                Phi_n_sed    = 1/(1 + 10**(pH_sed[i]-pKa_base))
                Phi_diss_sed = 1-Phi_n_sed
                KOC_sed_updated[i] = (KOC_sed_n * Phi_n_sed) + (Phi_diss_sed * KOC_sed_diss_acid)
                KOC_sedcorr[i] = KOC_sed_updated[i]/KOC_sed_initial

            elif diss=='amphoter':

                Phi_n_sed      = 1/(1 + 10**(pH_sed[i]-pKa_acid) + 10**(pKa_base))
                Phi_anion_sed  = Phi_n_sed * 10**(pH_sed[i]-pKa_acid)
                Phi_cation_sed = Phi_n_sed * 10**(pKa_base-pH_sed[i])

                KOC_sed_updated[i] = (KOC_sed_n * Phi_n_sed) + (Phi_anion_sed * KOC_sed_diss_acid) + (Phi_cation_sed * KOC_sed_diss_base)
                KOC_sedcorr[i]=KOC_sed_updated[i]/KOC_sed_initial

            elif diss=='undiss':
                KOC_sedcorr[i]=1

        return KOC_sedcorr

    def calc_KOC_watcorrSPM(self, KOC_SPM_initial, KOC_sed_n, pKa_acid, pKa_base, KOW, pH_water_SPM, diss):
        ''' Calculate correction of KOC due to pH of water for SPM
        '''
        # Estimate KOC for dissociated forms from KOW
        KOC_sed_diss_acid = (10**(0.11*np.log10(KOW)+1.54)) # KOC for dissociated acid species (L/kg_OC), from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202.

        KOC_sed_diss_base = 10**(pKa_acid**(0.65*((KOW/(KOW+1))**0.14))) # KOC for ionized form of base species (L/kg_OC) # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

        KOC_SPM_updated = np.empty_like(pH_water_SPM)
        KOC_SPMcorr = np.empty_like(pH_water_SPM)

        for i in (range(len(pH_water_SPM))):
            if diss=='acid':

                Phi_n_SPM    = 1/(1 + 10**(pH_water_SPM[i]-pKa_acid))
                Phi_diss_SPM = 1-Phi_n_SPM
                KOC_SPM_updated[i] = (KOC_sed_n * Phi_n_SPM) + (Phi_diss_SPM * KOC_sed_diss_acid)

                KOC_SPMcorr[i]=KOC_SPM_updated[i]/KOC_SPM_initial

            elif diss=='base':

                Phi_n_SPM    = 1/(1 + 10**(pH_water_SPM[i]-pKa_base))
                Phi_diss_SPM = 1-Phi_n_SPM
                KOC_SPM_updated[i] = (KOC_sed_n * Phi_n_SPM) + (Phi_diss_SPM * KOC_sed_diss_acid)

                KOC_SPMcorr[i]=KOC_SPM_updated[i]/KOC_SPM_initial

            elif diss=='amphoter':

                Phi_n_SPM      = 1/(1 + 10**(pH_water_SPM[i]-pKa_acid) + 10**(pKa_base))
                Phi_anion_SPM  = Phi_n_SPM * 10**(pH_water_SPM[i]-pKa_acid)
                Phi_cation_SPM = Phi_n_SPM * 10**(pKa_base-pH_water_SPM[i])
                KOC_SPM_updated[i] = (KOC_sed_n * Phi_n_SPM) + (Phi_anion_SPM * KOC_sed_diss_acid) + (Phi_cation_SPM * KOC_sed_diss_base)

                KOC_SPMcorr[i]=KOC_SPM_updated[i]/KOC_SPM_initial

            elif diss=='undiss':
                KOC_SPMcorr[i]=1

        return KOC_SPMcorr

    def calc_KOC_watcorrDOM(self, KOC_DOM_initial, KOC_DOM_n, pKa_acid, pKa_base, KOW, pH_water_DOM, diss):
        ''' Calculate correction of KOC due to pH of water for DOM
        '''

        Phi_n_DOM = np.empty_like(pH_water_DOM)
        Phi_diss_DOM = np.empty_like(pH_water_DOM)
        KOC_DOM_updated = np.empty_like(pH_water_DOM)
        KOC_DOMcorr = np.empty_like(pH_water_DOM)

        for i in (range(len(pH_water_DOM))):
            if diss=='acid':

                Phi_n_DOM    = 1/(1 + 10**(pH_water_DOM[i]-pKa_acid))
                Phi_diss_DOM    = 1-Phi_n_DOM
                KOC_DOM_updated[i] = (0.08 * ((Phi_n_DOM*(KOC_DOM_n)) + ((1 - Phi_diss_DOM)*10**(np.log10(KOW)-3.5))))/0.526 # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

                KOC_DOMcorr[i]=KOC_DOM_updated[i]/KOC_DOM_initial

            elif diss=='base':

                Phi_n_DOM    = 1/(1 + 10**(pH_water_DOM[i]-pKa_base))
                Phi_diss_DOM    = 1-Phi_n_DOM
                KOC_DOM_updated[i] = (0.08 * ((Phi_n_DOM*(KOC_DOM_n)) + ((1 - Phi_diss_DOM)*10**(np.log10(KOW)-3.5))))/0.526 # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

                KOC_DOMcorr[i]=KOC_DOM_updated[i]/KOC_DOM_initial

            elif diss=='amphoter':

                Phi_n_DOM      = 1/(1 + 10**(pH_water_DOM[i]-pKa_acid) + 10**(pKa_base))
                Phi_diss_DOM   = 1-Phi_n_DOM
                KOC_DOM_updated[i] = (0.08 * ((Phi_n_DOM*(KOC_DOM_n)) + ((1 - Phi_diss_DOM)*10**(np.log10(KOW)-3.5))))/0.526 # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202
                KOC_DOMcorr[i]=KOC_DOM_updated[i]/KOC_DOM_initial

            elif diss=='undiss':
                    KOC_DOMcorr[i]=1

        return KOC_DOMcorr


    def init_transfer_rates(self):
        ''' Initialization of background values in the transfer rates 2D array.
        '''

        transfer_setup=self.get_config('chemical:transfer_setup')

#        logger.info( 'transfer setup: %s' % transfer_setup)


        self.transfer_rates = np.zeros([self.nspecies,self.nspecies])
        self.ntransformations = np.zeros([self.nspecies,self.nspecies])

        if transfer_setup == 'organics':

            self.num_lmm    = self.specie_name2num('LMM')
            self.num_humcol = self.specie_name2num('Humic colloid')
            self.num_prev   = self.specie_name2num('Particle reversible')
            self.num_srev   = self.specie_name2num('Sediment reversible')
            #self.num_psrev  = self.specie_name2num('Particle slowly reversible')
            self.num_ssrev  = self.specie_name2num('Sediment slowly reversible')

            # Values from EMERGE-Aquatox
            Org2C      = 0.526  # kgOC/KgOM
            #Kd         = self.get_config('chemical:transformations:Kd')
            KOW        = 10**self.get_config('chemical:transformations:LogKOW')
            KOWTref    = self.get_config('chemical:transformations:TrefKOW')
            DH_KOC_Sed = self.get_config('chemical:transformations:DeltaH_KOC_Sed')
            DH_KOC_DOM = self.get_config('chemical:transformations:DeltaH_KOC_DOM')
            Setchenow  = self.get_config('chemical:transformations:Setchenow')


            diss       = self.get_config('chemical:transformations:dissociation')
            pKa_acid   = self.get_config('chemical:transformations:pKa_acid')
            if pKa_acid < 0 and diss!='nondiss':
                raise ValueError("pKa_acid must be positive")
                # print("pKa_acid must be positive")
                # UserWarning(("pKa_acid must be positive"))

            else:
                pass

            pKa_base   = self.get_config('chemical:transformations:pKa_base')
            if pKa_base < 0 and diss!='nondiss':
                raise ValueError("pKa_base must be positive")
                # print("pKa_base must be positive")
                # UserWarning(("pKa_base must be positive"))

            else:
                pass

            if diss == 'amphoter' and abs(pKa_acid - pKa_base) < 2:
                raise ValueError("pKa_base and pKa_acid must differ of at least two units")
            else:
                pass

            # Read water pH to calculate dissociation
            # pH_water = self.environment.sea_water_ph_reported_on_total_scale
            pH_water   = 8.1 # 8.1

            pH_sed     = 6.9 # 6.9

            fOC_SPM    = self.get_config('chemical:transformations:fOC_SPM')       # typical values from 0.01 to 0.1 gOC/g
            fOC_sed    = self.get_config('chemical:transformations:fOC_sed')       # typical values from 0.01 to 0.1 gOC/g

            concDOM   = 1.e-3 / Org2C    # concentration of available dissolved organic matter (kg/m3)
                                         # rough initial estimate for coastal waters, doi: 10.1002/lom3.10118
            #concDOM   = 50.e-3     # HIGHER VALUE FOR TESTING!!!!!!!!!!!!

            # Values from Simonsen et al (2019a)
            slow_coeff  = self.get_config('chemical:transformations:slow_coeff')
            concSPM     = 50.e-3                                                # available SPM (kg/m3)
            sed_L       = self.get_config('chemical:sediment:mixing_depth')     # sediment mixing depth (m)
            sed_dens    = self.get_config('chemical:sediment:density')          # default particle density (kg/m3)
            sed_phi     = self.get_config('chemical:sediment:corr_factor')      # sediment correction factor
            sed_poro    = self.get_config('chemical:sediment:porosity')         # sediment porosity
            sed_H       = self.get_config('chemical:sediment:layer_thickness')  # thickness of seabed interaction layer (m)
            sed_burial  = self.get_config('chemical:sediment:burial_rate')      # sediment burial rate (m/y)
            sed_leaking_rate = self.get_config( 'chemical:sediment:buried_leaking_rate')

            if diss=='nondiss':
                KOC_DOM = self.get_config('chemical:transformations:KOC_DOM')
                if KOC_DOM < 0:
                    KOC_DOM = 2.88 * KOW**0.67   # (L/KgOC), Park and Clough, 2014

                KOC_sed = self.get_config('chemical:transformations:KOC_sed')
                if KOC_sed < 0:
                    KOC_sed = 2.62 * KOW**0.82   # (L/KgOC), Park and Clough, 2014 (334)/Org2C
                    #KOC_Sed    = 1.26 * kOW**0.81   # (L/KgOC),Ragas et al., 2019

                KOC_SPM = KOC_sed

            else:
                if diss=='acid':
                    # Dissociation in water
                    Phi_n_water    = 1/(1 + 10**(pH_water-pKa_acid))
                    Phi_diss_water = 1-Phi_n_water

                    KOC_sed_n = self.get_config('chemical:transformations:KOC_sed')

                    if KOC_sed_n < 0:
                        KOC_sed_n    = 2.62 * KOW**0.82   # (L/KgOC), Park and Clough, 2014 (334)/Org2C TO DO Add if choice between input and estimation
                    else:
                        pass

                    KOC_sed_diss_acid = (10**(0.11*np.log10(KOW)+1.54)) # KOC for dissociated acid species (L/kg_OC), from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202.

                    KOC_SPM = (KOC_sed_n * Phi_n_water) + (Phi_diss_water * KOC_sed_diss_acid)

                    KOC_DOM_n = self.get_config('chemical:transformations:KOC_DOM')

                    if KOC_DOM_n <0:
                        KOC_DOM_n   = 2.88 * KOW**0.67   # (L/KgOC), Park and Clough, 2014 TO DO Add if choice between input and estimation
                    else:
                        pass

                    KOC_DOM = (0.08 * ((Phi_n_water*(KOC_DOM_n)) + ((1 - Phi_diss_water)*10**(np.log10(KOW)-3.5))))/0.526 # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

                    # Dissociation in sediments
                    Phi_n_sed    = 1/(1 + 10**(pH_sed-pKa_acid))
                    Phi_diss_sed = 1-Phi_n_sed
                    KOC_sed = (KOC_sed_n * Phi_n_sed) + (Phi_diss_sed * KOC_sed_diss_acid)

                elif diss=='base':
                    # Dissociation in water
                    Phi_n_water    = 1/(1 + 10**(pH_water-pKa_base))
                    Phi_diss_water = 1-Phi_n_water

                    KOC_sed_n = self.get_config('chemical:transformations:KOC_sed')
                    if KOC_sed_n <0:
                        KOC_sed_n   = 2.62 * KOW**0.82   # (L/KgOC), Park and Clough, 2014 (334)/Org2C TO DO Add if choice between input and estimation
                    else:
                        pass

                    KOC_sed_diss_base = 10**(pKa_acid**(0.65*((KOW/(KOW+1))**0.14))) # KOC for ionized form of base species (L/kg_OC) # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

                    KOC_SPM = (KOC_sed_n * Phi_n_water) + (Phi_diss_water * KOC_sed_diss_base)

                    KOC_DOM_n = self.get_config('chemical:transformations:KOC_DOM')
                    if   KOC_DOM_n <0:
                        KOC_DOM_n   = 2.88 * KOW**0.67   # (L/KgOC), Park and Clough, 2014 TO DO Add if choice between input and estimation
                    else:
                        pass

                    KOC_DOM = (0.08 * ((Phi_n_water*(KOC_DOM_n)) + ((1 - Phi_diss_water)*10**(np.log10(KOW)-3.5))))/0.526 # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

                    # Dissociation in sediments
                    Phi_n_sed    = 1/(1 + 10**(pH_sed-pKa_base))
                    Phi_diss_sed = 1-Phi_n_sed
                    KOC_sed = (KOC_sed_n * Phi_n_sed) + (Phi_diss_sed * KOC_sed_diss_base)

                elif diss=='amphoter':

                    # Dissociation in water # This approach ignores the zwitterionic fraction. 10.1002/etc.115
                    Phi_n_water      = 1/(1 + 10**(pH_water-pKa_acid) + 10**(pKa_base))
                    Phi_anion_water  = Phi_n_water * 10**(pH_water-pKa_acid)
                    Phi_cation_water = Phi_n_water * 10**(pKa_base-pH_water)
                    Phi_diss_water   = 1 - Phi_n_water

                    KOC_sed_n = self.get_config('chemical:transformations:KOC_sed')
                    if KOC_sed_n < 0:
                        KOC_sed_n =  KOC_sed_n    = 2.62 * KOW**0.82   # (L/KgOC), Park and Clough, 2014 (334)/Org2C TO DO Add if choice between input and estimation
                    else:
                        pass

                    KOC_sed_diss_base = 10**(pKa_acid**(0.65*((KOW/(KOW+1))**0.14))) # KOC for ionized form of base species (L/kg_OC) # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202
                    KOC_sed_diss_acid = (10**(0.11*np.log10(KOW)+1.54)) # KOC for dissociated acid species (L/kg_OC), from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202.
                    KOC_SPM = (KOC_sed_n * Phi_n_water) + (Phi_anion_water * KOC_sed_diss_acid) + (Phi_cation_water * KOC_sed_diss_base)

                    KOC_DOM_n = self.get_config('chemical:transformations:KOC_DOM')
                    if KOC_DOM_n <0:
                        KOC_DOM_n   = 2.88 * KOW**0.67   # (L/KgOC), Park and Clough, 2014 TO DO Add if choice between input and estimation
                    else:
                        pass

                    KOC_DOM = (0.08 * ((Phi_n_water*(KOC_DOM_n)) + ((1 - Phi_diss_water)*10**(np.log10(KOW)-3.5))))/0.526 # from  http://i-pie.org/wp-content/uploads/2019/12/ePiE_Technical_Manual-Final_Version_20191202

                    # Dissociation in sediments
                    Phi_n_sed      = 1/(1 + 10**(pH_sed-pKa_acid) + 10**(pKa_base))
                    Phi_anion_sed  = Phi_n_sed * 10**(pH_sed-pKa_acid)
                    Phi_cation_sed = Phi_n_sed * 10**(pKa_base-pH_sed)

                    KOC_sed = (KOC_sed_n * Phi_n_sed) + (Phi_anion_sed * KOC_sed_diss_acid) + (Phi_cation_sed * KOC_sed_diss_base)

            logger.debug('Partitioning coefficients (Tref,freshwater)')
            logger.debug('KOC_sed: %s L/KgOC' % KOC_sed)
            logger.debug('KOC_SPM: %s L/KgOC' % KOC_SPM)
            logger.debug('KOC_DOM: %s L/KgOC' % KOC_DOM)

            #KOM_sed = KOC_sed * Org2C #  L/KgOC * KgOC/KgOM = L/KgOM
            #KOM_SPM = KOC_sed * Org2C #  L/KgOC * KgOC/KgOM = L/KgOM
            #KOM_DOM = KOC_DOM * Org2C #  L/KgOC * KgOC/KgOM = L/KgOM

            # to be calculated separately for sed, SPM, dom (different KOC, pH, fOC)
            self.Kd_sed = Kd_sed = KOC_sed * fOC_sed    # L/KgOC * KgOC/KG = L/Kg
            self.Kd_SPM = Kd_SPM = KOC_SPM * fOC_SPM    # L/KgOC * KgOC/KG = L/Kg
            self.Kd_DOM = Kd_DOM = KOC_DOM * Org2C      # L/KgOC * KgOC/KgOM * 1KgOM/Kg = L/Kg (=KOM_DOM)
            # TODO Use setconfig() to store these?

            logger.debug('Kd_sed: %s L/Kg' % Kd_sed)
            logger.debug('Kd_SPM: %s L/Kg' % Kd_SPM)
            logger.debug('Kd_DOM: %s L/Kg' % Kd_DOM)

            # From Karickhoff and Morris 1985
            k_ads = 33.3 / (60*60) # L/(Kg*s) = 33 L/(kgOM*h)

            k_des_sed = k_ads / Kd_sed # 1/s
            k_des_SPM = k_ads / Kd_SPM # 1/s
            k_des_DOM = k_ads / Kd_DOM # 1/s

            # Default corrections, assuming temperature 25 salinity 35
            TcorrSed = self.tempcorr("Arrhenius",DH_KOC_Sed,25,KOWTref)
            TcorrDOM = self.tempcorr("Arrhenius",DH_KOC_DOM,25,KOWTref)
            Scorr    = self.salinitycorr(Setchenow,KOWTref,35)

            concSPM = concSPM * 1e-3 # (Kg/L)
            concDOM = concDOM * 1e-3 # (Kg/L)

            self.k_ads = k_ads
            self.k21_0 = k_des_DOM
            self.k31_0 = k_des_SPM
            self.k41_0 = k_des_sed * sed_phi
            # TODO Use setconfig() to store these?

            self.transfer_rates[self.num_lmm,self.num_humcol] = k_ads * concDOM             # k12
            self.transfer_rates[self.num_humcol,self.num_lmm] = k_des_DOM / TcorrDOM / Scorr# k21

            self.transfer_rates[self.num_lmm,self.num_prev] = k_ads * concSPM               # k13
            self.transfer_rates[self.num_prev,self.num_lmm] = k_des_SPM / TcorrSed / Scorr  # k31

            self.transfer_rates[self.num_lmm,self.num_srev] = \
                k_ads * sed_L * sed_dens * (1.-sed_poro) * sed_phi / sed_H                  # k14
                # TODO CHECK DIMENSIONS!!!!! L-m3 !!!!

            self.transfer_rates[self.num_srev,self.num_lmm] = \
                k_des_sed * sed_phi / TcorrSed / Scorr                                      # k41

            #self.transfer_rates[self.num_srev,self.num_ssrev] = slow_coeff                  # k46
            #self.transfer_rates[self.num_ssrev,self.num_srev] = slow_coeff*.1               # k64

            # Using slowly reversible specie for burial - TODO buried sediment should be a new specie
            self.transfer_rates[self.num_srev,self.num_ssrev] = sed_burial / sed_L / 31556926 # k46 (m/y) / m / (s/y) = s-1
            self.transfer_rates[self.num_ssrev,self.num_srev] = sed_leaking_rate                # k64


            self.transfer_rates[self.num_humcol,self.num_prev] = self.get_config('chemical:transformations:aggregation_rate')
            self.transfer_rates[self.num_prev,self.num_humcol] = 0          # TODO check if valid for organics

        elif transfer_setup == 'metals':                                # renamed from radionuclides Bokna_137Cs

            self.num_lmm    = self.specie_name2num('LMM')
            self.num_prev   = self.specie_name2num('Particle reversible')
            self.num_srev   = self.specie_name2num('Sediment reversible')
            self.num_psrev  = self.specie_name2num('Particle slowly reversible')
            self.num_ssrev  = self.specie_name2num('Sediment slowly reversible')


            # Values from Simonsen et al (2019a)
            Kd         = self.get_config('chemical:transformations:Kd') # (m3/Kg)
            Dc         = self.get_config('chemical:transformations:Dc') # (1/s)
            slow_coeff = self.get_config('chemical:transformations:slow_coeff')
            concSPM    = 1.e-3   # concentration of available suspended particulate matter (kg/m3)
            sed_L = self.get_config('chemical:sediment:mixing_depth')     # sediment mixing depth (m)
            sed_dens =  self.get_config('chemical:sediment:density') # default particle density (kg/m3)
            sed_f           =  self.get_config('chemical:sediment:effective_fraction')      # fraction of effective sorbents
            sed_phi         =  self.get_config('chemical:sediment:corr_factor')      # sediment correction factor
            sed_poro        =  self.get_config('chemical:sediment:porosity')      # sediment porosity
            sed_H =  self.get_config('chemical:sediment:layer_thickness')      # thickness of seabed interaction layer (m)

            #self.k_ads = Dc * Kd * 1e3 # L/(Kg*s)
            self.transfer_rates[self.num_lmm,self.num_prev] = Dc * Kd * concSPM
            self.transfer_rates[self.num_prev,self.num_lmm] = Dc
            self.transfer_rates[self.num_lmm,self.num_srev] = \
                Dc * Kd * sed_L * sed_dens * (1.-sed_poro) * sed_f * sed_phi / sed_H
            self.transfer_rates[self.num_srev,self.num_lmm] = Dc * sed_phi
            self.transfer_rates[self.num_srev,self.num_ssrev] = slow_coeff
            self.transfer_rates[self.num_prev,self.num_psrev] = slow_coeff
            self.transfer_rates[self.num_ssrev,self.num_srev] = slow_coeff*.1
            self.transfer_rates[self.num_psrev,self.num_prev] = slow_coeff*.1


        elif transfer_setup == '137Cs_rev':

            self.num_lmm    = self.specie_name2num('LMM')
            self.num_prev   = self.specie_name2num('Particle reversible')
            self.num_srev   = self.specie_name2num('Sediment reversible')


            # Simpler version of Values from Simonsen et al (2019a)
            # Only consider the reversible fraction
            Kd         = self.get_config('chemical:transformations:Kd')
            Dc         = self.get_config('chemical:transformations:Dc')
            concSPM    = 1.e-3   # concentration of available suspended particulate matter (kg/m3)
            sed_L           = self.get_config('chemical:sediment:mixing_depth')     # sediment mixing depth (m)
            sed_dens        = self.get_config('chemical:sediment:density') # default particle density (kg/m3)
            sed_f           = self.get_config('chemical:sediment:effective_fraction')      # fraction of effective sorbents
            sed_phi         = self.get_config('chemical:sediment:corr_factor')      # sediment correction factor
            sed_poro        = self.get_config('chemical:sediment:porosity')      # sediment porosity
            sed_H =  self.get_config('chemical:sediment:layer_thickness')      # thickness of seabed interaction layer (m)

            self.transfer_rates[self.num_lmm,self.num_prev] = Dc * Kd * concSPM
            self.transfer_rates[self.num_prev,self.num_lmm] = Dc
            self.transfer_rates[self.num_lmm,self.num_srev] = \
                Dc * Kd * sed_L * sed_dens * (1.-sed_poro) * sed_f * sed_phi / sed_H
            self.transfer_rates[self.num_srev,self.num_lmm] = Dc * sed_phi

        elif transfer_setup=='custom':
        # Set of custom values for testing/development

            self.num_lmm   = self.specie_name2num('LMM')
            if self.get_config('chemical:species:Colloid'):
                self.num_col = self.specie_name2num('Colloid')
            if self.get_config('chemical:species:Particle_reversible'):
                self.num_prev  = self.specie_name2num('Particle reversible')
            if self.get_config('chemical:species:Sediment_reversible'):
                self.num_srev  = self.specie_name2num('Sediment reversible')
            if self.get_config('chemical:slowly_fraction'):
                self.num_psrev  = self.specie_name2num('Particle slowly reversible')
                self.num_ssrev  = self.specie_name2num('Sediment slowly reversible')
            if self.get_config('chemical:irreversible_fraction'):
                self.num_pirrev  = self.specie_name2num('Particle irreversible')
                self.num_sirrev  = self.specie_name2num('Sediment irreversible')


            if self.get_config('chemical:species:Particle_reversible'):
                self.transfer_rates[self.num_lmm,self.num_prev] = 5.e-6 #*0.
                self.transfer_rates[self.num_prev,self.num_lmm] = \
                    self.get_config('chemical:transformations:Dc')
            if self.get_config('chemical:species:Sediment_reversible'):
                self.transfer_rates[self.num_lmm,self.num_srev] = 1.e-5 #*0.
                self.transfer_rates[self.num_srev,self.num_lmm] = \
                    self.get_config('chemical:transformations:Dc') * self.get_config('chemical:sediment:corr_factor')
#                self.transfer_rates[self.num_srev,self.num_lmm] = 5.e-6

            if self.get_config('chemical:slowly_fraction'):
                self.transfer_rates[self.num_prev,self.num_psrev] = 2.e-6
                self.transfer_rates[self.num_srev,self.num_ssrev] = 2.e-6
                self.transfer_rates[self.num_psrev,self.num_prev] = 2.e-7
                self.transfer_rates[self.num_ssrev,self.num_srev] = 2.e-7

        elif transfer_setup=='Sandnesfj_Al':
            # Use values from Simonsen et al (2019b)
            self.num_lmmanion    = self.specie_name2num('LMManion')
            self.num_lmmcation   = self.specie_name2num('LMMcation')
            self.num_humcol      = self.specie_name2num('Humic colloid')
            self.num_polymer     = self.specie_name2num('Polymer')
            self.num_prev        = self.specie_name2num('Particle reversible')
            self.num_srev        = self.specie_name2num('Sediment reversible')

            Dc         = self.get_config('chemical:transformations:Dc')

            self.salinity_intervals = [0,1,10,20]

            # Resize transfer rates array
            self.transfer_rates = np.zeros([len(self.salinity_intervals),self.transfer_rates.shape[0],self.transfer_rates.shape[1]])

            # Salinity interval 0-1 psu
            self.transfer_rates[0,self.num_lmmcation, self.num_humcol]    = 1.2e-5
            self.transfer_rates[0,self.num_lmmcation, self.num_prev]      = 4.e-6
            self.transfer_rates[0,self.num_humcol,    self.num_lmmcation] = .3*Dc
            self.transfer_rates[0,self.num_humcol,    self.num_prev]      = 2.e-6
            self.transfer_rates[0,self.num_prev,      self.num_lmmcation] = .3*Dc
            self.transfer_rates[0,self.num_srev,      self.num_lmmcation] = .03*Dc

            # Salinity interval 1-10 psu
            self.transfer_rates[1,self.num_lmmcation, self.num_humcol]    = 1.e-5
            self.transfer_rates[1,self.num_lmmcation, self.num_prev]      = 3.e-6
            self.transfer_rates[1,self.num_lmmcation, self.num_polymer]   = 1.2e-4
            self.transfer_rates[1,self.num_humcol,    self.num_lmmcation] = 7.*Dc
            self.transfer_rates[1,self.num_humcol,    self.num_prev]      = 4.e-6
            self.transfer_rates[1,self.num_prev,      self.num_lmmcation] = .5*Dc
            self.transfer_rates[1,self.num_srev,      self.num_lmmcation] = .05*Dc
            self.transfer_rates[1,self.num_lmmanion,  self.num_polymer]   = 5.e-6
            self.transfer_rates[1,self.num_polymer,   self.num_lmmanion]  = 12.*Dc
            self.transfer_rates[1,self.num_polymer,   self.num_prev]      = 2.4e-5

            # Salinity interval 10-20 psu
            self.transfer_rates[2,self.num_lmmcation, self.num_humcol]    = 8.e-6
            self.transfer_rates[2,self.num_lmmcation, self.num_prev]      = 2.e-6
            self.transfer_rates[2,self.num_lmmcation, self.num_polymer]   = 1.4e-4
            self.transfer_rates[2,self.num_humcol,    self.num_lmmcation] = 7.*Dc
            self.transfer_rates[2,self.num_humcol,    self.num_prev]      = 6.e-6
            self.transfer_rates[2,self.num_prev,      self.num_lmmcation] = .6*Dc
            self.transfer_rates[2,self.num_srev,      self.num_lmmcation] = .06*Dc
            self.transfer_rates[2,self.num_lmmanion,  self.num_polymer]   = 5.e-6
            self.transfer_rates[2,self.num_polymer,   self.num_lmmanion]  = 12.*Dc
            self.transfer_rates[2,self.num_polymer,   self.num_prev]      = 6.e-5

            # Salinity interval >20 psu
            self.transfer_rates[3,self.num_lmmcation, self.num_humcol]    = 6.e-6
            self.transfer_rates[3,self.num_lmmcation, self.num_prev]      = 1.8e-6
            self.transfer_rates[3,self.num_lmmcation, self.num_polymer]   = 1.5e-4
            self.transfer_rates[3,self.num_humcol,    self.num_lmmcation] = 7.*Dc
            self.transfer_rates[3,self.num_humcol,    self.num_prev]      = 1.e-5
            self.transfer_rates[3,self.num_prev,      self.num_lmmcation] = .8*Dc
            self.transfer_rates[3,self.num_srev,      self.num_lmmcation] = .08*Dc
            self.transfer_rates[3,self.num_lmmanion,  self.num_polymer]   = 5.e-6
            self.transfer_rates[3,self.num_polymer,   self.num_lmmanion]  = 12.*Dc
            self.transfer_rates[3,self.num_polymer,   self.num_prev]      = 8.e-5




        else:
            logger.ERROR('No transfer setup available')


        # Set diagonal to 0. (not possible to transform to present specie)
        if len(self.transfer_rates.shape) == 3:
            for ii in range(self.transfer_rates.shape[0]):
                np.fill_diagonal(self.transfer_rates[ii,:,:],0.)
        else:
            np.fill_diagonal(self.transfer_rates,0.)

#         # HACK :
#         self.transfer_rates[:] = 0.
#         print ('\n ###### \n IMPORTANT:: \n transfer rates have been hacked! \n#### \n ')

        logger.debug('nspecies: %s' % self.nspecies)
        logger.debug('Transfer rates:\n %s' % self.transfer_rates)








    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for Pelagic Egg

        according to
        S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81  # ms-2

        # Particle properties that determine settling velocity
        partsize = self.elements.diameter
        # prepare interpolation of temp, salt
        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.uint8), 0)
            lower = np.minimum(upper+1, Tprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                weight_upper + \
                Tprofiles[lower, range(Tprofiles.shape[1])] * \
                (1-weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                weight_upper + \
                Sprofiles[lower, range(Sprofiles.shape[1])] * \
                (1-weight_upper)

        DENSw = self.sea_water_density(T=T0, S=S0)
        DENSpart = self.elements.density
        dr = DENSw-DENSpart  # density difference

        # water viscosity
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0)
        # ~0.0014 kg m-1 s-1

        # terminal velocity for low Reynolds numbers
        W = (1.0/my_w)*(1.0/18.0)*g*partsize**2 * dr

        #W=np.zeros_like(W) #Setting to zero for debugging

        self.elements.terminal_velocity = W

        self.elements.terminal_velocity = W * self.elements.moving


    def update_transfer_rates(self):
        '''Pick out the correct row from transfer_rates for each element. Modify the
        transfer rates according to local environmental conditions '''

        transfer_setup=self.get_config('chemical:transfer_setup')
        if transfer_setup == 'metals' or \
         transfer_setup=='custom' or \
         transfer_setup=='137Cs_rev'or \
         transfer_setup=='organics':
            self.elements.transfer_rates1D = self.transfer_rates[self.elements.specie,:]
            diss       = self.get_config('chemical:transformations:dissociation')

            # Updating desorption rates according to local temperature, salinity, pH

            if transfer_setup=='organics' and diss=='nondiss':
                # filtering out zero values from temperature and salinity
                # TODO: Find out if problem is in the reader or in the data
                temperature=self.environment.sea_water_temperature
                temperature[temperature==0]=np.median(temperature)

                salinity=self.environment.sea_water_salinity
                salinity[salinity==0]=np.median(salinity)

                KOWTref    = self.get_config('chemical:transformations:TrefKOW')
                DH_KOC_Sed = self.get_config('chemical:transformations:DeltaH_KOC_Sed')
                DH_KOC_DOM = self.get_config('chemical:transformations:DeltaH_KOC_DOM')
                Setchenow  = self.get_config('chemical:transformations:Setchenow')

                tempcorrSed = self.tempcorr("Arrhenius",DH_KOC_Sed,temperature,KOWTref)
                tempcorrDOM = self.tempcorr("Arrhenius",DH_KOC_DOM,temperature,KOWTref)
                salinitycorr = self.salinitycorr(Setchenow,temperature,salinity)

                # Temperature and salinity correction for desorption rates (inversely proportional to Kd)

                self.elements.transfer_rates1D[self.elements.specie==self.num_humcol,self.num_lmm] = \
                    self.k21_0 / tempcorrDOM[self.elements.specie==self.num_humcol] / salinitycorr[self.elements.specie==self.num_humcol]

                self.elements.transfer_rates1D[self.elements.specie==self.num_prev,self.num_lmm] = \
                    self.k31_0 / tempcorrSed[self.elements.specie==self.num_prev] / salinitycorr[self.elements.specie==self.num_prev]

                self.elements.transfer_rates1D[self.elements.specie==self.num_srev,self.num_lmm] = \
                    self.k41_0 / tempcorrSed[self.elements.specie==self.num_srev] / salinitycorr[self.elements.specie==self.num_srev]

            elif transfer_setup=='organics' and diss!='nondiss':
                # Select elements for updating trasfer rates in sediments, SPM, and DOM

                #Sediments
                S =   (self.elements.specie == self.num_srev) \
                    + (self.elements.specie == self.num_ssrev)

                SPM = (self.elements.specie == self.num_prev)

                DOM = (self.elements.specie == self.num_humcol)

                pH_sed = self.environment.pH_sediment[S]
                # pH_sed[pH_sed==0]=np.median(pH_sed)

                pH_water_SPM=self.environment.sea_water_ph_reported_on_total_scale[SPM]
                # pH_water_SPM[pH_water_SPM==0]=np.median(TW)

                pH_water_DOM=self.environment.sea_water_ph_reported_on_total_scale[DOM]
                # pH_water_DOM[pH_water_DOM==0]=np.median(pH_water_DOM)

                pKa_acid   = self.get_config('chemical:transformations:pKa_acid')
                if pKa_acid < 0:
                    raise ValueError("pKa_acid must be positive")
                    # print("pKa_acid must be positive")
                    # UserWarning(("pKa_acid must be positive"))
                else:
                    pass

                pKa_base   = self.get_config('chemical:transformations:pKa_base')
                if pKa_base < 0:
                    raise ValueError("pKa_base must be positive")
                    # print("pKa_base must be positive")
                    # UserWarning(("pKa_base must be positive"))
                else:
                    pass

                KOW = 10**self.get_config('chemical:transformations:LogKOW')

                KOC_sed_n = self.get_config('chemical:transformations:KOC_sed')
                if KOC_sed_n < 0:
                    KOC_sed_n = 2.62 * KOW**0.82   # (L/KgOC), Park and Clough, 2014 (334)/Org2C TO DO Add if choice between input and estimation
                else:
                    pass

                KOC_DOM_n = self.get_config('chemical:transformations:KOC_DOM')
                if KOC_DOM_n < 0:
                    KOC_DOM_n = 2.88 * KOW**0.67   # (L/KgOC), Park and Clough, 2014 TO DO Add if choice between input and estimation
                else:
                    pass

                fOC_SPM    = self.get_config('chemical:transformations:fOC_SPM')       # typical values from 0.01 to 0.1 gOC/g
                fOC_sed    = self.get_config('chemical:transformations:fOC_sed')
                Org2C      = 0.526  # kgOC/KgOM

                # Calculate original KOC_Values
                # TO DO: Store directly KOC values

                KOC_sed_initial = (self.Kd_sed)/fOC_sed # L/Kg / KgOC/Kg = L/KgOC
                KOC_SPM_initial = (self.Kd_SPM)/fOC_SPM # L/Kg / KgOC/Kg = L/KgOC
                KOC_DOM_initial = (self.Kd_DOM)/Org2C

                # filtering out zero values from temperature and salinity
                # TODO: Find out if problem is in the reader or in the data
                temperature=self.environment.sea_water_temperature
                temperature[temperature==0]=np.median(temperature)

                salinity=self.environment.sea_water_salinity
                salinity[salinity==0]=np.median(salinity)

                KOWTref    = self.get_config('chemical:transformations:TrefKOW')
                DH_KOC_Sed = self.get_config('chemical:transformations:DeltaH_KOC_Sed')
                DH_KOC_DOM = self.get_config('chemical:transformations:DeltaH_KOC_DOM')
                Setchenow  = self.get_config('chemical:transformations:Setchenow')

                tempcorrSed = self.tempcorr("Arrhenius",DH_KOC_Sed,temperature,KOWTref)
                tempcorrDOM = self.tempcorr("Arrhenius",DH_KOC_DOM,temperature,KOWTref)
                salinitycorr = self.salinitycorr(Setchenow,temperature,salinity)

                KOC_sedcorr = self.calc_KOC_sedcorr(KOC_sed_initial, KOC_sed_n, pKa_acid, pKa_base, KOW, pH_sed, diss)
                KOC_watcorrSPM = self.calc_KOC_watcorrSPM(KOC_SPM_initial, KOC_sed_n, pKa_acid, pKa_base, KOW, pH_water_SPM, diss)
                KOC_watcorrDOM = self.calc_KOC_watcorrDOM(KOC_DOM_initial, KOC_DOM_n, pKa_acid, pKa_base, KOW, pH_water_DOM, diss)

                # Temperature and salinity correction for desorption rates (inversely proportional to Kd)

                ####

                self.elements.transfer_rates1D[self.elements.specie==self.num_humcol,self.num_lmm] = \
                    self.k21_0 * KOC_watcorrDOM / tempcorrDOM[self.elements.specie==self.num_humcol] / salinitycorr[self.elements.specie==self.num_humcol]

                self.elements.transfer_rates1D[self.elements.specie==self.num_prev,self.num_lmm] = \
                    self.k31_0 * KOC_watcorrSPM / tempcorrSed[self.elements.specie==self.num_prev] / salinitycorr[self.elements.specie==self.num_prev]

                self.elements.transfer_rates1D[self.elements.specie==self.num_srev,self.num_lmm] = \
                    self.k41_0 * KOC_sedcorr / tempcorrSed[self.elements.specie==self.num_srev] / salinitycorr[self.elements.specie==self.num_srev]

            # Updating sorption rates

            if transfer_setup=='organics':

                # Updating sorption rates according to local SPM concentration

                concSPM=self.environment.spm * 1e-6 # (Kg/L) from (g/m3)

                # Apply SPM concentration profile if SPM reader has not depth coordinate
                # SPM concentration is kept constant to surface value in the mixed layer
                # Exponentially decreasing with depth below the mixed layers
                if not self.SPM_vertical_levels_given:
                    lowerMLD = self.elements.z < -self.environment.ocean_mixed_layer_thickness
                    #concSPM[lowerMLD] = concSPM[lowerMLD]/2
                    concSPM[lowerMLD] = concSPM[lowerMLD] * np.exp(
                        -(self.elements.z[lowerMLD]+self.environment.ocean_mixed_layer_thickness[lowerMLD])
                        *np.log(0.5)/self.get_config('chemical:particle_concentration_half_depth')
                        )

                self.elements.transfer_rates1D[self.elements.specie==self.num_lmm,self.num_prev] = \
                    self.k_ads * concSPM[self.elements.specie==self.num_lmm]      # k13

            if transfer_setup == 'metals':

                # Updating sorption rates according to local SPM concentration and salinity

                concSPM=self.environment.spm * 1e-3 # (Kg/m3) from (g/m3)

                salinity=self.environment.sea_water_salinity

                # Apply SPM concentration profile if SPM reader has not depth coordinate
                # SPM concentration is kept constant to surface value in the mixed layer
                # Exponentially decreasing with depth below the mixed layers
                if not self.SPM_vertical_levels_given:
                    lowerMLD = self.elements.z < -self.environment.ocean_mixed_layer_thickness
                    #concSPM[lowerMLD] = concSPM[lowerMLD]/2
                    concSPM[lowerMLD] = concSPM[lowerMLD] * np.exp(
                        -(self.elements.z[lowerMLD]+self.environment.ocean_mixed_layer_thickness[lowerMLD])
                        *np.log(0.5)/self.get_config('chemical:particle_concentration_half_depth')
                        )

                Kd0         = self.get_config('chemical:transformations:Kd') # (m3/Kg)
                S0          = self.get_config('chemical:transformations:S0') # (PSU)
                Dc          = self.get_config('chemical:transformations:Dc') # (1/s)
                sed_L       = self.get_config('chemical:sediment:mixing_depth')     # sediment mixing depth (m)
                sed_dens    = self.get_config('chemical:sediment:density') # default particle density (kg/m3)
                sed_f       = self.get_config('chemical:sediment:effective_fraction')      # fraction of effective sorbents
                sed_phi     = self.get_config('chemical:sediment:corr_factor')      # sediment correction factor
                sed_poro    = self.get_config('chemical:sediment:porosity')      # sediment porosity
                sed_H       = self.get_config('chemical:sediment:layer_thickness')      # thickness of seabed interaction layer (m)

                # Adjust Kd for salinity according to Perianez 2018 https://doi.org/10.1016/j.jenvrad.2018.02.014
                if S0>0:
                    Kd=Kd0*(S0+salinity[self.elements.specie==self.num_lmm])/S0

                self.elements.transfer_rates1D[self.elements.specie==self.num_lmm,self.num_prev] = \
                    Dc * Kd * concSPM[self.elements.specie==self.num_lmm]      # k13

                self.elements.transfer_rates1D[self.elements.specie==self.num_lmm,self.num_srev] = \
                    Dc * Kd * sed_L * sed_dens * (1.-sed_poro) * sed_f * sed_phi / sed_H


            if transfer_setup=='organics':

                # Updating sorption rates according to local DOC concentration

                concDOM = self.environment.doc * 12e-6 / 1.025 / 0.526 * 1e-3 # (Kg[OM]/L) from (umol[C]/Kg)

                # Apply DOC concentration profile if DOC reader has not depth coordinate
                # DOC concentration is kept constant to surface value in the mixed layer
                # Exponentially decreasing with depth below the mixed layers

                if not self.DOC_vertical_levels_given:
                    lowerMLD = self.elements.z < -self.environment.ocean_mixed_layer_thickness
                    #concDOM[lowerMLD] = concDOM[lowerMLD]/2
                    concDOM[lowerMLD] = concDOM[lowerMLD] * np.exp(
                        -(self.elements.z[lowerMLD]+self.environment.ocean_mixed_layer_thickness[lowerMLD])
                        *np.log(0.5)/self.get_config('chemical:doc_concentration_half_depth')
                        )

                self.elements.transfer_rates1D[self.elements.specie==self.num_lmm,self.num_humcol] = \
                    self.k_ads * concDOM[self.elements.specie==self.num_lmm]      # k12

            if self.get_config('chemical:species:Sediment_reversible'):
                # Only LMM chemicals close to seabed are allowed to interact with sediments
                # minimum height/maximum depth for each particle
                Zmin = -1.*self.environment.sea_floor_depth_below_sea_level
                interaction_thick = self.get_config('chemical:sediment:layer_thickness')      # thickness of seabed interaction layer (m)
                dist_to_seabed = self.elements.z - Zmin
                self.elements.transfer_rates1D[(self.elements.specie == self.num_lmm) &
                                 (dist_to_seabed > interaction_thick), self.num_srev] = 0.

        elif transfer_setup=='Sandnesfj_Al':
            sal = self.environment.sea_water_salinity
            sali = np.searchsorted(self.salinity_intervals, sal) - 1
            self.elements.transfer_rates1D = self.transfer_rates[sali,self.elements.specie,:]



    def update_partitioning(self):
        '''Check if transformation processes shall occur
        Do transformation (change value of self.elements.specie)
        Update element properties for the transformed elements
        '''

        specie_in  = self.elements.specie.copy()    # for storage of the initial partitioning
        specie_out = self.elements.specie.copy()    # for storage of the final partitioning
        deltat = self.time_step.total_seconds()     # length of a time step
        phaseshift = np.array(self.num_elements_active()*[False])  # Denotes which trajectory that shall be transformed

        p = 1. - np.exp(-self.elements.transfer_rates1D*deltat)  # Probability for transformation

        psum = np.sum(p,axis=1)

        ran1=np.random.random(self.num_elements_active())

        # Transformation where ran1 < total probability for transformation
        phaseshift[ ran1 < psum ] = True

        logger.info('Number of transformations: %s' % sum(phaseshift))
        if sum(phaseshift) == 0:
            return

        ran4 = np.random.random(sum(phaseshift)) # New random number to decide which specie to end up in

        ttmp=[]  # list for storing the out specie
        # Loop through each trajectory
        for ii in range(sum(phaseshift)):
            # Compare random number to the relative probability for each transfer process
            ttmp.append(np.searchsorted(np.cumsum(p[phaseshift][ii]/psum[phaseshift][ii]),ran4[ii]))
        specie_out[phaseshift] = np.array(ttmp)


        # Set the new partitioning
        self.elements.specie = specie_out

        logger.debug('old species: %s' % specie_in[phaseshift])
        logger.debug('new species: %s' % specie_out[phaseshift])


        for iin in range(self.nspecies):
            for iout in range(self.nspecies):
                self.ntransformations[iin,iout]+=sum((specie_in[phaseshift]==iin) & (specie_out[phaseshift]==iout))

        logger.debug('Number of transformations total:\n %s' % self.ntransformations )


        # Update Chemical properties after transformations
        self.update_chemical_diameter(specie_in, specie_out)
        self.sorption_to_sediments(specie_in, specie_out)
        self.desorption_from_sediments(specie_in, specie_out)





    def sorption_to_sediments(self,sp_in=None,sp_out=None):
        '''Update Chemical properties  when sorption to sediments occurs'''


        # Set z to local sea depth
        if self.get_config('chemical:species:LMM'):
            self.elements.z[(sp_out==self.num_srev) & (sp_in==self.num_lmm)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_srev) & (sp_in==self.num_lmm)]
            self.elements.moving[(sp_out==self.num_srev) & (sp_in==self.num_lmm)] = 0
        if self.get_config('chemical:species:LMMcation'):
            self.elements.z[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)]
            self.elements.moving[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)] = 0
        # avoid setting positive z values
        if np.nansum(self.elements.z>0):
            logger.debug('Number of elements lowered down to sea surface: %s' % np.nansum(self.elements.z>0))
        self.elements.z[self.elements.z > 0] = 0



    def desorption_from_sediments(self,sp_in=None,sp_out=None):
        '''Update Chemical properties when desorption from sediments occurs'''

        desorption_depth = self.get_config('chemical:sediment:desorption_depth')
        std = self.get_config('chemical:sediment:desorption_depth_uncert')


        if self.get_config('chemical:species:LMM'):
            self.elements.z[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] + desorption_depth
            self.elements.moving[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] = 1
            if std > 0:
                logger.debug('Adding uncertainty for desorption from sediments: %s m' % std)
                self.elements.z[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] += np.random.normal(
                        0, std, sum((sp_out==self.num_lmm) & (sp_in==self.num_srev)))
        if self.get_config('chemical:species:LMMcation'):
            self.elements.z[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] + desorption_depth
            self.elements.moving[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] = 1
            if std > 0:
                logger.debug('Adding uncertainty for desorption from sediments: %s m' % std)
                self.elements.z[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] += np.random.normal(
                        0, std, sum((sp_out==self.num_lmmcation) & (sp_in==self.num_srev)))
        # avoid setting positive z values
        if np.nansum(self.elements.z>0):
            logger.debug('Number of elements lowered down to sea surface: %s' % np.nansum(self.elements.z>0))
        self.elements.z[self.elements.z > 0] = 0





    def update_chemical_diameter(self,sp_in=None,sp_out=None):
        '''Update the diameter of the chemicals when specie is changed'''


        dia_part=self.get_config('chemical:particle_diameter')
        dia_diss=self.get_config('chemical:dissolved_diameter')


        # Transfer to reversible particles
        self.elements.diameter[(sp_out==self.num_prev) & (sp_in!=self.num_prev)] = dia_part

        # TODO Choose a proper diameter for aggregated particles
        if self.get_config('chemical:species:Humic_colloid'):
            self.elements.diameter[(sp_out==self.num_prev) & (sp_in==self.num_humcol)] = dia_part/2

        logger.debug('Updated particle diameter for %s elements' % len(self.elements.diameter[(sp_out==self.num_prev) & (sp_in!=self.num_prev)]))

        std = self.get_config('chemical:particle_diameter_uncertainty')
        if std > 0:
            logger.debug('Adding uncertainty for particle diameter: %s m' % std)
            self.elements.diameter[(sp_out==self.num_prev) & (sp_in!=self.num_prev)] += np.random.normal(
                    0, std, sum((sp_out==self.num_prev) & (sp_in!=self.num_prev)))
        # Transfer to slowly reversible particles
        if self.get_config('chemical:slowly_fraction'):
            self.elements.diameter[(sp_out==self.num_psrev) & (sp_in!=self.num_psrev)] = dia_part
            if std > 0:
                logger.debug('Adding uncertainty for slowly rev particle diameter: %s m' % std)
                self.elements.diameter[(sp_out==self.num_psrev) & (sp_in!=self.num_psrev)] += np.random.normal(
                    0, std, sum((sp_out==self.num_psrev) & (sp_in!=self.num_psrev)))

        # Transfer to irreversible particles
        if self.get_config('chemical:irreversible_fraction'):
            self.elements.diameter[(sp_out==self.num_pirrev) & (sp_in!=self.num_pirrev)] = dia_part
            if std > 0:
                logger.debug('Adding uncertainty for irrev particle diameter: %s m' % std)
                self.elements.diameter[(sp_out==self.num_pirrev) & (sp_in!=self.num_pirrev)] += np.random.normal(
                    0, std, sum((sp_out==self.num_pirrev) & (sp_in!=self.num_pirrev)))

        # Transfer to LMM
        if self.get_config('chemical:species:LMM'):
            self.elements.diameter[(sp_out==self.num_lmm) & (sp_in!=self.num_lmm)] = dia_diss
        if self.get_config('chemical:species:LMManion'):
            self.elements.diameter[(sp_out==self.num_lmmanion) & (sp_in!=self.num_lmmanion)] = dia_diss
        if self.get_config('chemical:species:LMMcation'):
            self.elements.diameter[(sp_out==self.num_lmmcation) & (sp_in!=self.num_lmmcation)] = dia_diss

        # Transfer to colloids
        if self.get_config('chemical:species:Colloid'):
            self.elements.diameter[(sp_out==self.num_col) & (sp_in!=self.num_col)] = dia_diss
        if self.get_config('chemical:species:Humic_colloid'):
            self.elements.diameter[(sp_out==self.num_humcol) & (sp_in!=self.num_humcol)] = dia_diss
        if self.get_config('chemical:species:Polymer'):
            self.elements.diameter[(sp_out==self.num_polymer) & (sp_in!=self.num_polymer)] = dia_diss




    def bottom_interaction(self,Zmin=None):
        ''' Change partitioning of chemicals that reach bottom due to settling.
        particle specie -> sediment specie '''
        if not  ((self.get_config('chemical:species:Particle_reversible')) &
                  (self.get_config('chemical:species:Sediment_reversible')) or
                  (self.get_config('chemical:slowly_fraction')) or
                  (self.get_config('chemical:irreversible_fraction'))):
            return

        bottom = np.array(np.where(self.elements.z <= Zmin)[0])
        kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_prev)[0])
        self.elements.specie[bottom[kktmp]] = self.num_srev
        self.ntransformations[self.num_prev,self.num_srev]+=len(kktmp)
        self.elements.moving[bottom[kktmp]] = 0
        if self.get_config('chemical:slowly_fraction'):
            kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_psrev)[0])
            self.elements.specie[bottom[kktmp]] = self.num_ssrev
            self.ntransformations[self.num_psrev,self.num_ssrev]+=len(kktmp)
            self.elements.moving[bottom[kktmp]] = 0
        if self.get_config('chemical:irreversible_fraction'):
            kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_pirrev)[0])
            self.elements.specie[bottom[kktmp]] = self.num_sirrev
            self.ntransformations[self.num_pirrev,self.num_sirrev]+=len(kktmp)
            self.elements.moving[bottom[kktmp]] = 0


    def resuspension(self):
        """ Simple method to estimate the resuspension of sedimented particles,
        checking whether the current speed near the bottom is above a critical velocity
        Sediment species -> Particle specie
        """
        # Exit function if particles and sediments not are present
        if not  ((self.get_config('chemical:species:Particle_reversible')) &
                  (self.get_config('chemical:species:Sediment_reversible'))):
            return

        specie_in = self.elements.specie.copy()

        critvel = self.get_config('chemical:sediment:resuspension_critvel')
        resusp_depth = self.get_config('chemical:sediment:resuspension_depth')
        std = self.get_config('chemical:sediment:resuspension_depth_uncert')

        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level
        x_vel = self.environment.x_sea_water_velocity
        y_vel = self.environment.y_sea_water_velocity
        speed = np.sqrt(x_vel*x_vel + y_vel*y_vel)
        bottom = (self.elements.z <= Zmin)

        resusp = ( (bottom) & (speed >= critvel) )
        if self.get_config('chemical:slowly_fraction'):
            resusp = ( resusp & (self.elements.specie!=self.num_ssrev) )    # Prevent ssrev (buried) to be resuspended
                                                                        # TODO buried sediment should be a new specie
        logger.info('Number of resuspended particles: {}'.format(np.sum(resusp)))
        self.elements.moving[resusp] = 1

        self.elements.z[resusp] = Zmin[resusp] + resusp_depth
        if std > 0:
            logger.debug('Adding uncertainty for resuspension from sediments: %s m' % std)
            self.elements.z[resusp] += np.random.normal(
                        0, std, sum(resusp))
        # avoid setting positive z values
        if np.nansum(self.elements.z>0):
            logger.debug('Number of elements lowered down to sea surface: %s' % np.nansum(self.elements.z>0))
        self.elements.z[self.elements.z > 0] = 0

        self.ntransformations[self.num_srev,self.num_prev]+=sum((resusp) & (self.elements.specie==self.num_srev))
        self.elements.specie[(resusp) & (self.elements.specie==self.num_srev)] = self.num_prev
        if self.get_config('chemical:slowly_fraction'):
            self.ntransformations[self.num_ssrev,self.num_psrev]+=sum((resusp) & (self.elements.specie==self.num_ssrev))
            self.elements.specie[(resusp) & (self.elements.specie==self.num_ssrev)] = self.num_psrev

        if self.get_config('chemical:irreversible_fraction'):
            self.ntransformations[self.num_sirrev,self.num_pirrev]+=sum((resusp) & (self.elements.specie==self.num_sirrev))
            self.elements.specie[(resusp) & (self.elements.specie==self.num_sirrev)] = self.num_pirrev

        specie_out = self.elements.specie.copy()
        self.update_chemical_diameter(specie_in, specie_out)

    def degradation(self):
        '''degradation.'''

        if self.get_config('chemical:transformations:degradation') is True:
            if self.get_config('chemical:transformations:degradation_mode')=='OverallRateConstants':
                # TODO: Rearrange code. Calculations here are for overall degradation including
                # degradation, photodegradation, and hydrolysys

                logger.debug('Calculating overall degradation using overall rate constants')

                degraded_now = np.zeros(self.num_elements_active())

                # Degradation in the water
                k_W_tot = -np.log(0.5)/(self.get_config('chemical:transformations:t12_W_tot')*(60*60)) # (1/s)
                Tref_kWt = self.get_config('chemical:transformations:Tref_kWt')
                DH_kWt = self.get_config('chemical:transformations:DeltaH_kWt')

                W =   (self.elements.specie == self.num_lmm) \
                    + (self.elements.specie == self.num_humcol)

                TW=self.environment.sea_water_temperature[W]
                TW[TW==0]=np.median(TW)

                k_W_fin = k_W_tot * self.tempcorr("Arrhenius",DH_kWt,TW,Tref_kWt)

                degraded_now[W] = self.elements.mass[W] * (1-np.exp(-k_W_fin * self.time_step.total_seconds()))

                # Degradation in the sediments

                k_S_tot = -np.log(0.5)/(self.get_config('chemical:transformations:t12_S_tot')*(60*60)) # (1/s)
                Tref_kSt = self.get_config('chemical:transformations:Tref_kSt')
                DH_kSt = self.get_config('chemical:transformations:DeltaH_kSt')

                S =   (self.elements.specie == self.num_srev) \
                    + (self.elements.specie == self.num_ssrev)

                TS=self.environment.sea_water_temperature[S]
                TS[TS==0]=np.median(TS)

                k_S_fin = k_S_tot * self.tempcorr("Arrhenius",DH_kSt,TS,Tref_kSt)

                degraded_now[S] = self.elements.mass[S] * (1-np.exp(-k_S_fin * self.time_step.total_seconds()))

            self.elements.mass_degraded_water[W] = self.elements.mass_degraded_water[W] + degraded_now[W]
            self.elements.mass_degraded_sediment[S] = self.elements.mass_degraded_sediment[S] + degraded_now[S]

            self.elements.mass_degraded = self.elements.mass_degraded + degraded_now
            self.elements.mass = self.elements.mass - degraded_now
            self.deactivate_elements(self.elements.mass < (self.elements.mass + self.elements.mass_degraded + self.elements.mass_volatilized)/100,
                                     reason='removed')

            #to_deactivate = self.elements.mass < (self.elements.mass + self.elements.mass_degraded + self.elements.mass_volatilized)/100
            #vol_morethan_degr = self.elements.mass_degraded >= self.elements.mass_volatilized
            #
            #self.deactivate_elements(to_deactivate +  vol_morethan_degr, reason='volatilized')
            #self.deactivate_elements(to_deactivate + ~vol_morethan_degr, reason='degraded')

        else:
            pass


    def volatilization(self):
        if self.get_config('chemical:transformations:volatilization') is True:
            logger.debug('Calculating: volatilization')
            volatilized_now = np.zeros(self.num_elements_active())

            MolWtCO2=44
            MolWtH2O=18
            MolWt=self.get_config('chemical:transformations:MolWt')
            wind=5                  # (m/s) (to read from atmosferic forcing)
            mixedlayerdepth=50      # m     (to read from ocean forcing)
            Undiss_n=1              # 1 for PAHs

            Henry=self.get_config('chemical:transformations:Henry') # (atm m3/mol)

            Vp=self.get_config('chemical:transformations:Vpress')
            Tref_Vp=self.get_config('chemical:transformations:Tref_Vpress')
            DH_Vp=self.get_config('chemical:transformations:DeltaH_Vpress')

            Slb=self.get_config('chemical:transformations:Solub')
            Tref_Slb=self.get_config('chemical:transformations:Tref_Solub')
            DH_Slb=self.get_config('chemical:transformations:DeltaH_Solub')

            R=8.206e-05 #(atm m3)/(mol K)

            mixedlayerdepth = self.environment.ocean_mixed_layer_thickness

            # mask of dissolved elements within mixed layer
            W =     (self.elements.specie == self.num_lmm) \
                  * (-self.elements.z <= mixedlayerdepth)
                    # does volatilization apply only to num_lmm?
                    # check

            mixedlayerdepth = mixedlayerdepth[W]

            T=self.environment.sea_water_temperature[W]
            T[T==0]=np.median(T)                            # temporary fix for missing values

            S=self.environment.sea_water_salinity[W]

            wind=(self.environment.x_wind[W]**2 + self.environment.y_wind[W]**2)**.5

            Henry=(      (Vp * self.tempcorr("Arrhenius",DH_Vp,T,Tref_Vp)))   \
                       / (Slb *  self.tempcorr("Arrhenius",DH_Slb,T,Tref_Slb))  \
                       * MolWt / 101325.    # atm m3 mol-1

            #k_S_fin = k_S_tot * self.tempcorr("Arrhenius",DH_kSt,TS,Tref_kSt)

            # Calculate mass transfer coefficient water side
            # Schwarzenbach et al., 2016 Eq.(19-20)

            MTCw = ((9e-4)+(7.2e-6*wind**3)) * (MolWtCO2/MolWt)**0.25 / Undiss_n

            # Calculate mass transfer coefficient air side
            # Schwarzenbach et al., 2016 Eq.(19-17)(19-18)(19-19)

            # Simple
            #MTCaH2O = 0.1 + 0.11 * wind

            # More complex
            Sca_H2O = 0.62                                  # 0.6 in the book. check
            MTCaH2O = 0.1 + wind*(6.1+0.63*wind)**0.5 \
                /(13.3*(Sca_H2O)**0.5 + (6.1e-4+(6.3e-5)*wind)**-0.5 -5 + 1.25*np.log(Sca_H2O) )

            MTCa = MTCaH2O * (MolWtH2O/MolWt)**(1/3)

            # Calculate overall volatilization mass tansfer coefficient

            HenryLaw = Henry * (1 + 0.01143 * S) / ( R * (T+273.15) )

            MTCvol = 1 / ( 1/MTCw + 1/(MTCa * HenryLaw))     # (cm/s)
            #mixedlayerdepth = self.environment.ocean_mixed_layer_thickness[W]
            #Thick = np.clip(self.environment.sea_floor_depth_below_sea_level[W],0,mixedlayerdepth) # (m)
            Thick = mixedlayerdepth

            # Degubbing information to screen
            #print('################### Volatilization-info ##################')
            #print('Mixed Layer   ',len(mixedlayerdepth),min(mixedlayerdepth),max(mixedlayerdepth),'m')
            #print('Temperature   ',len(T),min(T),max(T),'C')
            #print('Salinity      ',len(S),min(S),max(S))
            #print('Henry         ',len(Henry),min(Henry),max(Henry),'atm m3 / mol')
            #print('HenryLaw      ',len(HenryLaw),min(HenryLaw),max(HenryLaw))
            #print('wind          ',len(wind),min(wind),max(wind), 'm/s')
            #print('MTCa          ',len(MTCa),min(MTCa),max(MTCa),'cm/s')
            #print('MTCw          ',len(MTCw),min(MTCw),max(MTCw),'cm/s')
            #print('MTCa*HenryLaw ',len(MTCa*HenryLaw),min(MTCa*HenryLaw),max(MTCa*HenryLaw),'cm/s')
            #print('MTCvol        ',len(MTCvol),min(MTCvol),max(MTCvol),'cm/s')

            K_volatilization = 0.01 * MTCvol / Thick # (1/s)

            #logger.debug('MTCa: %s cm/s' % MTCa)
            #logger.debug('MTCw: %s cm/s' % MTCw)
            #logger.debug('Henry: %s ' % HenryLaw)
            #logger.debug('MTCvol: %s cm/s' % MTCvol)
            #logger.debug('T: %s C' % T)
            #logger.debug('S: %s ' % S)
            #logger.debug('Thick: %s ' % Thick)

            volatilized_now[W] = self.elements.mass[W] * (1-np.exp(-K_volatilization * self.time_step.total_seconds()))

            self.elements.mass_volatilized = self.elements.mass_volatilized + volatilized_now
            self.elements.mass = self.elements.mass - volatilized_now
            self.deactivate_elements(self.elements.mass < (self.elements.mass + self.elements.mass_degraded + self.elements.mass_volatilized)/100,
                                     reason='removed')

            #to_deactivate = self.elements.mass < (self.elements.mass + self.elements.mass_degraded + self.elements.mass_volatilized)/100
            #vol_morethan_degr = self.elements.mass_degraded >= self.elements.mass_volatilized
            #
            #self.deactivate_elements(to_deactivate +  vol_morethan_degr, reason='volatilized')
            #self.deactivate_elements(to_deactivate + ~vol_morethan_degr, reason='degraded')


        else:
            pass

    def update(self):
        """Update positions and properties of Chemical particles."""

        # Workaround due to conversion of datatype
        self.elements.specie = self.elements.specie.astype(np.int32)

        # Degradation and Volatilization
        if self.get_config('chemical:transfer_setup')=='organics':
            self.degradation()
            self.volatilization()

        # Dynamic Partitioning
        if self.get_config('chemical:dynamic_partitioning') is True:
            self.update_transfer_rates()
            self.update_partitioning()

        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
        else:
            self.update_terminal_velocity()
            self.vertical_buoyancy()


        # Resuspension
        self.resuspension()
        logger.info('partitioning: {} {}'.format([sum(self.elements.specie==ii) for ii in range(self.nspecies)],self.name_species))



        # Horizontal advection
        self.advect_ocean_current()

        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()

        # Update transfer rates after last time step
        if      self.time == (self.expected_end_time - self.time_step) or \
                self.time == (self.expected_end_time) or \
                self.num_elements_active() == 0 :
            self.update_transfer_rates()






# ################
# POSTPROCESSING
    def simulation_summary(self, chemical_compound):
        '''Print a summary of the simulation: number of elements, number of transformations
        and final speciation
        '''

        print(chemical_compound)

        print('Final speciation:')
        for isp,sp in enumerate(self.name_species):
            print ('{:32}: {:>6}'.format(sp,sum(self.elements.specie==isp)))

        print('Number of transformations:')
        for isp in range(self.nspecies):
            print('{}'.format(['{:>9}'.format(np.int32(item)) for item in self.ntransformations[isp,:]]) )

        m_pre = sum(self.elements.mass)+sum(self.elements_deactivated.mass)
        m_deg = sum(self.elements.mass_degraded)+sum(self.elements_deactivated.mass_degraded)
        m_deg_w = sum(self.elements.mass_degraded_water)+sum(self.elements_deactivated.mass_degraded_water)
        m_deg_s = sum(self.elements.mass_degraded_sediment)+sum(self.elements_deactivated.mass_degraded_sediment)
        m_vol = sum(self.elements.mass_volatilized)+sum(self.elements_deactivated.mass_volatilized)
        m_tot = m_pre + m_deg + m_vol

        print('Mass balance:')
        print('mass preserved       :', m_pre * 1e-6,' g   ',m_pre/m_tot*100,'%')
        print('mass degraded        :', m_deg * 1e-6,' g   ',m_deg/m_tot*100,'%')
        print('     in water column :', m_deg_w * 1e-6,' g   ',m_deg_w/m_tot*100,'%')
        print('     in sediments    :', m_deg_s * 1e-6,' g   ',m_deg_s/m_tot*100,'%')
        print('mass volatilized     :', m_vol * 1e-6,' g   ',m_vol/m_tot*100,'%')


    def write_netcdf_chemical_density_map(self, filename, pixelsize_m='auto', zlevels=None,
                                              deltat=None,
                                              density_proj=None,
                                              llcrnrlon=None, llcrnrlat=None,
                                              urcrnrlon=None, urcrnrlat=None,
                                              mass_unit=None,
                                              time_avg_conc=False,
                                              horizontal_smoothing=False,
                                              smoothing_cells=0,
                                              reader_sea_depth=None,
                                              landmask_shapefile=None,
                                              origin_marker=None):
        '''Write netCDF file with map of Chemical species densities and concentrations'''

        from netCDF4 import Dataset, date2num #, stringtochar

        if landmask_shapefile is not None:
            if 'shape' in self.env.readers.keys():
                # removing previously stored landmask
                del self.env.readers['shape']
            # Adding new landmask
            from opendrift.readers import reader_shape
            custom_landmask = reader_shape.Reader.from_shpfiles(landmask_shapefile)
            self.add_reader(custom_landmask)
        elif 'global_landmask' not in self.env.readers.keys():
            from opendrift.readers import reader_global_landmask
            global_landmask = reader_global_landmask.Reader()
            self.add_reader(global_landmask)

        if reader_sea_depth is not None:
            from opendrift.readers import reader_netCDF_CF_generic
            reader_sea_depth = reader_netCDF_CF_generic.Reader(reader_sea_depth)
        else:
            print('A reader for ''sea_floor_depth_below_sea_level'' must be specified')
            import sys
            sys.exit()

        # Temporary workaround if self.nspecies and self.name_species are not defined
        # TODO Make sure that these are saved when the simulation data is saved to the ncdf file
        # Then this workaround can be removed
        if not hasattr(self,'nspecies'):
            self.nspecies=4
        if not hasattr(self,'name_species'):
            self.name_species = ['dissolved',
                                 'DOC',
                                 'SPM',
                                 'sediment']

        logger.info('Postprocessing: Write density and concentration to netcdf file')

        # Default bathymetry resolution 500x500. Can be increased (carefully) if high-res data is available and needed
        grid=np.meshgrid(np.linspace(llcrnrlon,urcrnrlon,500), np.linspace(llcrnrlat,urcrnrlat,500))
        self.conc_lon=grid[0]
        self.conc_lat=grid[1]

        self.conc_topo=reader_sea_depth.get_variables_interpolated_xy(['sea_floor_depth_below_sea_level'],
                x = self.conc_lon.flatten(),
                y = self.conc_lat.flatten(),
                time = reader_sea_depth.times[0] if reader_sea_depth.times is not None else None
                )[0]['sea_floor_depth_below_sea_level'].reshape(self.conc_lon.shape)

        if pixelsize_m == 'auto':
            lon, lat = self.get_lonlats()
            latspan = lat.max()-lat.min()
            pixelsize_m=30
            if latspan > .05:
                pixelsize_m = 50
            if latspan > .1:
                pixelsize_m = 300
            if latspan > .3:
                pixelsize_m = 500
            if latspan > .7:
                pixelsize_m = 1000
            if latspan > 2:
                pixelsize_m = 2000
            if latspan > 5:
                pixelsize_m = 4000


        if density_proj is None: # add default projection with equal-area property
            density_proj = pyproj.Proj('+proj=moll +ellps=WGS84 +lon_0=0.0')



        if mass_unit==None:
            mass_unit='microgram'  # default unit for chemicals

        z = self.get_property('z')[0]
        if not zlevels==None:
            zlevels = np.sort(zlevels)
            z_array = np.append(np.append(-10000, zlevels) , max(0,np.nanmax(z)))
        else:
            z_array = [min(-10000,np.nanmin(z)), max(0,np.nanmax(z))]
        logger.info('vertical grid boundaries: {}'.format(  [str(item) for item in z_array] ) )

        #
        # H is array containing number of elements within each box defined by lon_array, lat_array and z_array

        H, lon_array, lat_array = \
            self.get_chemical_density_array(pixelsize_m, z_array,
                                                density_proj=density_proj,
                                                llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                                                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                                                weight='mass',origin_marker=origin_marker)

        # calculating center point for eacxh pixel
        lon_array = (lon_array[:-1,:-1] + lon_array[1:,1:])/2
        lat_array = (lat_array[:-1,:-1] + lat_array[1:,1:])/2

        landmask = np.zeros_like(H[0,0,0,:,:])
        if landmask_shapefile is not None:
            landmask = self.env.readers['shape'].__on_land__(lon_array,lat_array)
        else:
            landmask = self.env.readers['global_landmask'].__on_land__(lon_array,lat_array)

        if horizontal_smoothing:
            # Compute horizontally smoother field
            logger.debug('H.shape: ' + str(H.shape))
            Hsm = np.zeros_like(H)
            for zi in range(len(z_array)-1):
                for sp in range(self.nspecies):
                    for ti in range(H.shape[0]):
                        Hsm[ti,sp,zi,:,:] = self.horizontal_smooth(H[ti,sp,zi,:,:],n=smoothing_cells)

        # Compute mean depth and volume in each pixel grid cell
        pixel_mean_depth  =  self.get_pixel_mean_depth(lon_array, lat_array)


        pixel_volume = np.zeros_like(H[0,0,:,:,:])
        for zi,zz in enumerate(z_array[:-1]):
            topotmp = -pixel_mean_depth.copy()
            topotmp[np.where(topotmp < zz)] = zz
            topotmp = z_array[zi+1] - topotmp
            topotmp[np.where(topotmp < .1)] = 0.

            pixel_volume[zi,:,:] = topotmp * pixelsize_m**2


        pixel_volume[np.where(pixel_volume==0.)] = np.nan

        # Compute mass of dry sediment in each pixel grid cell
        sed_L       = self.get_config('chemical:sediment:mixing_depth')
        sed_dens    = self.get_config('chemical:sediment:density')
        sed_poro    = self.get_config('chemical:sediment:porosity')
        pixel_sed_mass = (pixelsize_m**2 *sed_L)*(1-sed_poro)*sed_dens      # mass in kg

        # TODO this should be multiplied for the fraction og grid cell are that is not on land

        #conc = np.zeros_like(H)
        #if horizontal_smoothing:
        #    conc_sm = np.zeros_like(Hsm)
        for ti in range(H.shape[0]):
            for sp in range(self.nspecies):
                if not self.name_species[sp].lower().startswith('sed'):
                    #print('divide by volume')
                    H[ti,sp,:,:,:] = H[ti,sp,:,:,:] / pixel_volume
                    if horizontal_smoothing:
                        Hsm[ti,sp,:,:,:] = Hsm[ti,sp,:,:,:] / pixel_volume
                elif self.name_species[sp].lower().startswith('sed'):
                    #print('divide by mass')
                    #print(pixel_sed_mass)
                    H[ti,sp,:,:,:] = H[ti,sp,:,:,:] / pixel_sed_mass
                    if horizontal_smoothing:
                        Hsm[ti,sp,:,:,:] = Hsm[ti,sp,:,:,:] / pixel_sed_mass

        times = np.array( self.get_time_array()[0] )
        if time_avg_conc:
            conctmp = H[:-1,:,:,:,:]
            cshape = conctmp.shape
            mdt =    np.mean(times[1:] - times[:-1])    # output frequency in opendrift output file
            if deltat==None:
                ndt = 1
            else:
                ndt = int( deltat / (mdt.total_seconds()/3600.) )
            times2 = times[::ndt]
            times2 = times2[1:]
            odt = int(cshape[0]/ndt)
            logger.debug ('ndt '+ str(ndt))   # number of time steps over which to average in conc file
            logger.debug ('odt '+ str(odt))   # number of average slices

            Landmask=np.zeros_like(H[0:odt,:,:,:,:])
            for zi in range(len(z_array)-1):
                for sp in range(self.nspecies):
                    for ti in range(odt):
                        Landmask[ti,sp,zi,:,:] = landmask


            # This may probably be written more efficiently!
            mean_conc = np.zeros( [odt,cshape[1],cshape[2],cshape[3],cshape[4]] )
            for ii in range(odt):
                meantmp  = np.mean(conctmp[(ii*ndt):(ii+1)*ndt,:,:,:,:],axis=0)
                mean_conc[ii,:,:,:,:] = meantmp

        nc = Dataset(filename, 'w')
        nc.createDimension('x', lon_array.shape[0])
        nc.createDimension('y', lon_array.shape[1])
        nc.createDimension('depth', len(z_array)-1)
        nc.createDimension('specie', self.nspecies)
        nc.createDimension('time', H.shape[0])

#        times = self.get_time_array()[0]
        timestr = 'seconds since 1970-01-01 00:00:00'
        nc.createVariable('time', 'f8', ('time',))
        nc.variables['time'][:] = date2num(times, timestr)
        nc.variables['time'].units = timestr
        nc.variables['time'].standard_name = 'time'

        if time_avg_conc:
            nc.createDimension('avg_time', odt)
            nc.createVariable('avg_time', 'f8', ('avg_time',))
            nc.variables['avg_time'][:] = date2num(times2, timestr) # np.arange(mean_conc.shape[0])
            nc.variables['avg_time'].units = timestr

        # Projection
        nc.createVariable('projection', 'i8')
        nc.variables['projection'].proj4 = density_proj.definition_string()

        # Cell size
        nc.createVariable('cell_size','f8')
        nc.variables['cell_size'][:] = pixelsize_m
        nc.variables['cell_size'].long_name = 'Length of cell'
        nc.variables['cell_size'].unit = 'm'

        nc.createVariable('smoothing_cells','i8')
        nc.variables['smoothing_cells'][:] = smoothing_cells
        nc.variables['smoothing_cells'].long_name = 'Number of cells in each direction for horizontal smoothing'
        nc.variables['smoothing_cells'].units = '1'

        # Coordinates
        nc.createVariable('lon', 'f8', ('y','x'))
        nc.createVariable('lat', 'f8', ('y','x'))
        nc.createVariable('depth', 'f8', ('depth',))
        nc.createVariable('specie', 'i4', ('specie',))
        nc.variables['lon'][:] = lon_array.T
        nc.variables['lon'].long_name = 'longitude'
        nc.variables['lon'].short_name = 'longitude'
        nc.variables['lon'].units = 'degrees_east'
        nc.variables['lat'][:] = lat_array.T
        nc.variables['lat'].long_name = 'latitude'
        nc.variables['lat'].short_name = 'latitude'
        nc.variables['lat'].units = 'degrees_north'
        nc.variables['depth'][:] = z_array[1:]
        nc.variables['specie'][:] = np.arange(self.nspecies)
        nc.variables['specie'].long_name = ' '.join(['{}:{}'.format(isp,sp) for isp,sp in enumerate(self.name_species)])



        # Density
        #nc.createVariable('density', 'i4',
        #                  ('time','specie','depth','y', 'x'),fill_value=99999)
        #H = np.swapaxes(H, 3, 4).astype('i4')
        ##H = np.ma.masked_where(H==0, H)
        #nc.variables['density'][:] = H
        #nc.variables['density'].long_name = 'Number of elements in grid cell'
        #nc.variables['density'].grid_mapping = 'projection'
        #nc.variables['density'].units = '1'


        #if horizontal_smoothing:
        #    nc.createVariable('density_smooth', 'f8',
        #                      ('time','specie','depth','y', 'x'),fill_value=1.e36)
        #    Hsm = np.swapaxes(Hsm, 3, 4).astype('f8')
        #    #Hsm = np.ma.masked_where(Hsm==0, Hsm)
        #    nc.variables['density_smooth'][:] = Hsm
        #    nc.variables['density_smooth'].long_name = 'Horizontally smoothed number of elements in grid cell'
        #    nc.variables['density_smooth'].comment = 'Smoothed over '+str(smoothing_cells)+' grid points in all horizontal directions'



        # Chemical concentration
        if 0:
            nc.createVariable('concentration', 'f8',
                          ('time','specie','depth','y', 'x'),fill_value=1.e36)
            H = np.ma.masked_where(Landmask==1,H)
            H = np.swapaxes(H, 3, 4) #.astype('i4')
            nc.variables['concentration'][:] = H
            nc.variables['concentration'].long_name = self.get_config('chemical:compound') +' concentration ' + '\n' + 'specie '+ \
                                                            ' '.join(['{}:{}'.format(isp,sp) for isp,sp in enumerate(self.name_species)])
            nc.variables['concentration'].grid_mapping = 'projection_lonlat'
            nc.variables['concentration'].units = mass_unit+'/m3'+' (sed '+mass_unit+'/Kg)'


        # Chemical concentration, horizontally smoothed
        if horizontal_smoothing:
            nc.createVariable('concentration_smooth', 'f8',
                              ('time','specie','depth','y', 'x'),fill_value=1.e36)
            Hsm = np.ma.masked_where(Landmask==1, Hsm)
            Hsm = np.swapaxes(Hsm, 3, 4) #.astype('i4')
            nc.variables['concentration_smooth'][:] = Hsm
            nc.variables['concentration_smooth'].long_name = self.get_config('chemical:compound') +' horizontally smoothed concentration '  + '\n' + 'specie '+ \
                                                            ' '.join(['{}:{}'.format(isp,sp) for isp,sp in enumerate(self.name_species)])
            nc.variables['concentration_smooth'].grid_mapping = 'projection_lonlat'
            nc.variables['concentration_smooth'].units = mass_unit+'/m3'+' (sed '+mass_unit+'/Kg)'
            nc.variables['concentration_smooth'].comment = 'Smoothed over '+str(smoothing_cells)+' grid points in all horizontal directions'


        # Chemical concentration, time averaged
        if time_avg_conc:
            nc.createVariable('concentration_avg', 'f8',
                              ('avg_time','specie','depth','y', 'x'),fill_value=+1.e36)
            mean_conc = np.ma.masked_where(Landmask[0:odt,:,:,:,:]==1, mean_conc)
            conc2 = np.swapaxes(mean_conc, 3, 4) #.astype('i4')
            #conc2 = np.ma.masked_where(landmask==1, conc2)
            nc.variables['concentration_avg'][:] = conc2
            nc.variables['concentration_avg'].long_name = self.get_config('chemical:compound') + ' time averaged concentration ' + '\n' + 'specie '+ \
                                                            ' '.join(['{}:{}'.format(isp,sp) for isp,sp in enumerate(self.name_species)])
            nc.variables['concentration_avg'].grid_mapping = 'projection_lonlat'
            nc.variables['concentration_avg'].units = mass_unit+'/m3'+' (sed '+mass_unit+'/Kg)'


        # Volume of boxes
        nc.createVariable('volume', 'f8',
                          ('depth','y', 'x'),fill_value=0)
        pixel_volume = np.swapaxes(pixel_volume, 1, 2) #.astype('i4')
        pixel_volume = np.ma.masked_where(pixel_volume==0, pixel_volume)
        nc.variables['volume'][:] = pixel_volume
        nc.variables['volume'].long_name = 'Volume of grid cell (' + str(pixelsize_m)+'x'+str(pixelsize_m)+'m)'
        nc.variables['volume'].grid_mapping = 'projection_lonlat'
        nc.variables['volume'].units = 'm3'


        # Topography
        nc.createVariable('topo', 'f8', ('y', 'x'),fill_value=0)
        pixel_mean_depth = np.ma.masked_where(landmask==1, pixel_mean_depth)
        nc.variables['topo'][:] = pixel_mean_depth.T
        nc.variables['topo'].long_name = 'Depth of grid point'
        nc.variables['topo'].grid_mapping = 'projection_lonlat'
        nc.variables['topo'].units = 'm'

        # Binary mask
        nc.createVariable('land', 'i4', ('y', 'x'),fill_value=-1)
        #landmask = np.ma.masked_where(landmask==0, landmask)
        nc.variables['land'][:] = np.swapaxes(landmask,0,1).astype('i4')
        nc.variables['land'].long_name = 'Binary land mask'
        nc.variables['land'].grid_mapping = 'projection_lonlat'
        nc.variables['land'].units = 'm'

        nc.close()
        logger.info('Wrote to '+filename)




    def get_chemical_density_array(self, pixelsize_m, z_array,
                                       density_proj=None, llcrnrlon=None,llcrnrlat=None,
                                       urcrnrlon=None,urcrnrlat=None,
                                       weight=None, origin_marker=None):
        '''
        compute a particle concentration map from particle positions
        Use user defined projection (density_proj=<proj4_string>)
        or create a lon/lat grid (density_proj=None)
        '''
        lon = self.get_property('lon')[0]
        lat = self.get_property('lat')[0]
        times = self.get_time_array()[0]

        # Redundant ??
        if density_proj is None: # add default projection with equal-area property
            density_proj = pyproj.Proj('+proj=moll +ellps=WGS84 +lon_0=0.0')


        # create a grid in the specified projection
        x,y = density_proj(lon, lat)
        if llcrnrlon is not None:
            llcrnrx,llcrnry = density_proj(llcrnrlon,llcrnrlat)
            urcrnrx,urcrnry = density_proj(urcrnrlon,urcrnrlat)
        else:
            llcrnrx,llcrnry = x.min()-pixelsize_m, y.min()-pixelsize_m
            urcrnrx,urcrnry = x.max()+pixelsize_m, y.max()+pixelsize_m

        x_array = np.arange(llcrnrx,urcrnrx, pixelsize_m)
        y_array = np.arange(llcrnry,urcrnry, pixelsize_m)
        bins=(x_array, y_array)
        outsidex, outsidey = max(x_array)*1.5, max(y_array)*1.5
        z = self.get_property('z')[0]
        if weight is not None:
            weight_array = self.get_property(weight)[0]

        status = self.get_property('status')[0]
        specie = self.get_property('specie')[0]
        if origin_marker is not None:
            originmarker = self.get_property('origin_marker')[0]
        Nspecies = self.nspecies
        H = np.zeros((len(times),
                      Nspecies,
                      len(z_array) - 1,
                      len(x_array) - 1,
                      len(y_array) - 1
                      ))

        for sp in range(Nspecies):
            for i in range(len(times)):
                if weight is not None:
                    weights = weight_array[i,:]
                    if origin_marker is not None:
                        weight_array[i,:] = weight_array[i,:] * (originmarker[i,:]==origin_marker)
                else:
                    weights = None
                for zi in range(len(z_array)-1):
                    kktmp = ( (specie[i,:]==sp) & (z[i,:]>z_array[zi]) & (z[i,:]<=z_array[zi+1]) )
                    H[i,sp,zi,:,:], dummy, dummy = \
                        np.histogram2d(x[i,kktmp], y[i,kktmp],
                                   weights=weight_array[i,kktmp], bins=bins)

        if density_proj is not None:
            Y,X = np.meshgrid(y_array, x_array)
            lon_array, lat_array = density_proj(X,Y,inverse=True)

        return H, lon_array, lat_array





    def get_pixel_mean_depth(self,lons,lats):
        from scipy import interpolate
        # Ocean model depth and lat/lon
        h_grd = self.conc_topo
        h_grd[np.isnan(h_grd)] = 0.
        nx = h_grd.shape[0]
        ny = h_grd.shape[1]

        lat_grd = self.conc_lat[:nx,:ny]
        lon_grd = self.conc_lon[:nx,:ny]

        # Interpolate topography to new grid
        h = interpolate.griddata((lon_grd.flatten(),lat_grd.flatten()), h_grd.flatten(), (lons, lats), method='linear')

        return h



    def horizontal_smooth(self, a, n=0):
        if n==0:
            num_coarse=a
            return num_coarse


        xdm=a.shape[1]
        ydm=a.shape[0]
        #msk = self.conc_mask
        b=np.zeros([ydm+2*n,xdm+2*n],dtype=int)
        b[n:-n,n:-n]=a


        num_coarse = np.zeros([ydm,xdm],dtype=float)
        smo_tmp1=np.zeros([ydm,xdm])
        #smo_msk1=np.zeros([ydm-2*n,xdm-2*n],dtype=float)
        nlayers = 0
        for ism in np.arange(-n,n+1):
            for jsm in np.arange(-n,n+1):
                smo_tmp = b[n+jsm:ydm+n+jsm, n+ism:xdm+n+ism]
                smo_tmp1+=smo_tmp
                # Must preferrably take care of land points
#                smo_msk = msk[n+jsm:ydm-n+jsm, n+ism:xdm-n+ism]
#                smo_msk1+=smo_msk
                nlayers+=1

        if n>0:
#            num_coarse[n:-n,n:-n] = smo_tmp1 / smo_msk1
            num_coarse[:,:] = smo_tmp1 / nlayers
        else:
            num_coarse = smo_tmp1
#        num_coarse = num_coarse*msk

        return num_coarse

    def emission_factors(self, scrubber_type, chemical_compound):
        """Emission factors for heavy metals and PAHs in
            open loop and closed loop scrubbers

            Hermansson et al 2021
            https://doi.org/10.1016/j.trd.2021.102912

            bilge water, gray water, anti fouling paint,
            sewage, food waster

            from EMERGE Deliverable 2.1

            ash (atmospheric depositions)
            from EMERGE Deliverable 3.2

        """
        emission_factors_open_loop = {
            #                           mean    +/-95%
            #                           ug/L    ug/L
            "Arsenic":                  [6.8,    3.4],
            "Cadmium":                  [0.8,    0.3],
            "Chromium":                 [15.,    6.5],
            "Copper":                   [36.,    12.],
            "Iron":                     [260.,   250.],
            "Lead":                     [8.8,    4.4],
            "Mercury":                  [0.09,   0.01],
            "Nickel":                   [48.,    12.],
            "Vanadium":                 [170.,   49.],
            "Zinc":                     [110.,   59.],
            "Cobalt":                   [0.17,   0.14],
            "Selenium":                 [97.,    38],
            #
            "Naphthalene":              [2.81,   0.77],
            "Phenanthrene":             [1.51,   0.29],
            "Fluoranthene":             [0.16,   0.04],
            "Benzo-a-anthracene":       [0.12,   0.05],
            "Benzo-a-pyrene":           [0.05,   0.02],
            "Dibenzo-ah-anthracene":    [0.03,   0.01],
            #
            "Acenaphthylene":           [0.12,   0.07],
            "Acenaphthene":             [0.19,   0.07],
            "Fluorene":                 [0.46,   0.10],
            "Anthracene":               [0.08,   0.04],
            "Pyrene":                   [0.31,   0.11],
            "Chrysene":                 [0.19,   0.07],
            "Benzo-b-fluoranthene":     [0.04,   0.02],
            "Benzo-k-fluoranthene":     [0.01,   0.01],
            "Indeno-123cd-pyrene":      [0.07,   0.06],
            "Benzo-ghi-perylene":       [0.02,   0.01],
            #
            "Nitrate":                  [2830.,    2060.],
            "Nitrite":                  [760.,     680.],
            "Ammonium":                 [730.,     30.],
            "Sulphur":                  [2200000., 446000.],
            }

        emission_factors_closed_loop = {
            #                           mean    +/-95%
            #                           ug/L    ug/L
            "Arsenic":                  [22.,    9.4],
            "Cadmium":                  [0.55,   0.19],
            "Chromium":                 [1300.,  1700.],
            "Copper":                   [480.,   230.],
            "Iron":                     [490.,   82.],
            "Lead":                     [7.7,    3.1],
            "Mercury":                  [0.07,   0.02],
            "Nickel":                   [2700.,  860.],
            "Vanadium":                 [9100.,  3200.],
            "Zinc":                     [370.,   200.],
            "Cobalt":                   [0.,     0.],
            "Selenium":                 [0.,     0.],
            #
            "Naphthalene":              [2.08,   1.05],
            "Phenanthrene":             [5.00,   2.30],
            "Fluoranthene":             [0.63,	 0.41],
            "Benzo-a-anthracene":       [0.30,	 0.29],
            "Benzo-a-pyrene":           [0.06,	 0.05],
            "Dibenzo-ah-anthracene":    [0.03,	 0.02],
            #
            "Acenaphthylene":           [0.09,   0.06],
            "Acenaphthene":             [0.47,   0.31],
            "Fluorene":                 [1.32,   0.54],
            "Anthracene":               [1.55,   2.00],
            "Pyrene":                   [0.76,   0.59],
            "Chrysene":                 [0.50,   0.45],
            "Benzo-b-fluoranthene":     [0.14,   0.12],
            "Benzo-k-fluoranthene":     [0.02,   0.02],
            "Indeno-123-cd-pyrene":     [0.04,   0.03],
            "Benzo-ghi-perylene":       [0.07,   0.07],
            #
            "Nitrate":                  [110980.,   100000.],
            "Nitrite":                  [55760.,    55000.],
            "Ammonium":                 [0.,        0.],
            "Sulphur":                  [12280000., 10104000.],
            }

        emission_factors_grey_water = {
            #                           mean    +/-95%
            #                           ug/L    ug/L
            "Arsenic":                  [5.98,    3.17],
            "Cadmium":                  [0.16,    0.09],
            "Chromium":                 [7.28,    2.06],
            "Copper":                   [267.,    97.],
            "Lead":                     [25.6,    21.01],
            "Mercury":                  [0.16,    0.09],
            "Nickel":                   [25.0,    19.36],
            "Selenium":                 [16.1,    10.64],
            "Zinc":                     [517.,    112.],
            #
            "Nitrogen":                  [28900., 0.0],
         }

        emission_factors_bilge_water = {
            #                           mean    +/-95%
            #                           ug/L    ug/L
            "Arsenic":                  [35.9,    33.2],
            "Cadmium":                  [0.32,   0.07],
            "Chromium":                 [16.3,   15.4],
            "Copper":                   [49.7,   22.9],
            "Lead":                     [3.0,    1.24],
            "Nickel":                   [71.1,   11.8],
            "Selenium":                 [2.95,     1.01],
            "Vanadium":                 [76.5,   22.4],
            "Zinc":                     [949.,   660.],
            #
            "Nitrate":                  [110980.,  100000.],
            "Nitrite":                  [55760.,   55000.],
            "Ammonium":                 [0.,     0.],
            "Sulphur":                  [12280000., 10104000.],
            #
            "Naphthalene":              [50.6,   34.3],
            "Phenanthrene":             [3.67,   2.51],
            "Fluoranthene":             [0.60,   0.96],
            "Benzo(a)anthracene":       [0.10,   0.18],
            "Benzo(a)pyrene":           [0.10,   0.15],
            "Dibenzo(a,h)anthracene":   [0.02,   0.01],
            #
            "Acenaphthylene":           [0.29,   0.17],
            "Acenaphthene":             [1.42,   0.86],
            "Fluorene":                 [3.33,   2.43],
            "Anthracene":               [0.22,   0.14],
            "Pyrene":                   [1.23,   1.33],
            "Chrysene":                 [0.17,   0.25],
            "Benzo(b)fluoranthene":     [0.09,   0.13],
            "Benzo(k)fluoranthene":     [0.03,   0.00],
            "Indeno(1,2,3-cd)pyrene":   [0.05,   0.06],
            "Benzo(ghi)perylene":       [0.13,   0.16],
         }

        emission_factors_sewage_water = {
            #                           mean    +/-95%
            #                           ug/L    ug/L
            "Arsenic":                  [22.9,    7.4],
            "Cadmium":                  [0.12,   0.10],
            "Chromium":                 [11.9,    8.2],
            "Copper":                   [319,     190],
            "Lead":                     [6.5,     3.1],
            "Mercury":                  [0.22,   0.12],
            "Nickel":                   [32.3,   21.3],
            "Selenium":                 [43.7,   18.3],
            "Zinc":                     [395.,   174.],
            #
            "Nitrogen":                  [430.,  0.],
         }

        emission_factors_AFP = {
            # Copper = 63.546 g/mol
            # Zinc = 65.38 g/mol
            # CuPyr = 315.86 g/mol = Copper(II) pyrithione = 0.2112 of Cu
            # CuO = 79.55 g/mol = Copper(II) oxide = 0.7989 of Cu
            # Zineb = 275.7 g/mol = Zinc ethylenebis(dithiocarbamate) = 0.2371 of Zn
            # ZnO = 81.38 g/mol = Zinc(II) oxide = 0.8033 of Zn
            # ZPyr = 317.70 g/mol = Zinc(II) pyrithione = 0.2058 of Zn

            #                           mean    +/-95%
            #                           ug/L    ug/L
            "CuO_AFP":                  [0.7989,    0.],
            "CuPyr_AFP":                [0.2112,    0.],
            "Zineb_AFP":                [0.2371,    0.],
            "ZnO_AFP":                  [0.8033,    0.],
            "ZnPyr_AFP":                [0.2058,    0.],
         }

        emission_factors_SILAM_ash = {
            #                           g/g
            "Aresenic":                 [8.09E-5],
            "Cadmium":                  [6.30E-6],
            "Chromium":                 [2.10E-4],
            "Copper":                   [2.52E-4],
            "Iron":                     [2.52E-2],
            "Mercury":                  [6.30E-6],
            "Nickel":                   [4.10E-2],
            "Lead":                     [1.16E-4],
            "Vanadium":                 [8.30E-2],
            "Zinc":                     [2.42E-3],
         }

        if scrubber_type=="open_loop":
            Emission_factors = emission_factors_open_loop.get(chemical_compound)[0]
        elif scrubber_type=="closed_loop":
            Emission_factors = emission_factors_closed_loop.get(chemical_compound)[0]
        elif scrubber_type=="bilge_water":
            Emission_factors = emission_factors_bilge_water.get(chemical_compound)[0]
        elif scrubber_type=="grey_water":
            Emission_factors = emission_factors_grey_water.get(chemical_compound)[0]
        elif scrubber_type=="sewage_water":
            Emission_factors = emission_factors_sewage_water.get(chemical_compound)[0]
        elif scrubber_type=="AFP": # Copper and Zinc from antifouling paint
            Emission_factors = 1e6*emission_factors_AFP.get(chemical_compound)[0]  # 1g = 1e6 ug: AFP is expressed as g
        elif scrubber_type=="AFP_metals_total":
            Emission_factors = 1e6 # g to ug
        elif scrubber_type=="N_sewage": # Nitrogen from sewage
            Emission_factors = 1e9  # 1kg = 1e9 ug: N_sewage is expressed as kg
        elif scrubber_type=="N_foodwaste": # Nitrogen from foodwaste
            Emission_factors = 1e9  # 1kg = 1e9 ug: N_sewage is expressed as kg
        elif scrubber_type=="SILAM_metals":
            Emission_factors = 1e9  #+ 1kg = 1e9 ug: Lead and Cadmium depositions given in kg
        elif scrubber_type=="SILAM_metals_from_ash":
            Emission_factors = 1e9*emission_factors_SILAM_ash.get(chemical_compound)[0] # 1kg=1e9ug: Ash depositions given in kg

        return Emission_factors
        # TODO: Add emission uncertainty based on 95% confidence interval

    def seed_from_DataArray(self, steam, lowerbound=0, higherbound=np.inf, radius=0, scrubber_type="open_loop", chemical_compound="Copper", mass_element_ug=100e3, number_of_elements=None, **kwargs):
            """Seed elements based on a dataarray with STEAM emission data

            Arguments:
                steam: dataarray with steam emission data, with coordinates
                    * latitude   (latitude) float32
                    * longitude  (longitude) float32
                    * time       (time) datetime64[ns]


                radius:      scalar, unit: meters
                lowerbound:  scalar, elements with lower values are discarded
            """

            if chemical_compound is None:
                chemical_compound = self.get_config('chemical:compound')

            #mass_element_ug=1e3      # 1e3 - 1 element is 1mg chemical
            #mass_element_ug=20e3      # 100e3 - 1 element is 100mg chemical
            #mass_element_ug=100e3      # 100e3 - 1 element is 100mg chemical
            #mass_element_ug=1e6     # 1e6 - 1 element is 1g chemical

            sel=np.where((steam > lowerbound) & (steam < higherbound))
            t=steam.time[sel[0]].data
            la=steam.latitude[sel[1]].data
            lo=steam.longitude[sel[2]].data

            data=np.array(steam.data)

            if number_of_elements is not None:
                total_volume = np.sum(data[sel])
                total_mass = total_volume * self.emission_factors(scrubber_type, chemical_compound)
                mass_element_ug = total_mass / number_of_elements
                mass_element_ug_0 = total_mass / number_of_elements

            for i in range(0,t.size):
                scrubberwater_vol_l=data[sel][i]
                mass_ug=scrubberwater_vol_l * self.emission_factors(scrubber_type, chemical_compound)

                if number_of_elements is None:
                    number=np.array(mass_ug / mass_element_ug).astype('int')
                else:
                    number=np.ceil(np.array(mass_ug / mass_element_ug_0)).astype('int')
                    mass_element_ug=mass_ug/number

                time = datetime.utcfromtimestamp((t[i] - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's'))

                if number>0:
                    z = -1*np.random.uniform(0, 1, number)
                    self.seed_elements(lon=lo[i]*np.ones(number), lat=la[i]*np.ones(number),
                                radius=radius, number=number, time=time,
                                mass=mass_element_ug,mass_degraded=0,mass_volatilized=0, z=z, origin_marker=1)

                mass_residual = mass_ug - number*mass_element_ug

                if mass_residual>0 and number_of_elements is None:
                    z = -1*np.random.uniform(0, 1, 1)
                    self.seed_elements(lon=lo[i], lat=la[i],
                                radius=radius, number=1, time=time,
                                mass=mass_residual,mass_degraded=0,mass_volatilized=0, z=z, origin_marker=1)

    seed_from_STEAM = seed_from_DataArray
    ''' Alias of seed_from_DataArray method for backward compatibility
    '''

    def init_chemical_compound(self, chemical_compound = None):
        ''' Chemical parameters for a selection of PAHs:
            Naphthalene, Phenanthrene, Fluorene,
            Benzo-a-anthracene, Benzo-a-pyrene, Dibenzo-ah-anthracene

            Data collected from literature by
            Isabel Hanstein (University of Heidelberg / Norwegian Meteorological Insitute)
            Mattia Boscherini, Loris Calgaro (University Ca' Foscari, Venezia)
            Manuel Aghito (Norwegian Meteorological Institute / University of Bergen)
        '''

        if chemical_compound is not None:
            self.set_config('chemical:compound',chemical_compound)

        if self.get_config('chemical:compound') is None:
            raise ValueError("Chemical compound not defined")

        if  self.get_config('chemical:compound') in ["Naphthalene",'C1-Naphthalene','Acenaphthene','Acenaphthylene','Fluorene']:

            #partitioning
            self.set_config('chemical:transfer_setup','organics')
            self.set_config('chemical:transformations:dissociation','nondiss')
            self.set_config('chemical:transformations:LogKOW',3.361)
            self.set_config('chemical:transformations:TrefKOW',25)
            self.set_config('chemical:transformations:DeltaH_KOC_Sed',-21036)
            self.set_config('chemical:transformations:DeltaH_KOC_DOM',-25900)                 ### Phenanthrene value
            self.set_config('chemical:transformations:Setchenow', 0.2503)

            #degradation
            self.set_config('chemical:transformations:t12_W_tot', 224.08)
            self.set_config('chemical:transformations:Tref_kWt', 25)
            self.set_config('chemical:transformations:DeltaH_kWt', 50000)                     ### generic
            self.set_config('chemical:transformations:t12_S_tot', 5012.4)
            self.set_config('chemical:transformations:Tref_kSt', 25)
            self.set_config('chemical:transformations:DeltaH_kSt', 50000)                     ### generic

            #volatilization
            self.set_config('chemical:transformations:MolWt', 128.1705)
            self.set_config('chemical:transformations:Henry', 4.551e-4)

            self.set_config('chemical:transformations:Vpress', 11.2)
            self.set_config('chemical:transformations:Tref_Vpress', 25)
            self.set_config('chemical:transformations:DeltaH_Vpress', 55925)

            self.set_config('chemical:transformations:Solub', 31.4)
            self.set_config('chemical:transformations:Tref_Solub', 25)
            self.set_config('chemical:transformations:DeltaH_Solub', 25300)

        elif self.get_config('chemical:compound') in ["Phenanthrene",'Dibenzothiophene','C2-Naphthalene','Anthracene','C3-Naphthalene','C1-Dibenzothiophene']:

            #partitioning
            self.set_config('chemical:transfer_setup','organics')
            self.set_config('chemical:transformations:dissociation','nondiss')
            self.set_config('chemical:transformations:LogKOW',4.505)
            self.set_config('chemical:transformations:TrefKOW',25)
            self.set_config('chemical:transformations:DeltaH_KOC_Sed',-24900)
            self.set_config('chemical:transformations:DeltaH_KOC_DOM',-25900)
            self.set_config('chemical:transformations:Setchenow', 0.3026)

            #degradation
            self.set_config('chemical:transformations:t12_W_tot', 1125.79)
            self.set_config('chemical:transformations:Tref_kWt', 25)
            self.set_config('chemical:transformations:DeltaH_kWt', 50000)                     ### generic
            self.set_config('chemical:transformations:t12_S_tot', 29124.96)
            self.set_config('chemical:transformations:Tref_kSt', 25)
            self.set_config('chemical:transformations:DeltaH_kSt', 50000)                     ### generic

            #volatilization
            self.set_config('chemical:transformations:MolWt', 178.226)
            self.set_config('chemical:transformations:Henry', 4.294e-5)

            self.set_config('chemical:transformations:Vpress', 0.0222)
            self.set_config('chemical:transformations:Tref_Vpress', 25)
            self.set_config('chemical:transformations:DeltaH_Vpress', 71733)

            self.set_config('chemical:transformations:Solub', 1.09)
            self.set_config('chemical:transformations:Tref_Solub', 25)
            self.set_config('chemical:transformations:DeltaH_Solub', 34800)

        elif self.get_config('chemical:compound') in ["Fluoranthene",'Pyrene','C1-Phenanthrene','C2-Dibenzothiophene']:

            #partitioning
            self.set_config('chemical:transfer_setup','organics')
            self.set_config('chemical:transformations:dissociation','nondiss')
            self.set_config('chemical:transformations:LogKOW',5.089)
            self.set_config('chemical:transformations:TrefKOW',25)
            self.set_config('chemical:transformations:DeltaH_KOC_Sed',-47413)
            self.set_config('chemical:transformations:DeltaH_KOC_DOM',-27900)
            self.set_config('chemical:transformations:Setchenow', 0.2885)

            #degradation
            self.set_config('chemical:transformations:t12_W_tot', 1705.02)
            self.set_config('chemical:transformations:Tref_kWt', 25)
            self.set_config('chemical:transformations:DeltaH_kWt', 50000)                     ### generic
            self.set_config('chemical:transformations:t12_S_tot', 55000)
            self.set_config('chemical:transformations:Tref_kSt', 25)
            self.set_config('chemical:transformations:DeltaH_kSt', 50000)                     ### generic

            #volatilization
            self.set_config('chemical:transformations:MolWt', 202.26)
            self.set_config('chemical:transformations:Henry', 1.439e-5)

            self.set_config('chemical:transformations:Vpress', 0.00167)
            self.set_config('chemical:transformations:Tref_Vpress', 25)
            self.set_config('chemical:transformations:DeltaH_Vpress', 79581)

            self.set_config('chemical:transformations:Solub', 0.231)
            self.set_config('chemical:transformations:Tref_Solub', 25)
            self.set_config('chemical:transformations:DeltaH_Solub', 30315)

        elif self.get_config('chemical:compound') in ["Benzo-a-anthracene",'C2-Phenanthrene','Benzo-b-fluoranthene','Chrysene']:

            #partitioning
            self.set_config('chemical:transfer_setup','organics')
            self.set_config('chemical:transformations:dissociation','nondiss')
            self.set_config('chemical:transformations:LogKOW',5.724)
            self.set_config('chemical:transformations:TrefKOW',25)
            self.set_config('chemical:transformations:DeltaH_KOC_Sed', -38000)                ### Pyrene value
            self.set_config('chemical:transformations:DeltaH_KOC_DOM', -25400)                ### Pyrene value
            self.set_config('chemical:transformations:Setchenow', 0.3605)

            #degradation
            self.set_config('chemical:transformations:t12_W_tot', 1467.62)
            self.set_config('chemical:transformations:Tref_kWt', 25)
            self.set_config('chemical:transformations:DeltaH_kWt', 50000)                     ### generic
            self.set_config('chemical:transformations:t12_S_tot', 46600)
            self.set_config('chemical:transformations:Tref_kSt', 25)
            self.set_config('chemical:transformations:DeltaH_kSt', 50000)                     ### generic

            #volatilization
            self.set_config('chemical:transformations:MolWt', 228.29)
            self.set_config('chemical:transformations:Henry', 6.149e-6)

            self.set_config('chemical:transformations:Vpress', 0.0000204)
            self.set_config('chemical:transformations:Tref_Vpress', 25)
            self.set_config('chemical:transformations:DeltaH_Vpress', 100680)

            self.set_config('chemical:transformations:Solub', 0.011)
            self.set_config('chemical:transformations:Tref_Solub', 25)
            self.set_config('chemical:transformations:DeltaH_Solub', 46200)

        elif self.get_config('chemical:compound') in ["Benzo-a-pyrene",'C3-Dibenzothiophene','C3-Phenanthrene']:

            #partitioning
            self.set_config('chemical:transfer_setup','organics')
            self.set_config('chemical:transformations:dissociation','nondiss')
            self.set_config('chemical:transformations:LogKOW', 6.124)
            self.set_config('chemical:transformations:TrefKOW',25)
            self.set_config('chemical:transformations:DeltaH_KOC_Sed', -43700)                ### mean value 16 PAHs
            self.set_config('chemical:transformations:DeltaH_KOC_DOM', -31280)
            self.set_config('chemical:transformations:Setchenow', 0.171)

            #degradation
            self.set_config('chemical:transformations:t12_W_tot', 1491.42)
            self.set_config('chemical:transformations:Tref_kWt', 25)
            self.set_config('chemical:transformations:DeltaH_kWt', 50000)                     ### generic
            self.set_config('chemical:transformations:t12_S_tot', 44934.76)
            self.set_config('chemical:transformations:Tref_kSt', 25)
            self.set_config('chemical:transformations:DeltaH_kSt', 50000)                     ### generic

            #volatilization
            self.set_config('chemical:transformations:MolWt', 252.32)
            self.set_config('chemical:transformations:Henry', 6.634e-7)

            self.set_config('chemical:transformations:Vpress', 0.00000136)
            self.set_config('chemical:transformations:Tref_Vpress', 25)
            self.set_config('chemical:transformations:DeltaH_Vpress', 107887)

            self.set_config('chemical:transformations:Solub', 0.00229)
            self.set_config('chemical:transformations:Tref_Solub', 25)
            self.set_config('chemical:transformations:DeltaH_Solub', 38000)

        elif self.get_config('chemical:compound') in ["Dibenzo-ah-anthracene",'Benzo-k-fluoranthene','Benzo-ghi-perylene','Indeno-123cd-pyrene']:

            #partitioning
            self.set_config('chemical:transfer_setup','organics')
            self.set_config('chemical:transformations:dissociation','nondiss')
            self.set_config('chemical:transformations:LogKOW', 6.618)
            self.set_config('chemical:transformations:TrefKOW',25)
            self.set_config('chemical:transformations:DeltaH_KOC_Sed', -43700)                ### mean value 16 PAHs
            self.set_config('chemical:transformations:DeltaH_KOC_DOM', -30900)
            self.set_config('chemical:transformations:Setchenow', 0.338)

            #degradation
            self.set_config('chemical:transformations:t12_W_tot', 1464.67)
            self.set_config('chemical:transformations:Tref_kWt', 25)
            self.set_config('chemical:transformations:DeltaH_kWt', 50000)                     ### generic
            self.set_config('chemical:transformations:t12_S_tot', 40890.08)
            self.set_config('chemical:transformations:Tref_kSt', 25)
            self.set_config('chemical:transformations:DeltaH_kSt', 50000)                     ### generic

            #volatilization
            self.set_config('chemical:transformations:MolWt', 278.35)
            self.set_config('chemical:transformations:Henry', 4.894e-8)

            self.set_config('chemical:transformations:Vpress', 0.0000000427)
            self.set_config('chemical:transformations:Tref_Vpress', 25)
            self.set_config('chemical:transformations:DeltaH_Vpress', 112220)

            self.set_config('chemical:transformations:Solub', 0.00142)
            self.set_config('chemical:transformations:Tref_Solub', 25)
            self.set_config('chemical:transformations:DeltaH_Solub', 38000)                   ### Benzo-a-pyrene value

        elif self.get_config('chemical:compound') == "Copper":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 60.1)            # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 50)             # Merlin Expo, high confidence
            self.set_config('chemical:transformations:S0', 17.0)            # note below

        elif self.get_config('chemical:compound') == "Zinc":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 173)             # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 100)            # Merlin Expo, high confidence
            self.set_config('chemical:transformations:S0', 17.0)            # note below

        elif self.get_config('chemical:compound') == "Lead":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 369)             # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 500)            # Merlin Expo, strong confidence
            self.set_config('chemical:transformations:S0', 17.0)            # note below

        elif self.get_config('chemical:compound') == "Vanadium":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 42.9)            # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 5)              # Merlin Expo, weak confidence
            self.set_config('chemical:transformations:S0', 17.0)            # note below

        elif self.get_config('chemical:compound') == "Cadmium":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 134)             # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 79)             # Merlin Expo, strong confidence
            #self.set_config('chemical:transformations:Kd', 6.6)            # Turner Millward 2002
            self.set_config('chemical:transformations:S0', 17.0)            # note below

        elif self.get_config('chemical:compound') == "Chromium":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 124)             # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 130)            # Cr(III) Merlin Expo, moderate confidence
            #self.set_config('chemical:transformations:Kd', 180)            # Turner Millward 2002
            self.set_config('chemical:transformations:S0', 17.0)            # note below

        elif self.get_config('chemical:compound') == "Nickel":
            self.set_config('chemical:transfer_setup','metals')
            self.set_config('chemical:transformations:Kd', 31.1)            # Tomczak et Al 2019
            #self.set_config('chemical:transformations:Kd', 25)             # Merlin Expo, strong confidence
            #self.set_config('chemical:transformations:Kd', 5.3)            # Turner Millward 2002
            self.set_config('chemical:transformations:S0', 17.0)            # note below

# Default value for S0 is set to 17.0. This correspond to a Kd at salinity 35 being 32.7%
# of the fresh water value, which was the average reduction obtained comparing the values
# in Tomczak et Al 2019 to the "ocean margins" recommended values in IAEA TRS no.422, for a
# selection of metals (Cd, Cr, Hg, Ni, Pb, Zn). This gives very similar results the value
# 15.8, suggested in Perianez 2018.
# https://doi.org/10.1016/j.apgeochem.2019.04.003
# https://www-pub.iaea.org/MTCD/Publications/PDF/TRS422_web.pdf
# https://doi.org/10.1016/j.jenvrad.2018.02.014
#
# Merlin Expo Kd values are mean values from Allison and Allison 2005
# https://cfpub.epa.gov/si/si_public_record_report.cfm?dirEntryId=135783

    def plot_mass(self,
                  legend=['dissolved','SPM','sediment'],
                  mass_unit='g',
                  time_unit='hours',
                  title=None,
                  filename=None,
                  start_date=None):
        """Plot chemical mass distribution between the different species
            legend      list of specie labels, for example ['dissolved','SPM','sediment']
            mass_unit   'g','mg','ug'
            time_unit   'seconds', 'minutes', 'hours' , 'days'
            title       figure title string
        """

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        if not title == []:
            fig.suptitle(title)

        mass=self.get_property('mass')
        sp=self.get_property('specie')

        steps=self.steps_output

        bars=np.zeros((steps,5))

        mass_conversion_factor=1e-6
        if mass_unit=='g' and self.elements.variables['mass']['units']=='ug':
            mass_conversion_factor=1e-6
        if mass_unit=='mg' and self.elements.variables['mass']['units']=='ug':
            mass_conversion_factor=1e-3
        if mass_unit=='ug' and self.elements.variables['mass']['units']=='ug':
            mass_conversion_factor=1
        if mass_unit=='kg' and self.elements.variables['mass']['units']=='ug':
            mass_conversion_factor=1e-9

        time_conversion_factor = self.time_step_output.total_seconds() / (60*60)
        if time_unit=='seconds':
            time_conversion_factor = self.time_step_output.total_seconds()
        if time_unit=='minutes':
            time_conversion_factor = self.time_step_output.total_seconds() / 60
        if time_unit=='hours':
            time_conversion_factor = self.time_step_output.total_seconds() / (60*60)
        if time_unit=='days':
            time_conversion_factor = self.time_step_output.total_seconds() / (24*60*60)

        for i in range(steps):

            bars[i]=[np.sum(mass[0][i,:]*(sp[0][i,:]==0))*mass_conversion_factor,
                     np.sum(mass[0][i,:]*(sp[0][i,:]==1))*mass_conversion_factor,
                     np.sum(mass[0][i,:]*(sp[0][i,:]==2))*mass_conversion_factor,
                     np.sum(mass[0][i,:]*(sp[0][i,:]==3))*mass_conversion_factor,
                     np.sum(mass[0][i,:]*(sp[0][i,:]==4))*mass_conversion_factor]
        bottom=np.zeros_like(bars[:,0])
        if 'dissolved' in legend:
            ax.bar(np.arange(steps),bars[:,self.num_lmm],width=1.25,color='midnightblue')
            bottom=bars[:,self.num_lmm]
            print('dissolved' + ' : ' + str(bars[-1,self.num_lmm]) + mass_unit +' ('+ str(100*bars[-1,self.num_lmm]/np.sum(bars[-1,:]))+'%)')
        if 'DOC' in legend:
            ax.bar(np.arange(steps),bars[:,self.num_humcol],bottom=bottom,width=1.25,color='royalblue')
            bottom=bottom+bars[:,self.num_humcol]
            print('DOC' + ' : ' + str(bars[-1,self.num_humcol]) + mass_unit +' ('+ str(100*bars[-1,self.num_humcol]/np.sum(bars[-1,:]))+'%)')
        if 'SPM' in legend:
            ax.bar(np.arange(steps),bars[:,self.num_prev],bottom=bottom,width=1.25,color='palegreen')
            bottom=bottom+bars[:,self.num_prev]
            print('SPM' + ' : ' + str(bars[-1,self.num_prev]) + mass_unit +' ('+ str(100*bars[-1,self.num_prev]/np.sum(bars[-1,:]))+'%)')
        if 'sediment' in legend:
            ax.bar(np.arange(steps),bars[:,self.num_srev],bottom=bottom,width=1.25,color='orange')
            bottom=bottom+bars[:,self.num_srev]
            print('sediment' + ' : ' + str(bars[-1,self.num_srev]) + mass_unit +' ('+ str(100*bars[-1,self.num_srev]/np.sum(bars[-1,:]))+'%)')

        ax.legend(list(filter(None, legend)))
        ax.set_ylabel('mass (' + mass_unit + ')')
        if start_date is None:
            ax.axes.get_xaxis().set_ticklabels(np.round(ax.axes.get_xticks() * time_conversion_factor))
            ax.set_xlabel('time (' + time_unit + ')')
        else:
            date_values = [datetime.strptime(start_date,"%Y-%m-%d") + i*self.time_step for i in range(steps)]

            # Set a fraction of datetime values as labels for the x-axis
            fraction = 24*7  # Show every second datetime value
            ax.set_xticks(np.arange(0, steps, fraction))
            ax.set_xticklabels([date.strftime('%m-%d') for date in date_values[::fraction]])
            ax.set_xlabel('time (month-day)')
        fig.show()

        if filename is not None:
            plt.savefig(filename, format=filename[-3:], transparent=True, bbox_inches="tight", dpi=300)
