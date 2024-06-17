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
# Copyright 2015, Magne Simonsen, MET Norway

import numpy as np
import logging; logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED
import pyproj

# Defining the radionuclide element properties
class Radionuclide(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for radionuclides
    """

    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.}),
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': '[]',
                                       'default': 31.25}),  # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 2650.}),  # Mineral particles
        ('specie', {'dtype': np.int32,
                    'units': '',
                    'default': 0}),
#         ('transfer_rates1D', {'dtype':np.array(3, dtype=np.float32),
#                     'units': '1/s',
#                     'default': 0.})
        # ('LMM_fraction', {'dtype':np.float32,
        #                   'units':'',
        #                   'default':0,
        #                   'seed':False}),
        # ('particle_fraction', {'dtype':np.float32,
        #                   'units':'',
        #                   'default':0,
        #                   'seed':False})
        ])


class RadionuclideDrift(OceanDrift):
    """Radionuclide particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Radionuclide functionality include interactions with solid matter
        (particles and sediments) through transformation processes, implemented
        with stochastic approach for speciation.

        Under construction.
    """

    ElementType = Radionuclide

    required_variables = {
        'x_sea_water_velocity': {'fallback': None},
        'y_sea_water_velocity': {'fallback': None},
        'sea_surface_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': None},
        'ocean_vertical_diffusivity': {'fallback': 0.0001, 'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
#        'surface_downward_x_stress': {'fallback': 0},
#        'surface_downward_y_stress': {'fallback': 0},
#        'turbulent_kinetic_energy': {'fallback': 0},
#        'turbulent_generic_length_scale': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'conc3': {'fallback': 1.e-3},
        }


    def specie_num2name(self,num):
        return self.name_species[num]

    def specie_name2num(self,name):
        num = self.name_species.index(name)
        return num

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(RadionuclideDrift, self).__init__(*args, **kwargs)

        # TODO: descriptions and units must be added in config setting below
        self._add_config({
            'radionuclide:dissolved_diameter': {'type': 'float', 'default': 0,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:particle_diameter': {'type': 'float', 'default': 5e-6,
                'min': .45e-6, 'max': 63.e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'Mean particle diameter. Determines the settling velocity.'},
            'radionuclide:particle_diameter_uncertainty': {'type': 'float', 'default': 1e-7,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'Standard deviation of particle size distribution.'},
            'radionuclide:particle_diameter_minimum': {'type': 'float', 'default': 0.45e-6,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'Mimimum particle size.'},
            'radionuclide:particle_diameter_maximum': {'type': 'float', 'default': 63.e-6,
                'min': 0, 'max': 100e-6, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': 'Maximum particle size.'},
            'radionuclide:particlesize_distribution': {'type':'enum',
                'enum': ['normal','lognormal'], 'default':'normal',
                'level':CONFIG_LEVEL_ADVANCED,
                'description':'Distribution of particle diameter around a mean value at seeding, NB: not at sorption!'},
            'seed:LMM_fraction': {'type': 'float','default': .1,
                'min': 0, 'max': 1, 'units': '1',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Fraction of initial discharge released as LMM species'},
            'seed:particle_fraction': {'type': 'float','default': 0.9,
                'min': 0, 'max': 1, 'units': '1',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Fraction of initial discharge released as particle species'},
            'seed:slowly_fraction': {'type': 'float','default': 0.,
                'min': 0, 'max': 1, 'units': '1',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Fraction of PARTICLE discharge released as slowly reversible particle species'},
            'seed:total_release': {'type': 'float','default': 100.e9,
                'min': 0, 'max': 1e36, 'units': 'Bq',
                'level': CONFIG_LEVEL_ESSENTIAL, 'description': 'Total release (Bq)'},

            # ISOTOPES
            'radionuclide:isotope': {'type': 'enum', 'default': '137Cs',
                                   'enum': ['Al', '137Cs', '129I', '241Am', 'manual'],
                                   'level':CONFIG_LEVEL_ESSENTIAL, 'description':'Isotope'},
            'radionuclide:specie_setup': {'type': 'enum', 'default': 'LMM + Rev',
                                   'enum': ['LMM + Rev', 'LMM + Rev + Slow rev', 'LMM + Rev + Irrev', 
                                            'LMM + Rev + Slow rev + Irrev', 'LMM + Colloid + Rev'],
                                   'level':CONFIG_LEVEL_ESSENTIAL, 'description':'Species enabled'},


            # Transformation parameters
            'radionuclide:transformations:Kd': {'type': 'float', 'default': 2.0,
                'min': 0, 'max': 1e9, 'units': 'm3/kg',
                'level': CONFIG_LEVEL_BASIC, 'description': ''},
            'radionuclide:transformations:Dc': {'type': 'float', 'default': 1.16e-5,
                'min': 0, 'max': 1e6, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:transformations:slow_coeff': {'type': 'float', 'default': 1.2e-7,
                'min': 0, 'max': 1e6, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            # Sediment
            'radionuclide:sediment:sedmixdepth': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:sediment_density': {'type': 'float', 'default': 2600,
                'min': 0, 'max': 10000, 'units': 'kg/m3',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:effective_fraction': {'type': 'float', 'default': 0.9,
                'min': 0, 'max': 1, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:corr_factor': {'type': 'float', 'default': 0.1,
                'min': 0, 'max': 10, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:porosity': {'type': 'float', 'default': 0.6,
                'min': 0, 'max': 1, 'units': '',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:layer_thick': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:desorption_depth': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:desorption_depth_uncert': {'type': 'float', 'default': .5,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:resuspension_depth': {'type': 'float', 'default': 1,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:resuspension_depth_uncert': {'type': 'float', 'default': .5,
                'min': 0, 'max': 100, 'units': 'm',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            'radionuclide:sediment:resuspension_critvel': {'type': 'float', 'default': .01,
                'min': 0, 'max': 1, 'units': 'm/s',
                'level': CONFIG_LEVEL_ADVANCED, 'description': ''},
            # OUTPUT config
            'radionuclide:output:depthintervals': {'type':'str',
                                                   'default': '-25, -10., -5., -1.',
                                                   'min_length':0, 'max_length':60,
                'level':CONFIG_LEVEL_ESSENTIAL,
                'description':'Depth intervals for computation of concentration'},
                        })





    def prepare_run(self):

        logger.info( 'Number of species: {}'.format(self.nspecies) )
        for i,sp in enumerate(self.name_species):
            logger.info( '{:>3} {}'.format( i, sp ) )


        logger.info( 'Isotope: %s' % self.get_config('radionuclide:isotope'))
        logger.info( 'specie setup: %s' % self.get_config('radionuclide:specie_setup'))

        logger.info('nspecies: %s' % self.nspecies)
        logger.info('Transfer rates:\n %s' % self.transfer_rates)

        tmp = self.get_config('radionuclide:output:depthintervals')
        self.depthintervals = [float(item) for item in tmp.split(',')]
        logger.info('DEPTH_INTERVALS: %s',self.depthintervals)



        super(RadionuclideDrift, self).prepare_run()


    def init_species(self):
        
        # Initialize slowly and irrev fractions to False, set True later if enabled
        self.species_slowly_fraction = False
        self.species_irreversible_fraction = False


        # Initialize species, add enabled species to name_species list
        self.name_species=[]
        self.specie_setup = self.get_config('radionuclide:specie_setup')

        if 'lmm + rev' in self.specie_setup.casefold():
            self.name_species.append('LMM')
            self.name_species.append('Particle reversible')
            self.name_species.append('Sediment reversible')
        # else:
        #     if not self.isotope == 'Al':
        #         logger.error('No valid specie setup {}'.format(self.get_config('radionuclide:specie_setup')))
        #     else:
        #         pass


        if  '+ slow rev' in self.specie_setup.casefold():
            self.name_species.append('Particle slowly reversible')
            self.name_species.append('Sediment slowly reversible')
            self.species_slowly_fraction = True
        if '+ irrev' in self.specie_setup.casefold():
            self.name_species.append('Particle irreversible')
            self.name_species.append('Sediment irreversible')
            self.species_irreversible_fraction = True
        if self.specie_setup.casefold() == 'lmm + colloid + rev':
            self.name_species.append('LMMcation')
            self.name_species.append('LMManion')
            self.name_species.append('Humic colloid')
            self.name_species.append('Polymer')
            self.name_species.append('Particle reversible')
            self.name_species.append('Sediment reversible')


        self.nspecies      = len(self.name_species)
        logger.info( 'Number of species: {}'.format(self.nspecies) )
        for i,sp in enumerate(self.name_species):
            logger.info( '{:>3} {}'.format( i, sp ) )





    def set_init_diameter(self, num, idxs,diam):
        """ Initialize diameter for particles, according to size distribution """

        distr = self.get_config('radionuclide:particlesize_distribution')

        init_diam = np.zeros(num)
        npart = len(idxs) 

        mu = 0.
        uncert  = self.get_config('radionuclide:particle_diameter_uncertainty')
        if diam>0:
            uncert_ln = uncert/diam * 3.
        else:
            uncert_ln = 0.

        rng = np.random.default_rng()

        if distr=='normal':
            truncs = rng.normal(mu, uncert, npart)
            init_diam[idxs] = diam + truncs

        elif distr=='lognormal':
            truncs = rng.lognormal(mu, uncert_ln, size=npart ) # * diameter
            init_diam[idxs] = diam * truncs
        else:
            logger.error('Not available distribution: %s',distr)

        pmin = self.get_config('radionuclide:particle_diameter_minimum')
        pmax = self.get_config('radionuclide:particle_diameter_maximum')
        init_diam[idxs][init_diam[idxs]<pmin] = pmin 
        init_diam[idxs][init_diam[idxs]>pmax] = pmax


        return init_diam




    def check_speciation(self):
        isotop = self.get_config('radionuclide:isotope')
        specie_setup = self.get_config('radionuclide:specie_setup')

        legal_species = { '137Cs'   : ['LMM + Rev', 
                                       'LMM + Rev + Slow rev', 
                                       'LMM + Rev + Irrev', 
                                       'LMM + Rev + Slow rev + Irrev'],
                      '129I'    : ['LMM + Rev', 
                                   'LMM + Rev + Slow rev + Irrev'],
                      '241Am'   : ['LMM + Rev', 
                                   'LMM + Rev + Slow rev', 
                                    'LMM + Rev + Slow rev + Irrev'],
                      'Al'      : ['LMM + Colloid + Rev'],
                     }


        if specie_setup in legal_species[isotop]:
            logger.info(f'All good, isotop: {isotop}, species: {specie_setup}')
        else:
            logger.info(f'Not a legal combination of isotop: {isotop} and speciation: {specie_setup}')
            logger.error('Illegal speciation for %s: %s',isotop, specie_setup)
            exit()






    def seed_elements(self, *args, **kwargs):


        self.check_speciation()
        self.init_species()
        self.init_transfer_rates()


        if 'number' in kwargs:
            num_elements = kwargs['number']
        else:
            num_elements = self.get_config('seed:number')




        # Initialize species 
        if 'specie' in kwargs:
            print('num_elements', num_elements)
            try:
                print('len specie:',len(kwargs['specie']))
            except:
                print('specie:',kwargs['specie'])

            init_specie = np.ones(num_elements,dtype=int)
            init_specie[:] = kwargs['specie']
            kwargs['specie'] = init_specie[:]

        else:

            # Set initial speciation
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
            dissolved = np.random.rand(num_elements)<lmm_frac
            ndiss = np.sum(dissolved)
            if 'LMMcation' in self.name_species:
                init_specie[dissolved] = self.num_lmmcation
            elif 'LMM' in self.name_species:
                init_specie[dissolved] = self.num_lmm
            else:
                logger.error('No dissolved species present %s: ',self.name_species)
            
            if self.get_config('seed:slowly_fraction')>0 and 'Particle slowly reversible' in self.name_species:
                ndiss = np.sum(dissolved)
                prtidx = np.where(~dissolved)[0]
                slowprt  = np.random.rand(num_elements-ndiss)<self.get_config('seed:slowly_fraction')
                init_specie[prtidx[slowprt]] = self.num_psrev
                init_specie[prtidx[~slowprt]] = self.num_prev
            else:
                init_specie[~dissolved] = self.num_prev

            kwargs['specie'] = init_specie


        logger.info('Initial speciation:')
        for i,sp in enumerate(self.name_species):
            logger.info( '{:>9} {:>3} {:24} '.format(  np.sum(init_specie==i), i, sp ) )


        if all(init_specie == self.num_srev):
            print( 'ALL ELEMENTS ARE SEDIMENTS')
            kwargs['z'] = 'seafloor'
            # TODO SET MOVING = 0 

        elif any(init_specie == self.num_srev):
            print( 'SOME ELEMENTS ARE SEDIMENTS')
            exit()



        # Set initial particle size according to speciation
        if 'diameter' in kwargs:
            diameter = kwargs['diameter']
        else:
            diameter = self.get_config('radionuclide:particle_diameter')

        logger.info('PARTICLE DIAMETER %s',diameter)

        particles = []
        for i,sp in enumerate(self.name_species):
            if 'particle' in sp.casefold():
                prt = np.where(init_specie==i)[0]
                particles.extend(prt)

        init_diam = self.set_init_diameter(num_elements, particles, diameter)


        try:
            logger.info('Min: {:.2e}, Max: {:.2e}, Numer of particles seeded: {}, at {} distribution'.format(
                np.min(init_diam[init_diam>0.]), np.max(init_diam[init_diam>0.]), sum(init_diam>0.),
                self.get_config('radionuclide:particlesize_distribution')))
        except Exception:
            logger.info('All diameters are 0 (min:{}, max: {})'.format( np.min(init_diam), np.max(init_diam)) )

        kwargs['diameter'] = init_diam

        
        
        # Set z to seabed for sediments
        # sediments = []
        # for i,sp in enumerate(self.name_species):
        #     if 'sediment' in sp.casefold():
        #         sed = np.where(init_specie==i)[0]
        #         sediments.extend(sed)
        
        # kwargs['z'] = 'seafloor'
        # logger.debug('Seafloor is selected, neglecting z')
#         init_z = np.ones(num_elements,dtype=float)
#         init_z[:] = kwargs['z']
#         init_z[sediments] = -10000. #\
#                 #-1.*self.environment.sea_floor_depth_below_sea_level[sediments]
# #        self.elements.moving[sediments] = 0
#         kwargs['z'] = init_z


        super(RadionuclideDrift, self).seed_elements(*args, **kwargs)



    def init_kd(self):
        ''' Initialization of Kd value, dependent on simulated radionuclide
        '''

        # Values from IAEA (2004)
        kd_values = { '137Cs'   : 4.0e0,
                      '129I'    : 7.0e-2, 
                      '241Am'   : 2.0e3,
                      'Al'      : None,
                     }


        print('ISOTOPE',self.isotope)

        if self.isotope == 'manual':
            self.kd = self.get_config('radionuclide:transformations:Kd')
            return

        if not self.isotope in kd_values.keys():
            logger.error('No Kd value implemented for %s ', self.isotope)

        else:
            self.kd = kd_values[self.isotope]






    def init_transfer_rates(self):
        ''' Initialization of background values in the transfer rates 2D array.
        '''


#        logger.info( 'transfer setup: %s' % transfer_setup)

        self.isotope = self.get_config('radionuclide:isotope')

        self.init_kd()

        logger.info('ISOTOPE %s',self.isotope)
        logger.info('SPECIE SETUP %s',self.specie_setup)
        logger.info('Kd: %s',self.kd)

        self.transfer_rates = np.zeros([self.nspecies,self.nspecies])
        self.ntransformations = np.zeros([self.nspecies,self.nspecies])


        if 'lmm + rev' in self.specie_setup.casefold(): 

            self.num_lmm    = self.specie_name2num('LMM')
            self.num_prev   = self.specie_name2num('Particle reversible')
            self.num_srev   = self.specie_name2num('Sediment reversible')


            # Simpler version of Values from Simonsen et al (2019a)
            # Only consider the reversible fraction
            Kd        = self.kd     # depends on isotope
            Dc         = self.get_config('radionuclide:transformations:Dc')
            susp_mat    = 1.e-3   # concentration of available suspended particulate matter (kg/m3)
            sedmixdepth = self.get_config('radionuclide:sediment:sedmixdepth')     # sediment mixing depth (m)
            default_density =  self.get_config('radionuclide:sediment:sediment_density') # default particle density (kg/m3)
            f           =  self.get_config('radionuclide:sediment:effective_fraction')      # fraction of effective sorbents
            phi         =  self.get_config('radionuclide:sediment:corr_factor')      # sediment correction factor
            poro        =  self.get_config('radionuclide:sediment:porosity')      # sediment porosity
            layer_thick =  self.get_config('radionuclide:sediment:layer_thick')      # thickness of seabed interaction layer (m)

            self.transfer_rates[self.num_lmm,self.num_prev] = Dc * Kd * susp_mat
            self.transfer_rates[self.num_prev,self.num_lmm] = Dc
            self.transfer_rates[self.num_lmm,self.num_srev] = \
                Dc * Kd * sedmixdepth * default_density * (1.-poro) * f * phi / layer_thick
            self.transfer_rates[self.num_srev,self.num_lmm] = Dc * phi

        else:
            logger.error('No transfer setup available')


        if '+ slow rev' in self.specie_setup.casefold(): 

            self.num_psrev  = self.specie_name2num('Particle slowly reversible')
            self.num_ssrev  = self.specie_name2num('Sediment slowly reversible')

            slow_coeff = self.get_config('radionuclide:transformations:slow_coeff')

            self.transfer_rates[self.num_srev,self.num_ssrev] = slow_coeff
            self.transfer_rates[self.num_prev,self.num_psrev] = slow_coeff
            self.transfer_rates[self.num_ssrev,self.num_srev] = slow_coeff*.1
            self.transfer_rates[self.num_psrev,self.num_prev] = slow_coeff*.1

        if '+ irrev' in self.specie_setup.casefold():

            self.num_pirrev  = self.specie_name2num('Particle irreversible')
            self.num_sirrev  = self.specie_name2num('Sediment irreversible')

            slow_coeff = self.get_config('radionuclide:transformations:slow_coeff')
            
            self.transfer_rates[self.num_ssrev,self.num_sirrev] = slow_coeff
            self.transfer_rates[self.num_psrev,self.num_pirrev] = slow_coeff




        if self.specie_setup.casefold() == 'lmm + colloid + rev' and \
               self.isotope.casefold() == 'al':
            # Use values from Simonsen et al (2019b)
            self.num_lmmanion    = self.specie_name2num('LMManion')
            self.num_lmmcation   = self.specie_name2num('LMMcation')
            self.num_humcol      = self.specie_name2num('Humic colloid')
            self.num_polymer     = self.specie_name2num('Polymer')
            self.num_prev        = self.specie_name2num('Particle reversible')
            self.num_srev        = self.specie_name2num('Sediment reversible')

            Dc         = self.get_config('radionuclide:transformations:Dc')

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




        # Set diagonal to 0. (not possible to transform to present specie)
        if len(self.transfer_rates.shape) == 3:
            for ii in range(self.transfer_rates.shape[0]):
                np.fill_diagonal(self.transfer_rates[ii,:,:],0.)
        else:
            np.fill_diagonal(self.transfer_rates,0.)










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

        self.elements.terminal_velocity = W * self.elements.moving






    def update_transfer_rates(self):
        '''Pick out the correct row from transfer_rates for each element. Modify the
        transfer rates according to local environmental conditions '''

        if not self.specie_setup.casefold() == 'lmm + colloid + rev':
            self.elements.transfer_rates1D = self.transfer_rates[self.elements.specie,:]

            if 'Sediment reversible' in self.name_species:
                # Only LMM radionuclides close to seabed are allowed to interact with sediments
                # minimum height/maximum depth for each particle
                Zmin = -1.*self.environment.sea_floor_depth_below_sea_level
                interaction_thick = self.get_config('radionuclide:sediment:layer_thick')      # thickness of seabed interaction layer (m)
                dist_to_seabed = self.elements.z - Zmin
                self.elements.transfer_rates1D[(self.elements.specie == self.num_lmm) &
                                 (dist_to_seabed > interaction_thick), self.num_srev] = 0.

            if 'Particle reversible' in self.name_species:
                # Modify particle adsorption according to local particle concentration
                # (LMM -> reversible particles)
                kktmp = self.elements.specie == self.num_lmm
                self.elements.transfer_rates1D[kktmp, self.num_prev] = \
                            self.elements.transfer_rates1D[kktmp, self.num_prev] * \
                            self.environment.conc3[kktmp] / 1.e-3
        #                    self.environment.particle_conc[kktmp] / 1.e-3

        elif self.specie_setup.casefold() == 'lmm + colloid + rev' and self.isotope=='Al':
            sal = self.environment.sea_water_salinity
            sali = np.searchsorted(self.salinity_intervals, sal) - 1
            self.elements.transfer_rates1D = self.transfer_rates[sali,self.elements.specie,:]



    def update_speciation(self):
        '''Check if transformation processes shall occur
        Do transformation (change value of self.elements.specie)
        Update element properties for the transformed elements
        '''

        specie_in  = self.elements.specie.copy()    # for storage of the out speciation
        specie_out = self.elements.specie.copy()    # for storage of the out speciation
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


        # Set the new speciation
        self.elements.specie=specie_out

        logger.debug('old species: %s' % specie_in[phaseshift])
        logger.debug('new species: %s' % specie_out[phaseshift])


        for iin in range(self.nspecies):
            for iout in range(self.nspecies):
                self.ntransformations[iin,iout]+=sum((specie_in[phaseshift]==iin) & (specie_out[phaseshift]==iout))

        logger.debug('Number of transformations total:\n %s' % self.ntransformations )


        # Update radionuclide properties after transformations
        self.update_radionuclide_diameter(specie_in, specie_out)
        self.sorption_to_sediments(specie_in, specie_out)
        self.desorption_from_sediments(specie_in, specie_out)



    def sorption_to_sediments(self,sp_in=None,sp_out=None):
        '''Update radionuclide properties  when sorption to sediments occurs'''


        # Set z to local sea depth
        if 'LMM' in self.name_species:
            self.elements.z[(sp_out==self.num_srev) & (sp_in==self.num_lmm)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_srev) & (sp_in==self.num_lmm)]
            self.elements.moving[(sp_out==self.num_srev) & (sp_in==self.num_lmm)] = 0
        if 'LMMcation' in self.name_species:
            self.elements.z[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)]
            self.elements.moving[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)] = 0
        # avoid setting positive z values
        if np.nansum(self.elements.z>0):
            logger.debug('Number of elements lowered down to sea surface: %s' % np.nansum(self.elements.z>0))
        self.elements.z[self.elements.z > 0] = 0



    def desorption_from_sediments(self,sp_in=None,sp_out=None):
        '''Update radionuclide properties when desorption from sediments occurs'''

        desorption_depth = self.get_config('radionuclide:sediment:desorption_depth')
        std = self.get_config('radionuclide:sediment:desorption_depth_uncert')


        if 'LMM' in self.name_species:
            self.elements.z[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] + desorption_depth
            self.elements.moving[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] = 1
            if std > 0:
                logger.debug('Adding uncertainty for desorption from sediments: %s m' % std)
                self.elements.z[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] += np.random.normal(
                        0, std, sum((sp_out==self.num_lmm) & (sp_in==self.num_srev)))
        if 'LMMcation' in self.name_species:
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





    def update_radionuclide_diameter(self,sp_in=None,sp_out=None):
        '''Update the diameter of the radionuclides when specie is changed'''


        dia_part=self.get_config('radionuclide:particle_diameter')
        dia_diss=self.get_config('radionuclide:dissolved_diameter')



        # Transform to particle species
        to_particles = np.where( (sp_out==self.num_prev) & (sp_in!=self.num_prev) ) [0]
        if self.species_slowly_fraction:
            to_particles = np.append(to_particles, np.where((sp_out==self.num_psrev) & (sp_in==self.num_ssrev))[0] )
        if self.species_irreversible_fraction:
            to_particles = np.append(to_particles, np.where((sp_out==self.num_pirrev) & (sp_in==self.num_sirrev))[0] )

        new_diam = self.set_init_diameter(self.num_elements_total(), to_particles, dia_part)
        self.elements.diameter[to_particles] = new_diam[to_particles]




        # Transform to LMM and colloid species
        if 'LMM' in self.name_species:
            to_diss = np.where( (sp_out==self.num_lmm) & (sp_in!=self.num_lmm) )[0]
        elif 'LMMcation' and 'LMManion' in self.name_species: 
            to_diss = np.where( (sp_out==self.num_lmmcation) & (sp_in!=self.num_lmmcation) )[0]
            to_diss = np.append( to_diss, np.where( (sp_out==self.num_lmmanion) & (sp_in!=self.num_lmmanion) )[0] )
        if 'Colloid' in self.name_species:
            to_diss = np.append(to_diss, np.where( (sp_out==self.num_col) & (sp_in!=self.num_col) )[0] )
        if 'Humic_colloid' in self.name_species:
            to_diss = np.append( to_diss, np.where( (sp_out==self.num_humcol) & (sp_in!=self.num_humcol) )[0] )
        if 'Polymer' in self.name_species:
            to_diss = np.append( to_diss, np.where( (sp_out==self.num_polymer) & (sp_in!=self.num_polymer) )[0] )

        new_diam = self.set_init_diameter(self.num_elements_total(), to_diss, dia_diss)
        self.elements.diameter[to_diss] = new_diam[to_diss]









    def bottom_interaction(self,Zmin=None):
        ''' Change speciation of radionuclides that reach bottom due to settling.
        particle specie -> sediment specie '''

        if not ( ('Particle reversible' in self.name_species) & 
                ('Sediment reversible' in self.name_species) or 
                self.species_slowly_fraction or self.species_irreversible_fraction):
            print('RETURN: ', ('Particle reversible' in self.name_species),  
                ('Sediment reversible' in self.name_species),
                self.species_slowly_fraction,  
                self.species_irreversible_fraction)
            return

        bottom = np.array(np.where(self.elements.z <= Zmin)[0])
        kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_prev)[0])
        self.elements.specie[bottom[kktmp]] = self.num_srev
        self.ntransformations[self.num_prev,self.num_srev]+=len(kktmp)
        self.elements.moving[bottom[kktmp]] = 0
#        self.elements.moving[bottom] = 0

        if self.species_slowly_fraction:
            kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_psrev)[0])
            self.elements.specie[bottom[kktmp]] = self.num_ssrev
            self.ntransformations[self.num_psrev,self.num_ssrev]+=len(kktmp)
            self.elements.moving[bottom[kktmp]] = 0

        if self.species_irreversible_fraction:
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
        # if not  ((self.get_config('radionuclide:species:Particle_reversible')) &
        #           (self.get_config('radionuclide:species:Sediment_reversible'))):
        #     return
        if not ('Particle reversible' in self.name_species and 'Sediment reversible' in self.name_species):
            logger.info('No sediments and particles present...')
            return

        specie_in = self.elements.specie.copy()

        critvel = self.get_config('radionuclide:sediment:resuspension_critvel')
        resusp_depth = self.get_config('radionuclide:sediment:resuspension_depth')
        std = self.get_config('radionuclide:sediment:resuspension_depth_uncert')

        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level
        x_vel = self.environment.x_sea_water_velocity
        y_vel = self.environment.y_sea_water_velocity
        speed = np.sqrt(x_vel*x_vel + y_vel*y_vel)
        bottom = (self.elements.z <= Zmin)

        resusp = ( (bottom) & (speed >= critvel) )
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
#        if self.get_config('radionuclide:slowly_fraction'):
        if self.species_slowly_fraction:
            self.ntransformations[self.num_ssrev,self.num_psrev]+=sum((resusp) & (self.elements.specie==self.num_ssrev))
            self.elements.specie[(resusp) & (self.elements.specie==self.num_ssrev)] = self.num_psrev
        if self.species_irreversible_fraction:
#        if self.get_config('radionuclide:irreversible_fraction'):
            self.ntransformations[self.num_sirrev,self.num_pirrev]+=sum((resusp) & (self.elements.specie==self.num_sirrev))
            self.elements.specie[(resusp) & (self.elements.specie==self.num_sirrev)] = self.num_pirrev

        specie_out = self.elements.specie.copy()
        self.update_radionuclide_diameter(specie_in, specie_out)






    def update(self):
        """Update positions and properties of radionuclide particles."""

        # Workaround due to conversion of datatype
        self.elements.specie = self.elements.specie.astype(np.int32)

        # Radionuclide speciation
        self.update_transfer_rates()
        self.update_speciation()


        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
        else:
            self.update_terminal_velocity()
            self.vertical_buoyancy()


        # Resuspension
        self.resuspension()
        logger.info('Speciation: {} {}'.format([sum(self.elements.specie==ii) for ii in range(self.nspecies)],self.name_species))


        #self.elements.moving[self.elements.z <= -1*self.environment.sea_floor_depth_below_sea_level] = 0.

        
        # Horizontal advection
        self.advect_ocean_current()


        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()







# ################
# POSTPROCESSING


    def write_netcdf_radionuclide_density_map(self, filename, pixelsize_m='auto', zlevels=None,
                                              deltat=None,
                                              density_proj=None,
                                              llcrnrlon=None, llcrnrlat=None,
                                              urcrnrlon=None, urcrnrlat=None,
                                              activity_unit=None,
                                              time_avg_conc=False,
                                              horizontal_smoothing=False,
                                              smoothing_cells=0,
                                              reader_sea_depth=None,
                                              ):
        '''Write netCDF file with map of radionuclide species densities and concentrations'''

        from netCDF4 import Dataset, date2num #, stringtochar

        logger.info('Postprocessing: Write density and concentration to netcdf file')

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




        if reader_sea_depth is not None:
            from opendrift.readers import reader_netCDF_CF_generic,reader_ROMS_native
            reader = reader_sea_depth
        else:
            for readerName in self.env.readers:
                reader = self.env.readers[readerName]
                if 'sea_floor_depth_below_sea_level' in reader.variables:
                    break

       
        try:
            self.conc_lat  = reader.get_variables('latitude',
                                        x=[reader.xmin,reader.xmax],
                                        y=[reader.ymin,reader.ymax])['latitude'][:]
        except:
            self.conc_lat  = reader.get_variables('grid_latitude_at_cell_center',
                                        x=[reader.xmin,reader.xmax],
                                        y=[reader.ymin,reader.ymax])['grid_latitude_at_cell_center'][:]
        try:
            self.conc_lon  = reader.get_variables('longitude',
                                        x=[reader.xmin,reader.xmax],
                                        y=[reader.ymin,reader.ymax])['longitude'][:]
        except:
            self.conc_lon  = reader.get_variables('grid_longitude_at_cell_center',
                                        x=[reader.xmin,reader.xmax],
                                        y=[reader.ymin,reader.ymax])['grid_longitude_at_cell_center'][:]

        self.conc_topo = reader.get_variables('sea_floor_depth_below_sea_level',
                                        x=[reader.xmin,reader.xmax],
                                        y=[reader.ymin,reader.ymax])['sea_floor_depth_below_sea_level'][:]




        if density_proj is None: # add default projection with equal-area property
            density_proj = pyproj.Proj('+proj=moll +ellps=WGS84 +lon_0=0.0')



        if activity_unit==None:
            activity_unit='Bq'  # default unit for radionuclides

        activity_per_element = self.get_config('seed:total_release') / self.num_elements_total()


        z = self.get_property('z')[0]
        if not zlevels==None:
            zlevels = np.sort(zlevels)
            z_array = np.append(np.append(-10000, zlevels) , max(0,np.nanmax(z)))
        else:
            z_array = [min(-10000,np.nanmin(z)), max(0,np.nanmax(z))]
        logger.info('z_array: {}'.format(  [str(item) for item in z_array] ) )




        #
        # H is array containing number of elements within each box defined by lon_array, lat_array and z_array

        H, lon_array, lat_array = \
            self.get_radionuclide_density_array(pixelsize_m, z_array,
                                                density_proj=density_proj,
                                                llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                                                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat
                                                )

        lon_array = (lon_array[:-1,:-1] + lon_array[1:,1:])/2
        lat_array = (lat_array[:-1,:-1] + lat_array[1:,1:])/2



        if horizontal_smoothing:
            # Compute horizontally smoother field
            logger.info('H.shape: ' + str(H.shape))
            Hsm = np.zeros_like(H)
            for zi in range(len(z_array)-1):
                for sp in range(self.nspecies):
                    for ti in range(H.shape[0]):
                        Hsm[ti,sp,zi,:,:] = self.horizontal_smooth(H[ti,sp,zi,:,:],n=smoothing_cells)



        # Convert from density to concentration
        logger.info('Activity: '+str(activity_per_element)+' '+ activity_unit+ ' per unit')

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





        conc = np.zeros_like(H)
        if horizontal_smoothing:
            conc_sm = np.zeros_like(Hsm)
        for ti in range(H.shape[0]):
            for sp in range(self.nspecies):
                conc[ti,sp,:,:,:] = H[ti,sp,:,:,:] / pixel_volume * activity_per_element
                if horizontal_smoothing:
                    conc_sm[ti,sp,:,:,:] = Hsm[ti,sp,:,:,:] / pixel_volume * activity_per_element




        times = np.array( self.get_time_array()[0] )
        if time_avg_conc:
            conctmp = conc[:-1,:,:,:,:]
            cshape = conctmp.shape
            mdt =    np.mean(times[1:] - times[:-1])    # output frequency in opendrift output file
            if deltat==None:
                ndt = 1
            else:
                ndt = int( deltat / (mdt.total_seconds()/3600.) )
            times2 = times[::ndt]
            times2 = times2[1:]
            odt = int(cshape[0]/ndt)
            logger.info ('ndt '+ str(ndt))   # number of time steps over which to average in conc file
            logger.info ('odt '+ str(odt))   # number of average slices


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

        #
        nc.createVariable('concfactor','f8')
        nc.variables['concfactor'][:] = activity_per_element
        nc.variables['concfactor'].long_name = 'Activity per unit element'
        nc.variables['concfactor'].unit = activity_unit

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
        outstr = ['{}:{}'.format(isp,sp) for isp,sp in enumerate(self.name_species)]
        nc.variables['specie'].names = outstr



        # Density
        nc.createVariable('density', 'i4',
                          ('time','specie','depth','y', 'x'),fill_value=-999)
        H = np.swapaxes(H, 3, 4).astype('i4')
        H = np.ma.masked_where(H==0, H)
        nc.variables['density'][:] = H
        nc.variables['density'].long_name = 'Number of elements in grid cell'
        nc.variables['density'].grid_mapping = 'projection'
        nc.variables['density'].units = '1'


        if horizontal_smoothing:
            nc.createVariable('density_smooth', 'f8',
                              ('time','specie','depth','y', 'x'),fill_value=1.e36)
            Hsm = np.swapaxes(Hsm, 3, 4).astype('f8')
            #Hsm = np.ma.masked_where(Hsm==0, Hsm)
            nc.variables['density_smooth'][:] = Hsm
            nc.variables['density_smooth'].long_name = 'Horizontally smoothed number of elements in grid cell'
            nc.variables['density_smooth'].comment = 'Smoothed over '+str(smoothing_cells)+' grid points in all horizontal directions'



        # Radionuclide concentration, horizontally smoothed
        nc.createVariable('concentration', 'f8',
                          ('time','specie','depth','y', 'x'),fill_value=1.e36)
        conc = np.swapaxes(conc, 3, 4) #.astype('i4')
        #conc = np.ma.masked_where(conc==0, conc)
        nc.variables['concentration'][:] = conc
        nc.variables['concentration'].long_name = 'Radionuclide concentration'
        nc.variables['concentration'].grid_mapping = 'projection_lonlat'
        nc.variables['concentration'].units = activity_unit+'/m3'



        if horizontal_smoothing:
            # Radionuclide concentration, horizontally smoothed
            nc.createVariable('concentration_smooth', 'f8',
                              ('time','specie','depth','y', 'x'),fill_value=1.e36)
            conc_sm = np.swapaxes(conc_sm, 3, 4) #.astype('i4')
          #  conc_sm = np.ma.masked_where(conc_sm==0, conc_sm)
            nc.variables['concentration_smooth'][:] = conc_sm
            nc.variables['concentration_smooth'].long_name = 'Horizontally smoothed radionuclide concentration'
            nc.variables['concentration_smooth'].grid_mapping = 'projection_lonlat'
            nc.variables['concentration_smooth'].units = activity_unit+'/m3'
            nc.variables['concentration_smooth'].comment = 'Smoothed over '+str(smoothing_cells)+' grid points in all horizontal directions'



        if time_avg_conc:
            nc.createVariable('concentration_avg', 'f8',
                              ('avg_time','specie','depth','y', 'x'),fill_value=0)
            conc2 = np.swapaxes(mean_conc, 3, 4) #.astype('i4')
            conc2 = np.ma.masked_where(conc2==0, conc2)
            nc.variables['concentration_avg'][:] = conc2
            nc.variables['concentration_avg'].long_name = 'Time averaged radionuclide concentration'
            nc.variables['concentration_avg'].grid_mapping = 'projection_lonlat'
            nc.variables['concentration_avg'].units = activity_unit+'/m3'


        # Volume of boxes
        nc.createVariable('volume', 'f8',
                          ('depth','y', 'x'),fill_value=0)
        pixel_volume = np.swapaxes(pixel_volume, 1, 2) #.astype('i4')
        pixel_volume = np.ma.masked_where(pixel_volume==0, pixel_volume)
        nc.variables['volume'][:] = pixel_volume
        nc.variables['volume'].long_name = 'Volume of grid cell'
        nc.variables['volume'].grid_mapping = 'projection_lonlat'
        nc.variables['volume'].units = 'm3'


        # Topography
        nc.createVariable('topo', 'f8', ('y', 'x'),fill_value=0)
        pixel_mean_depth = np.ma.masked_where(pixel_mean_depth==0, pixel_mean_depth)
        nc.variables['topo'][:] = pixel_mean_depth.T
        nc.variables['topo'].long_name = 'Depth of grid point'
        nc.variables['topo'].grid_mapping = 'projection_lonlat'
        nc.variables['topo'].units = 'm'


        nc.close()
        logger.info('Wrote to '+filename)





    def get_radionuclide_density_array(self, pixelsize_m, z_array,
                                       density_proj=None, llcrnrlon=None,llcrnrlat=None,
                                       urcrnrlon=None,urcrnrlat=None,
                                       weight=None):
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
                else:
                    weights = None
                for zi in range(len(z_array)-1):
                    kktmp = ( (specie[i,:]==sp) & (z[i,:]>z_array[zi]) & (z[i,:]<=z_array[zi+1]) )
                    H[i,sp,zi,:,:], dummy, dummy = \
                        np.histogram2d(x[i,kktmp], y[i,kktmp],
                                   weights=weights, bins=bins)

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
        h = interpolate.griddata((lon_grd[h_grd.mask==False].flatten(),lat_grd[h_grd.mask==False].flatten()), h_grd[h_grd.mask==False].flatten(), (lons, lats), method='linear')

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







    def gui_postproc(self):
        from os.path import expanduser

        logger.info('Postprocessing radionuclides')

        logger.info('Final speciation:')
        for isp,sp in enumerate(self.name_species):
            logger.info ('{:32}: {:>6}'.format(sp,sum(self.elements.specie==isp)))

        homefolder = expanduser("~")
        filename = homefolder+'/conc_radio_gui.nc'


        self.guipp_saveconcfile(filename)

    """
        zlayer   = [-1,-2]
        self.guipp_plotandsaveconc(filename=filename,
                                         outfilename=homefolder+'/radio_plots/RadioConc',
                                         zlayers=zlayer,
                                         specie= ['Total', 'LMM', 'Particle reversible'],
                                         )
    """





    def guipp_saveconcfile(self,filename='none', pixelsize_m=200., reader_sea_depth=None):



        if reader_sea_depth is not None:
            from opendrift.readers import reader_netCDF_CF_generic
            reader_sea_depth = reader_netCDF_CF_generic.Reader(reader_sea_depth)
        else:
            for readerName in self.env.readers:
                reader = self.env.readers[readerName]
                if 'sea_floor_depth_below_sea_level' in reader.variables:
                    break


        self.conc_lon  = reader.get_variables('longitude',
                                       x=[reader.xmin,reader.xmax],
                                       y=[reader.ymin,reader.ymax])['longitude'][:]

        self.conc_lat  = reader.get_variables('latitude',
                                       x=[reader.xmin,reader.xmax],
                                       y=[reader.ymin,reader.ymax])['latitude'][:]

        self.conc_topo  = reader.get_variables('sea_floor_depth_below_sea_level',
                                       x=[reader.xmin,reader.xmax],
                                       y=[reader.ymin,reader.ymax])['sea_floor_depth_below_sea_level'][:]


        self.write_netcdf_radionuclide_density_map(filename, pixelsize_m=pixelsize_m,
                                        zlevels=self.depthintervals,
                                        activity_unit='Bq',
                                        horizontal_smoothing=True,
                                        smoothing_cells=1,
                                        )


    def guipp_showanimationprofile(self):
        self.animation_profile(color='specie',
            vmin=0,vmax=self.nspecies-1,
            legend=[self.specie_num2name(i) for i in range(self.nspecies)],
            legend_loc =3,
#            markersize=10
            )






    def guipp_plotandsaveconc(self, filename=None, outfilename=None, zlayers=None, time=None, specie=None):
        import netCDF4
        from datetime import timedelta, datetime
        import cartopy.crs as ccrs
        import matplotlib.pyplot as plt
        from os.path import expanduser


        if zlayers==None:
            zlayers=[-1]


        logger.info('Plotting concentration for {} at layer {} at time {}'.format(specie, zlayers, time))

        homefolder = expanduser("~")
        outdir = homefolder+'/radio_plots'
        vcmin  = 0
        vcmaxt = self.get_config('seed:total_release') * 1.e-8
        cblabel= 'Bq/m$^3$'

        proj_pp=ccrs.PlateCarree()
        cmap = 'CMRmap_r'


        if specie==None:
            specie_arr = ['Total', 'LMM', 'Particle reversible', 'Sediment reversible']
        else:
            specie_arr = specie

        # for readerName in self.readers:
        #     reader = self.readers[readerName]
        #     if 'sea_floor_depth_below_sea_level' in reader.variables:
        #         break

        # gtopo  = reader.get_variables('sea_floor_depth_below_sea_level',
        #                                x=[reader.xmin,reader.xmax],
        #                                y=[reader.ymin,reader.ymax])['sea_floor_depth_below_sea_level'][:]



        logger.info (' READ DATA FILE: {}'.format(filename))
        nc=netCDF4.Dataset(filename,'r')
        time       = nc.variables['time'][:]
        t_units    = nc.variables['time'].units
        #conc       = nc.variables['concentration_smooth'][:]; print('concentration_smooth')
        specie     = nc.variables['specie'][:]
        sp_names   = nc.variables['specie'].names
        lat        = nc.variables['lat'][:]
        lon        = nc.variables['lon'][:]
        depth      = nc.variables['depth'][:]
        topo       = nc.variables['topo'][:]
        #nc.close()



        Nspecies = len(specie[:])
        logger.info ('Nspecies: {}'.format(Nspecies))
        if Nspecies>1:
            specie_names = [i.split(':')[1] for i in sp_names]
        else:
            specie_names = [sp_names.split(':')[1]]
        logger.info(specie_names)




#        Ndepths = len(depth)
#        logger.info ('Ndepths: {}'.format(Ndepths) )
#        [print(item) for item in depth]
#         if Ndepths==1:
#             figsize=[5,6]
#         elif Ndepths==2:
#             figsize=[10,8]
#         elif Ndepths==3:
#             figsize=[13,8]


        ntime=len(time)
        time = netCDF4.num2date(time,units=t_units,only_use_cftime_datetimes=False,only_use_python_datetimes=True)


        for sp in specie_arr:

            fig, ax, crs, x, y, index_of_first, index_of_last = \
                self.set_up_map()

            cb=False

            if not sp=='Total':
                try:
                    spi = specie_names.index(sp)
                    logger.debug('{} {}'.format(sp, spi))
                except Exception:
                    logger.debug(sp)
                    break


            for date_idx in np.arange(0,len(time)):

                idx1=date_idx
                t1 = time[idx1]

                for zlayer in zlayers:
                    if sp=='Total':
                        tmp = np.sum(nc.variables['concentration_smooth'][idx1,:,zlayer,:,:],axis=0)
                    else:
                        tmp = nc.variables['concentration_smooth'][idx1,spi,zlayer,:,:]
                    pcm=ax.pcolormesh(lon,lat,tmp,
                                      vmin=vcmin,vmax=vcmaxt,
                                      #norm=colors.LogNorm(vmin=vcmin,vmax=vcmax),
                                      cmap=cmap, zorder=-1,
                                      transform=proj_pp
                                      )
                    c1 = ax.contour(lon,lat,tmp,
                                   [1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.e0,1.e1,1.e2,1.e3,1.e4],
                                        colors='darkgrey',linewidths=.6,zorder=4,
                                         transform=proj_pp)
                    if not cb:
                        plt.colorbar(pcm,extend='both',label=cblabel)
                        cb=True
                    title=sp+' concentration, layer '+str(depth[zlayer])+' '+t1.strftime('%Y%m%d %H:%M')
                    fig.suptitle(title)
                    title = title.replace(' ','').replace('\n','_').replace(':','')

                    if outfilename:
                        ofn = outfilename+'_'+sp.replace(' ','')+'L'+str(depth[zlayer])+'_'+t1.strftime('%Y%m%dT%H%M')+'.png'
                    else:
                        ofn = outdir+'/'+title+'_'+'.png'
                    plt.savefig(ofn, dpi=300,bbox_inches='tight')

                    # remove contour lines and shading
                    for level in c1.collections:
                        level.remove()
                    pcm.remove()
                    logger.info(ofn)
            plt.close()
        logger.info('DONE')



