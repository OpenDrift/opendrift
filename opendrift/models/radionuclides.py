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

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
#from opendrift.elements import LagrangianArray

from opendrift.readers.basereader import pyproj

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
                    'default': 0})#,
#         ('transfer_rates1D', {'dtype':np.array(3, dtype=np.float32),
#                     'units': '1/s',
#                     'default': 0.})
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

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
#                          'sea_surface_wave_significant_height',
#                          'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature',
                          'sea_water_salinity',
                          'surface_downward_x_stress',
                          'surface_downward_y_stress',
                          'turbulent_kinetic_energy',
                          'turbulent_generic_length_scale',
                          'upward_sea_water_velocity',
                          'conc3'
                          ]

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented)
    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-20, 0]

    fallback_values = {#'x_sea_water_velocity': 0,
                       #'y_sea_water_velocity': 0,
                       #'sea_surface_wave_significant_height': 0,
                       #'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       #'sea_floor_depth_below_sea_level': 100,
#                       'ocean_vertical_diffusivity': 0.02,  # m2s-1
                       'ocean_vertical_diffusivity': 0.0001,  # m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity': 34.,
                       'surface_downward_x_stress': 0,
                       'surface_downward_y_stress': 0,
                       'turbulent_kinetic_energy': 0,
                       'turbulent_generic_length_scale': 0,
                       'upward_sea_water_velocity': 0,
                       'conc3': 1.e-3
                       }

    # Default colors for plotting
#    status_colors = {'initial': 'green', 'active': 'blue',
#                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta', 'sedimented':'sandybrown'}

#    specie_colors = {'LMM':'blue',
#                      'Colloid':'Sandybrown',
#                      'Particle reversible': 'grey',
#                      'Particle irreversible':'darkgrey',
#                      'Sediment reversible':'red',
#                      'Sediment irreversible':'darkred'
#                      }


#     species_all = {
#                     0:{'name':'LMM','color':'blue'},
#                     1:{'name':'Colloid','color':'Sandybrown'},
#                     2:{'name':'Particle reversible','color':'grey'},
#                     3:{'name':'Particle slowly reversible','color':'darkgrey'},
#                     4:{'name':'Particle irreversible','color':'black'},
#                     5:{'name':'Sediment reversible','color':'red'},
#                     6:{'name':'Sediment slowly reversible','color':'lightred'},
#                     7:{'name':'Sediment irreversible','color':'darkred'}
#                     }

    configspec_radionuclidedrift = '''
        [drift]
            vertical_mixing = boolean(default=True)
        [radionuclide]
            transfer_setup = option('Sandnesfj_Al','Bokna_137Cs', 'dummy', default='dummy')
            slowly_fraction = boolean(default=False)
            irreversible_fraction = boolean(default=False)
            dissolved_diameter = float(min=0., max=100.e-6,default=0.)
            particle_diameter = float(min=0., max=100.e-6,default=5.e-6)
            particle_diameter_uncertainty = float(min=0., max=100.e-6, default=1.e-7)
            activity_per_element = float(min=0., max=1.e18, default=1.)
            [[species]]
                LMM                        = boolean(default=False)
                LMMcation                  = boolean(default=False)
                LMManion                   = boolean(default=False)
                Colloid                    = boolean(default=False)
                Humic_colloid              = boolean(default=False)
                Polymer                    = boolean(default=False)
                Particle_reversible        = boolean(default=False)
                Particle_slowly_reversible = boolean(default=False)
                Particle_irreversible      = boolean(default=False)
                Sediment_reversible        = boolean(default=False)
                Sediment_slowly_reversible = boolean(default=False)
                Sediment_irreversible      = boolean(default=False)
            [[transformations]]
                Kd         = float(min=0., max=1.e9, default=0.)
                Dc         = float(min=0., max=1.e6, default=1.16e-5)
                slow_coeff = float(min=0., max=1.e6, default=1.2e-7)
            [[sediment]]
                sedmixdepth        = float(min=0., max=100., default=1.)
                sediment_density   = float(min=0., max=10000., default=2600.)
                effective_fraction = float(min=0., max=1., default=0.9)
                corr_factor        = float(min=0., max=10., default=0.1)
                porosity           = float(min=0., max=1., default=0.6)
                layer_thick        = float(min=0., max=100., default=1.)
                desorption_depth          = float(min=0., max=100., default=1.)
                desorption_depth_uncert   = float(min=0., max=100., default=.5)
                resuspension_depth        = float(min=0., max=100., default=1.)
                resuspension_depth_uncert = float(min=0., max=100., default=.5)
                resuspension_critvel      = float(min=0., max=1., default=.01)
            '''


    def specie_num2name(self,num):
        return self.name_species[num]



    def specie_name2num(self,name):
        num = self.name_species.index(name)
        return num
#



    def __init__(self, *args, **kwargs):

        self._add_configstring(self.configspec_radionuclidedrift)

        # Calling general constructor of parent class
        super(RadionuclideDrift, self).__init__(*args, **kwargs)









    def prepare_run(self):

        self.name_species=[]
        if self.get_config('radionuclide:species:LMM'):
            self.name_species.append('LMM')
        if self.get_config('radionuclide:species:LMMcation'):
            self.name_species.append('LMMcation')
        if self.get_config('radionuclide:species:LMManion'):
            self.name_species.append('LMManion')
        if self.get_config('radionuclide:species:Colloid'):
            self.name_species.append('Colloid')
        if self.get_config('radionuclide:species:Humic_colloid'):
            self.name_species.append('Humic colloid')
        if self.get_config('radionuclide:species:Polymer'):
            self.name_species.append('Polymer')
        if self.get_config('radionuclide:species:Particle_reversible'):
            self.name_species.append('Particle reversible')
        if self.get_config('radionuclide:species:Particle_slowly_reversible'):
            self.name_species.append('Particle slowly reversible')
        if self.get_config('radionuclide:species:Particle_irreversible'):
            self.name_species.append('Particle irreversible')
        if self.get_config('radionuclide:species:Sediment_reversible'):
            self.name_species.append('Sediment reversible')
        if self.get_config('radionuclide:species:Sediment_slowly_reversible'):
            self.name_species.append('Sediment slowly reversible')
        if self.get_config('radionuclide:species:Sediment_irreversible'):
            self.name_species.append('Sediment irreversible')


        if self.get_config('radionuclide:species:Sediment_slowly_reversible') and \
                    self.get_config('radionuclide:species:Particle_slowly_reversible'):
            self.set_config('radionuclide:slowly_fraction', True)
        if self.get_config('radionuclide:species:Sediment_irreversible') and \
                    self.get_config('radionuclide:species:Particle_irreversible'):
            self.set_config('radionuclide:irreversible_fraction', True)


        self.nspecies      = len(self.name_species)
        self.logger.info( 'Number of species: {}'.format(self.nspecies) )
        for i,sp in enumerate(self.name_species):
            self.logger.info( '{:>3} {}'.format( i, sp ) )

        self.init_transfer_rates()




    def init_transfer_rates(self):
        ''' Initialization of background values in the transfer rates 2D array.
        '''

        transfer_setup=self.get_config('radionuclide:transfer_setup')

        self.logger.info( 'transfer setup: %s' % transfer_setup)


        self.transfer_rates = np.zeros([self.nspecies,self.nspecies])
        self.ntransformations = np.zeros([self.nspecies,self.nspecies])

        if transfer_setup == 'Bokna_137Cs':

            self.num_lmm    = self.specie_name2num('LMM')
            self.num_prev   = self.specie_name2num('Particle reversible')
            self.num_srev   = self.specie_name2num('Sediment reversible')
            self.num_psrev  = self.specie_name2num('Particle slowly reversible')
            self.num_ssrev  = self.specie_name2num('Sediment slowly reversible')


            # Values from Simonsen et al (2019a)
            Kd         = self.get_config('radionuclide:transformations:Kd')
            Dc         = self.get_config('radionuclide:transformations:Dc')
            slow_coeff = self.get_config('radionuclide:transformations:slow_coeff')
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
            self.transfer_rates[self.num_srev,self.num_ssrev] = slow_coeff
            self.transfer_rates[self.num_prev,self.num_psrev] = slow_coeff
            self.transfer_rates[self.num_ssrev,self.num_srev] = slow_coeff*.1
            self.transfer_rates[self.num_psrev,self.num_prev] = slow_coeff*.1

        elif transfer_setup=='dummy':
        # Set of dummy values for testing/development

            self.num_lmm   = self.specie_name2num('LMM')
            if self.get_config('radionuclide:species:Colloid'):
                self.num_col = self.specie_name2num('Colloid')
            if self.get_config('radionuclide:species:Particle_reversible'):
                self.num_prev  = self.specie_name2num('Particle reversible')
            if self.get_config('radionuclide:species:Sediment_reversible'):
                self.num_srev  = self.specie_name2num('Sediment reversible')
            if self.get_config('radionuclide:slowly_fraction'):
                self.num_psrev  = self.specie_name2num('Particle slowly reversible')
                self.num_ssrev  = self.specie_name2num('Sediment slowly reversible')
            if self.get_config('radionuclide:irreversible_fraction'):
                self.num_pirrev  = self.specie_name2num('Particle irreversible')
                self.num_sirrev  = self.specie_name2num('Sediment irreversible')


            if self.get_config('radionuclide:species:Particle_reversible'):
                self.transfer_rates[self.num_lmm,self.num_prev] = 5.e-6 #*0.
                self.transfer_rates[self.num_lmm,self.num_srev] = 1.e-5 #*0.
            if self.get_config('radionuclide:species:Sediment_reversible'):
                self.transfer_rates[self.num_srev,self.num_lmm] = 5.e-6

            if self.get_config('radionuclide:slowly_fraction'):
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




        else:
            self.logger.ERROR('No transfer setup available')


        # Set diagonal to 0. (not possible to transform to present specie)
        if len(self.transfer_rates.shape) == 3:
            for ii in range(self.transfer_rates.shape[0]):
                np.fill_diagonal(self.transfer_rates[ii,:,:],0.)
        else:
            np.fill_diagonal(self.transfer_rates,0.)

#         # HACK :
#         self.transfer_rates[:] = 0.
#         print ('\n ###### \n IMPORTANT:: \n transfer rates have been hacked! \n#### \n ')

        self.logger.info('nspecies: %s' % self.nspecies)
        self.logger.info('Transfer rates:\n %s' % self.transfer_rates)








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
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
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

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W*1000*partsize/my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * T0)  # in cm2/s
        d0 = (partsize * 100) - 0.4 * \
            (9.0 * my_w**2 / (100 * g) * DENSw / dr)**(1.0 / 3.0)  # cm
        W2 = 19.0*d0*(0.001*dr)**(2.0/3.0)*(my_w*0.001*DENSw)**(-1.0/3.0)
        # cm/s
        W2 = W2/100.  # back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W


    def update_transfer_rates(self):
        '''Pick out the correct row from transfer_rates for each element. Modify the
        transfer rates according to local environmental conditions '''

        transfer_setup=self.get_config('radionuclide:transfer_setup')
        if transfer_setup == 'Bokna_137Cs' or transfer_setup=='dummy':
            self.elements.transfer_rates1D = self.transfer_rates[self.elements.specie,:]

            if self.get_config('radionuclide:species:Sediment_reversible'):
                # Only LMM radionuclides close to seabed are allowed to interact with sediments
                # minimum height/maximum depth for each particle
                Zmin = -1.*self.environment.sea_floor_depth_below_sea_level
                interaction_thick = self.get_config('radionuclide:sediment:layer_thick')      # thickness of seabed interaction layer (m)
                dist_to_seabed = self.elements.z - Zmin
                self.elements.transfer_rates1D[(self.elements.specie == self.num_lmm) &
                                 (dist_to_seabed > interaction_thick), self.num_srev] = 0.


            if self.get_config('radionuclide:species:Particle_reversible'):
                # Modify particle adsorption according to local particle concentration
                # (LMM -> reversible particles)
                kktmp = self.elements.specie == self.num_lmm
                self.elements.transfer_rates1D[kktmp, self.num_prev] = \
                            self.elements.transfer_rates1D[kktmp, self.num_prev] * \
                            self.environment.conc3[kktmp] / 1.e-3
        #                    self.environment.particle_conc[kktmp] / 1.e-3


        elif transfer_setup=='Sandnesfj_Al':
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
        deltat = self.time_step.seconds             # length of a time step
        phaseshift = np.array(self.num_elements_active()*[False])  # Denotes which trajectory that shall be transformed

        p = 1. - np.exp(-self.elements.transfer_rates1D*deltat)  # Probability for transformation
        psum = np.sum(p,axis=1)

        ran1=np.random.random(self.num_elements_active())

        # Transformation where ran1 < total probability for transformation
        phaseshift[ ran1 < psum ] = True

        self.logger.info('Number of transformations: %s' % sum(phaseshift))
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

        self.logger.debug('old species: %s' % specie_in[phaseshift])
        self.logger.debug('new species: %s' % specie_out[phaseshift])


        for iin in range(self.nspecies):
            for iout in range(self.nspecies):
                self.ntransformations[iin,iout]+=sum((specie_in[phaseshift]==iin) & (specie_out[phaseshift]==iout))

        self.logger.debug('Number of transformations total:\n %s' % self.ntransformations )


        # Update radionuclide properties after transformations
        self.update_radionuclide_diameter(specie_in, specie_out)
        self.sorption_to_sediments(specie_in, specie_out)
        self.desorption_from_sediments(specie_in, specie_out)





    def sorption_to_sediments(self,sp_in=None,sp_out=None):
        '''Update radionuclide properties  when sorption to sediments occurs'''


        # Set z to local sea depth
        if self.get_config('radionuclide:species:LMM'):
            self.elements.z[(sp_out==self.num_srev) & (sp_in==self.num_lmm)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_srev) & (sp_in==self.num_lmm)]
        if self.get_config('radionuclide:species:LMMcation'):
            self.elements.z[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)] = \
                -1.*self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_srev) & (sp_in==self.num_lmmcation)]



    def desorption_from_sediments(self,sp_in=None,sp_out=None):
        '''Update radionuclide properties when desorption from sediments occurs'''

        desorption_depth = self.get_config('radionuclide:sediment:desorption_depth')
        std = self.get_config('radionuclide:sediment:desorption_depth_uncert')


        if self.get_config('radionuclide:species:LMM'):
            self.elements.z[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] = \
                self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] + desorption_depth
            if std > 0:
                self.logger.debug('Adding uncertainty for desorption from sediments: %s m' % std)
                self.elements.z[(sp_out==self.num_lmm) & (sp_in==self.num_srev)] += np.random.normal(
                        0, std, sum((sp_out==self.num_lmm) & (sp_in==self.num_srev)))
        if self.get_config('radionuclide:species:LMMcation'):
            self.elements.z[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] = \
                self.environment.sea_floor_depth_below_sea_level[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] + desorption_depth
            if std > 0:
                self.logger.debug('Adding uncertainty for desorption from sediments: %s m' % std)
                self.elements.z[(sp_out==self.num_lmmcation) & (sp_in==self.num_srev)] += np.random.normal(
                        0, std, sum((sp_out==self.num_lmmcation) & (sp_in==self.num_srev)))





    def update_radionuclide_diameter(self,sp_in=None,sp_out=None):
        '''Update the diameter of the radionuclides when specie is changed'''


        dia_part=self.get_config('radionuclide:particle_diameter')
        dia_diss=self.get_config('radionuclide:dissolved_diameter')


        # Transfer to reversible particles
        self.elements.diameter[(sp_out==self.num_prev) & (sp_in!=self.num_prev)] = dia_part
        std = self.get_config('radionuclide:particle_diameter_uncertainty')
        if std > 0:
            self.logger.debug('Adding uncertainty for particle diameter: %s m' % std)
            self.elements.diameter[(sp_out==self.num_prev) & (sp_in!=self.num_prev)] += np.random.normal(
                    0, std, sum((sp_out==self.num_prev) & (sp_in!=self.num_prev)))

        # Transfer to slowly reversible particles
        if self.get_config('radionuclide:slowly_fraction'):
            self.elements.diameter[(sp_out==self.num_psrev) & (sp_in!=self.num_psrev)] = dia_part
            if std > 0:
                self.logger.debug('Adding uncertainty for slowly rev particle diameter: %s m' % std)
                self.elements.diameter[(sp_out==self.num_psrev) & (sp_in!=self.num_psrev)] += np.random.normal(
                    0, std, sum((sp_out==self.num_psrev) & (sp_in!=self.num_psrev)))

        # Transfer to irreversible particles
        if self.get_config('radionuclide:irreversible_fraction'):
            self.elements.diameter[(sp_out==self.num_pirrev) & (sp_in!=self.num_pirrev)] = dia_part
            if std > 0:
                self.logger.debug('Adding uncertainty for irrev particle diameter: %s m' % std)
                self.elements.diameter[(sp_out==self.num_pirrev) & (sp_in!=self.num_pirrev)] += np.random.normal(
                    0, std, sum((sp_out==self.num_pirrev) & (sp_in!=self.num_pirrev)))

        # Transfer to LMM
        if self.get_config('radionuclide:species:LMM'):
            self.elements.diameter[(sp_out==self.num_lmm) & (sp_in!=self.num_lmm)] = dia_diss
        if self.get_config('radionuclide:species:LMManion'):
            self.elements.diameter[(sp_out==self.num_lmmanion) & (sp_in!=self.num_lmmanion)] = dia_diss
        if self.get_config('radionuclide:species:LMMcation'):
            self.elements.diameter[(sp_out==self.num_lmmcation) & (sp_in!=self.num_lmmcation)] = dia_diss

        # Transfer to colloids
        if self.get_config('radionuclide:species:Colloid'):
            self.elements.diameter[(sp_out==self.num_col) & (sp_in!=self.num_col)] = dia_diss
        if self.get_config('radionuclide:species:Humic_colloid'):
            self.elements.diameter[(sp_out==self.num_humcol) & (sp_in!=self.num_humcol)] = dia_diss
        if self.get_config('radionuclide:species:Polymer'):
            self.elements.diameter[(sp_out==self.num_polymer) & (sp_in!=self.num_polymer)] = dia_diss




    def bottom_interaction(self,Zmin=None):
        ''' Change speciation of radionuclides that reach bottom due to settling.
        particle specie -> sediment specie '''
        if not  ((self.get_config('radionuclide:species:Particle_reversible')) &
                  (self.get_config('radionuclide:species:Sediment_reversible')) or
                  (self.get_config('radionuclide:slowly_fraction')) or
                  (self.get_config('radionuclide:irreversible_fraction'))):
            return

        bottom = np.array(np.where(self.elements.z <= Zmin)[0])
        kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_prev)[0])
        self.elements.specie[bottom[kktmp]] = self.num_srev
        self.ntransformations[self.num_prev,self.num_srev]+=len(kktmp)
        if self.get_config('radionuclide:slowly_fraction'):
            kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_psrev)[0])
            self.elements.specie[bottom[kktmp]] = self.num_ssrev
            self.ntransformations[self.num_psrev,self.num_ssrev]+=len(kktmp)
        if self.get_config('radionuclide:irreversible_fraction'):
            kktmp = np.array(np.where(self.elements.specie[bottom] == self.num_pirrev)[0])
            self.elements.specie[bottom[kktmp]] = self.num_sirrev
            self.ntransformations[self.num_pirrev,self.num_sirrev]+=len(kktmp)



    def resuspension(self):
        """ Simple method to estimate the resuspension of sedimented particles,
        checking whether the current speed near the bottom is above a critical velocity
        Sediment species -> Particle specie
        """
        # Exit function if particles and sediments not are present
        if not  ((self.get_config('radionuclide:species:Particle_reversible')) &
                  (self.get_config('radionuclide:species:Sediment_reversible'))):
            return

        critvel = self.get_config('radionuclide:sediment:resuspension_critvel')
        resusp_depth = self.get_config('radionuclide:sediment:resuspension_depth')
        std = self.get_config('radionuclide:sediment:resuspension_depth_uncert')

        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level
        x_vel = self.environment.x_sea_water_velocity
        y_vel = self.environment.y_sea_water_velocity
        speed = np.sqrt(x_vel*x_vel + y_vel*y_vel)
        bottom = (self.elements.z <= Zmin)

        resusp = ( (bottom) & (speed >= critvel) )
        self.logger.info('Number of resuspended particles: {}'.format(np.sum(resusp)))

        self.elements.z[resusp] = Zmin[resusp] + resusp_depth
        if std > 0:
            self.logger.debug('Adding uncertainty for resuspension from sediments: %s m' % std)
            self.elements.z[resusp] += np.random.normal(
                        0, std, sum(resusp))
        self.elements.z[resusp] = [min(0,zz) for zz in self.elements.z[resusp]]
        self.ntransformations[self.num_srev,self.num_prev]+=sum((resusp) & (self.elements.specie==self.num_srev))
        self.elements.specie[(resusp) & (self.elements.specie==self.num_srev)] = self.num_prev
        if self.get_config('radionuclide:slowly_fraction'):
            self.ntransformations[self.num_ssrev,self.num_psrev]+=sum((resusp) & (self.elements.specie==self.num_ssrev))
            self.elements.specie[(resusp) & (self.elements.specie==self.num_ssrev)] = self.num_psrev
        if self.get_config('radionuclide:irreversible_fraction'):
            self.ntransformations[self.num_sirrev,self.num_pirrev]+=sum((resusp) & (self.elements.specie==self.num_sirrev))
            self.elements.specie[(resusp) & (self.elements.specie==self.num_sirrev)] = self.num_pirrev






    def update(self):
        """Update positions and properties of radionuclide particles."""

        # Workaround due to conversion of datatype
        self.elements.specie = self.elements.specie.astype(np.int32)

        # Radionuclide speciation
        self.update_transfer_rates()
        self.update_speciation()


        # Turbulent Mixing
        z_before = self.elements.z.copy()
        self.update_terminal_velocity()
        self.vertical_mixing()


        # Resuspension
        self.resuspension()
        self.logger.info('Speciation: {} {}'.format([sum(self.elements.specie==ii) for ii in range(self.nspecies)],self.name_species))



        # Horizontal advection
        lon, lat = self.elements.lon, self.elements.lat
        self.advect_ocean_current()

        # Reset lat lon for sedimented trajectories (reject hor. advection)
        if self.get_config('radionuclide:species:Sediment_reversible'):
            self.elements.lon[self.elements.specie==self.num_srev]   = lon[self.elements.specie==self.num_srev]
            self.elements.lat[self.elements.specie==self.num_srev]   = lat[self.elements.specie==self.num_srev]
        if self.get_config('radionuclide:slowly_fraction'):
            self.elements.lon[self.elements.specie==self.num_ssrev] = lon[self.elements.specie==self.num_ssrev]
            self.elements.lat[self.elements.specie==self.num_ssrev] = lat[self.elements.specie==self.num_ssrev]
        if self.get_config('radionuclide:irreversible_fraction'):
            self.elements.lon[self.elements.specie==self.num_sirrev] = lon[self.elements.specie==self.num_sirrev]
            self.elements.lat[self.elements.specie==self.num_sirrev] = lat[self.elements.specie==self.num_sirrev]

        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()

        # Reset z for sedimented trajectories (reject vertical advection and mixing)
        if self.get_config('radionuclide:species:Sediment_reversible'):
            self.elements.z[self.elements.specie==self.num_srev]   = z_before[self.elements.specie==self.num_srev]
        if self.get_config('radionuclide:slowly_fraction'):
            self.elements.z[self.elements.specie==self.num_ssrev] = z_before[self.elements.specie==self.num_ssrev]
        if self.get_config('radionuclide:irreversible_fraction'):
            self.elements.z[self.elements.specie==self.num_sirrev] = z_before[self.elements.specie==self.num_sirrev]










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
                                              ):
        '''Write netCDF file with map of radionuclide species densities and concentrations'''

        from netCDF4 import Dataset, date2num #, stringtochar

        self.logger.info('Postprocessing: Write density and concentration to netcdf file')

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



        if activity_unit==None:
            activity_unit='Bq'  # default unit for radionuclides

        activity_per_element = self.get_config('radionuclide:activity_per_element')


        z = self.get_property('z')[0]
        if not zlevels==None:
            zlevels = np.sort(zlevels)
            z_array = np.append(np.append(-10000, zlevels) , max(0,np.nanmax(z)))
        else:
            z_array = [min(-10000,np.nanmin(z)), max(0,np.nanmax(z))]
        self.logger.info('z_array: {}'.format(  [str(item) for item in z_array] ) )




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
            self.logger.info('H.shape: ' + str(H.shape))
            Hsm = np.zeros_like(H)
            for zi in range(len(z_array)-1):
                for sp in range(self.nspecies):
                    for ti in range(H.shape[0]):
                        Hsm[ti,sp,zi,:,:] = self.horizontal_smooth(H[ti,sp,zi,:,:],n=smoothing_cells)



        # Convert from density to concentration
        self.logger.info('Activity: '+str(activity_per_element)+' '+ activity_unit+ ' per unit')

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
                ndt = int( deltat / (mdt.seconds/3600.) )
            times2 = times[::ndt]
            times2 = times2[1:]
            odt = int(cshape[0]/ndt)
            self.logger.info ('ndt '+ str(ndt))   # number of time steps over which to average in conc file
            self.logger.info ('odt '+ str(odt))   # number of average slices


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
        self.logger.info('Wrote to '+filename)





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

