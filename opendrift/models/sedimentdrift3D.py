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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# Copyright 2018, Simon, MetOcean Solutions Ltd.

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift import Lagrangian3DArray

class SedimentElement(Lagrangian3DArray):
    # Lagrangian3DArray has already the variables terminal_velocity, and wind_drift_factor
    variables = Lagrangian3DArray.add_variables([
        ('settled', {'dtype': np.int16,  # 0 is active, 1 is settled
                     'units': '1',
                     'default': 0}),
        ('age_seconds', {'dtype': np.float32,
                 'units': 's',
                 'default': 0}),
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': -0.001}), # 1 mm/s negative buoyancy
        ('critical_shear_stress', {'dtype': np.float32,
                                  'units': 'N/m2',
                                  'default': 0.0}),
        ('d50', {'dtype': np.float32,  # median grain size
                     'units': 'm',
                     'default': 125.0*1e-6}),
        ('density', {'dtype': np.float32, # grain density = mass_solid/volume_solid
                     'units': 'kg/m^3',
                     'default': 2650.}),       
        ])

class SedimentDrift3D(OceanDrift): # based on OceanDrift base class
    """Trajectory model based on the OpenDrift framework using the OceanDrift baseclass

    Sediment 3D motion 
    Propagation with horizontal and vertical ocean currents, horizontal and 
    vertical diffusions (additional wind drag inherited from base class but probably not relevant here).
    Suitable for sediment tracers, e.g. for tracking sediment particles.
    
    Authors:
    Simon Weppe - MetOcean Solutions, NZ
    Remy Zyngfogel - Calypso Science, NZ

    """
    ElementType = SedimentElement 

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height' : {'fallback': 0},
        'sea_surface_wave_period_at_variance_spectral_density_maximum': {'fallback': 0},
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
        'sea_floor_depth_below_sea_level': {'fallback': 0},
        # 'sea_water_temperature': {'fallback': 15., 'profiles': True},
        # 'sea_water_salinity': {'fallback': 35., 'profiles': True},
        'sea_water_temperature': {'fallback': 15.},
        'sea_water_salinity': {'fallback': 35.},
        }
        
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-20, 0]

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(SedimentDrift3D, self).__init__(*args, **kwargs)
        
        # resuspension (switched off by default)
        self._add_config({
            'drift:resuspension': {'type': 'bool', 'default': False,
                # 'min': 0, 'max': 1e10, 'units': '-',
                'description': 'switch to activate/deactivate resuspension',
                'level': self.CONFIG_LEVEL_ESSENTIAL},

            'processes:update_terminal_velocity': {'type': 'bool', 'default': False,
                'description': 'switch to activate "online" terminal velocity calculation based on ambient settings at each time step.\
                 If set to False, terminal_velocity is kept constant with value provided at seeding.',
                'level': self.CONFIG_LEVEL_ADVANCED},

            'processes:terminal_velocity_method': {'type': 'enum', 'enum': ['dietrich','stokes','zhiyao','ahrens','komar','van_rijn1984','soulsby1997'], 'default': 'zhiyao',
                'description': 'selection of equations used for calculating the particle''s terminal_velocity',
                'level': self.CONFIG_LEVEL_ADVANCED},
                })

        # By default, sediments do strand at coastline
        self._set_config_default('general:coastline_action','stranding')
        # Vertical mixing is enabled as default
        self._set_config_default('drift:vertical_mixing', True)
        # Vertical advection switched off by default (if w is available)
        self._set_config_default('drift:vertical_advection', False)
        #Stokes drift probably not relevant here, expect maybe for slow-settling, surface sediment
        # keeping here as place holder
        self._set_config_default('drift:stokes_drift',False)

        # Settling on seafloor : if no resuspension (default) : deactivate settled particles 
        #                        if resuspension is on, this will be set to 'lift_to_seafloor'
        #                        other options : ['none', 'lift_to_seafloor', 'deactivate', 'previous']
        self._set_config_default('general:seafloor_action', 'deactivate') # deactivate by default. Particles will be labelled as 'seafloor'

        self.max_speed = 5.0

        # set some default values for specific parameters used in some terminal_velocity equations
        self.factor = 1 
        self.powers_roundness = 6

    def update_terminal_velocity(self, Tprofiles=None, Sprofiles=None,
                                 z_index=None):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables. Sub-modules should overload
        this method for particle-specific behaviour
        Code for terminal_velocity computations translated from R from https://github.com/Davidatlarge/sinking_velocity
        """

        if self.get_config('processes:update_terminal_velocity'): 
            logger.debug('processes:update_terminal_velocity = True : update sediment particles'' terminal_velocity' )
            g = 9.81  # ms-2
            factor=self.factor
            powers_roundness=self.powers_roundness 
            d50 = self.elements.d50  # NB: r is diameter, not radius
            method = self.get_config('processes:terminal_velocity_method')
            logger.debug('    using terminal_velocity_method = %s',method)
            
            # Prepare interpolation of temp, salt
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
            
            # Do interpolation of temp, salt if profiles were passed into
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

            # define some parameters required for computation of terminal_velocity
            pdensity = self.elements.density
            rho_water = self.sea_water_density(T=T0, S=S0) # from opendrift built-in function...could be sufficient for use in equations ? 
            part_depth = np.abs(self.elements.z) # particle depth below surface
            latitude = self.elements.lat 
            # water density at particle depths 
            water_density=sw_dens(s = S0, t = T0, p = np.array(d2p(part_depth, lat = latitude)/10+1)) # value of d2p is dbar of pressure exerted by water, without air, hence /10 and +1 ; method="Gibbs" by default  
            #viscosity at particle depths          
            dynamic_viscosity = viscosity(S = S0, t = T0, P = d2p(part_depth, lat = latitude)/10+1) # viscosity in centipoise (cP); 1 cP = 0.001 kg·m−1·s−1
            dynamic_viscosity = dynamic_viscosity / 1000 # factor 1e3 to convert viscosity to [kg·m−1·s−1]
            kinematic_viscosity = dynamic_viscosity / water_density 

            # compute terminal_velocity using correct method
            if method=="dietrich": # Dietrich 1982 (Water Resour. Res.)
                sv=dietrich(pdensity,water_density,g,d50,kinematic_viscosity,factor,powers_roundness)
            elif method=="stokes":
                sv=stokes(d50,g,pdensity,water_density,dynamic_viscosity)
            elif method == "zhiyao": # Zhiyao et al. 2008 (Water Science and Engineering)
                sv=zhiyao(pdensity,water_density,g,kinematic_viscosity,d50)
            elif method == "ahrens":
                sv=ahrens(pdensity,water_density,g,kinematic_viscosity,d50)
            elif method == "komar":
                sv=komar(dynamic_viscosity,pdensity,water_density,g,d50,factor)
            elif method == "soulsby1997":
                sv=soulsby1997(pdensity,water_density,g,kinematic_viscosity,d50)
            elif method == "van_rijn1984":
                sv=van_rijn1984(pdensity,water_density,g,kinematic_viscosity,d50)

            self.elements.terminal_velocity = sv*-1 # for negative buoyancy = settling, terminal_velocity<0
            logger.debug('    %.5f [m] =< d50 <= %.5f [m]' % (np.min(self.elements.d50),np.max(self.elements.d50)))
            logger.debug('    %.5f [m/s] =< terminal_velocity <=  %.5f [m/s]' % (np.min(self.elements.terminal_velocity),np.max(self.elements.terminal_velocity)))
        else:
            logger.debug('processes:update_terminal_velocity = False : using constant terminal_velocity set in seed_elements(), or default one' )
            logger.debug('    %.5f [m/s] =< terminal_velocity <=  %.5f [m/s]' % (np.min(self.elements.terminal_velocity),np.max(self.elements.terminal_velocity)))
            pass 

    def update(self):
        """Update positions and properties of elements."""
        self.elements.age_seconds += self.time_step.total_seconds()

        # Simply move particles with ambient current
        self.advect_ocean_current() # from physics_methods.py

        # Advect particles due to wind drag
        # (according to specified wind_drift_factor)
        self.advect_wind()

        # Stokes drift - deactivated for now
        self.stokes_drift()
        
        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()  
            self.vertical_mixing()
        else:  # Buoyancy
            self.update_terminal_velocity() 
            self.vertical_buoyancy()

        # Vertical advection
        self.vertical_advection()

        # Sediment resuspension physics
        self.sediment_resuspension() #-

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')
        
    def sediment_resuspension(self):
        """
        Compute ambient bed shear stresses at particle positions and determine if particles that settled 
        can be resuspensed or not based on their critical_shear_stress

        The function uses the element attribute self.elements.moving to be able to freeze/unfreeze particles 
        motions while still keeping them "active" (see in advect_ocean_current())

        Note the self.elements.moving is set to 0 in bottom_interaction(), called within vertical_mixing()
        when particles touch the seafloor
        """

        if self.get_config('drift:resuspension') is True:
            self.set_config('general:seafloor_action', 'lift_to_seafloor')
            logger.debug('Sediemnt resuspension physics included : drift:resuspension == True')
            # 1-compute bed shear stresses at particle locations
            tau_cw,tau_cw_max = self.bedshearstress_current_wave()
            # 2-compare ambient max bedshearstress with critical bed shearstress
            # Note the self.elements.moving is set to 0 in bottom_interaction(), called within vertical_mixing()
            # when particles touch the seafloor
            to_resuspend = (np.logical_and(self.elements.moving == 0,tau_cw_max>self.elements.critical_shear_stress))
            if np.sum(to_resuspend) > 0 :
                logger.debug('Resuspending %s elements where tau_cw_max>critical_shear_stress' % np.sum(to_resuspend))
                sea_floor_depth = self.sea_floor_depth()
                # Resuspend 1 cm above seafloor
                self.elements.z[to_resuspend] = -sea_floor_depth[to_resuspend] + 0.01
                # Allow moving again
                self.elements.moving[to_resuspend] = 1
            else:              
              logger.debug('No elements to resuspend (tau_cw_max < critical_shear_stress everywhere')
 
    def bottom_interaction(self, seafloor_depth):
        """Sub method of vertical_mixing, determines settling"""
        # Elements at or below seafloor are settled, by setting self.elements.moving to 0.
        # These elements will not move until possible later resuspension.

        sea_floor_depth =  self.sea_floor_depth()
        below = np.where(self.elements.z < -sea_floor_depth)[0]
        self.elements.z[below] = -sea_floor_depth[below]
        settling = np.logical_and(self.elements.z <= seafloor_depth, self.elements.moving==1)
        if np.sum(settling) > 0:
            logger.debug('Settling %s elements at seafloor' % np.sum(settling))
            self.elements.moving[settling] = 0

#######################################################################################################################
# General physics functions
# >> could be moved to physics_methods.py once cross-checked / accepted 
#######################################################################################################################

    def bedshearstress_current_wave(self):
        """
        Computation of bed shear stress due to current and waves
        current-related stress is computed following a drag-coefficient approach
        wave-related stress is computed following Van Rijn approach
        Combined wave-current mean and max stresses are computed following Soulsby(1995) approach

        https://odnature.naturalsciences.be/coherens/manual#manual
        https://odnature.naturalsciences.be/downloads/coherens/documentation/chapter7.pdf#nameddest=Bed_shear_stresses
        
        http://www.coastalwiki.org/wiki/Shallow-water_wave_theory#
        http://www.coastalwiki.org/wiki/Shallow-water_wave_theory#Seabed_Friction  

        General relationships obtained from :
        https://repository.tudelft.nl/islandora/object/uuid%3Aea12eb20-aee3-4f58-99fb-ebc216e98879
        Description of TRANSPOR2004 and Implementation in Delft3D-ONLINE

        see also :
        https://content.oss.deltares.nl/delft3d/manuals/D-Water_Quality_Processes_Technical_Reference_Manual.pdf
        >> Section 13.10
        """

        rho_water = 1027 # kg/m3
        z0 = 0.001 # roughness height - hard-coded constant for now 
        water_depth = np.abs(self.sea_floor_depth()) # water depth positive down
        current_speed = self.current_speed() # depth-averaged current 

        #######################################################
        # current-related bed shear stress
        #######################################################

        # depth-averaged current approach :
        if True : # current data fron reader is depth-averaged
            Cdrag=( 0.4 /(np.log(abs(water_depth/z0))-1) )**2
            #Now compute the bed shear stress [N/m2] 
            tau_cur=rho_water*Cdrag*current_speed**2 # eq. 7.1 in COHERENS Manual
        else:
            # 3D currents - to implement
            last_wet_bin_depth = 0.0 
            Cdrag=( 0.4 /(np.log(abs(last_wet_bin_depth/z0))-1) )**2
            #Now compute the bed shear stress [N/m2] 
            tau_cur=rho_water*Cdrag*current_speed**2        

        #######################################################
        # wave-related bed shear stress (if wave available)
        #######################################################
        hs = self.significant_wave_height()
        tp = self.wave_period()
        # wave-related roughness

        # see VanRijn 
        # https://tinyurl.com/nyjcss5w
        # SIMPLE GENERAL FORMULAE FOR SAND TRANSPORT IN RIVERS, ESTUARIES AND COASTAL WATERS
        # >> page 6
        # 
        # Note : VanRijn 2007 suggests same equations than for current-related roughness 
        # where 20*d50 <ksw<150*d50: Here we are using Nikuradse roughness for consistency 
        # with the use of z0 in the current-related shear stress 

        ksw=30*z0 # wave related bed roughness - taken as Nikuradse roughness 
        w=2*np.pi/tp # angular frequency
        kh = qkhfs( w, water_depth ) # from dispersion relationship 
        Adelta = hs/(2*np.sinh(kh)) # peak wave orbital excursion 
        Udelta = (np.pi*hs)/(tp*np.sinh(kh))  # peak wave orbital velocity linear theory 
        # wave-related friction coefficient (Swart,1974) and eq. 3.8 on VanRijn pdf
        # see also COHERENS manual eq. 7.17 which is equivalent since exp(a+b) =exp(a)*exp(b)
        fw_swart = np.exp(-5.977+5.213*(Adelta/ksw)**-0.194)  
        fw_swart = np.minimum(fw_swart,0.3)
        fw_soulsby = 0.237 * (Adelta/ksw)**-0.52 #eq. 7.18 COHERENS, not used for now

        tau_wave = 0.25 * rho_water * fw_swart * (Udelta)**2 # wave-related bed shear stress eq. 3.7 on VanRijn pdf
        #cycle mean bed shear stress according to Soulsby,1995, see also COHERENS manual eq. 7.14
        tau_cw=tau_cur*[1+1.2*(tau_wave/(tau_cur+tau_wave))**3.2]
        # max bed shear stress during wave cycle >> used for the resuspension criterion.
        theta_cur_dir = 0.0 #angle between direction of travel of wave and current, in radians, in practice rarely known so assume 0.0 for now
        # tau_max = ( (tau_cur + tau_wave*np.cos(theta_cur_dir))**2 + (tau_wave*np.sin(theta_cur_dir))**2 )**0.5 
        tau_cw_max = ( tau_cur**2 + tau_wave**2 + 2*tau_cur*tau_wave*np.cos(theta_cur_dir) )**0.5 # COHERENS eq. 7.15
        
        return tau_cw[0],tau_cw_max 


        # For more reference on bedshearstress and friction factor
        # >> https://www.researchgate.net/publication/233398452_A_fully_coupled_3D_wave-current_interaction_model_on_unstructured_grids
        # or even see SCHISM code ?
        # p4/18


#from  https://github.com/csherwood-usgs/crspy/blob/master/crspy.py
def qkhfs( w, h ):
    """
    Quick iterative calculation of kh in gravity-wave dispersion relationship
    kh = qkhfs(w, h )
    
    Input
        w - angular wave frequency = 2*pi/T where T = wave period [1/s]
        h - water depth [m]
    Returns
        kh - wavenumber * depth [ ]
    Orbital velocities from kh are accurate to 3e-12 !
    RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    """
    g = 9.81
    x = w**2.0 *h/g
    y = np.sqrt(x) * (x<1.) + x *(x>=1.)
    # is this faster than a loop?
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    kh = y
    return kh

#######################################################################################################################
# methods for computing terminal_velocity of particles
#######################################################################################################################

def  stokes(pdiameter,g,pdensity,wdensity,dynamic_viscosity):
    # ... according to Stokes' Law
    if(pdiameter > 0.0002):
        print("Particle diameter > 200 µm! \n Stokes' Law will overestimate sinking velocity! \n Use another method!")
    sinking = (pdiameter**2*g*(pdensity-wdensity)) / (18*dynamic_viscosity)
    if sinking<0:
        sinking = np.nan
    return sinking

def dietrich(pdensity,wdensity,g,pdiameter,kinematic_viscosity,factor,p):
    # ... according to Dietrich 1982
    # args :
    #   pdensity : grain density [kg/m3]
    #   wdensity : seawater density [kg/m3]
    #   g : 9.81 #  - acceleration due to gravity [m.s-2]
    #   kinematic_viscosity : kinematic viscosity of water  [m2 s-1 ]
    #   pdiameter : grain diameter [m]
    #   factor : used in equation, 1 by default
    #   p : powers_roundness 6 by default
    Dstar = ((pdensity - wdensity) * g * pdiameter**3) / (wdensity * kinematic_viscosity**2) # dimensionless size D*; introduces NaN for particles with lower density than water
    R1 = -3.76715 + 1.92944*(np.log10(Dstar)) - 0.09815*\
        (np.log10(Dstar))**2 - 0.00575*(np.log10(Dstar))**3 + 0.00056*(np.log10(Dstar))**4 # size and denisty effect
    R2 = (np.log10(1-((1-factor)/0.85))) - (1-factor)**2.3*np.tanh(np.log10(Dstar)-4.6) +\
         0.3*(0.5-factor)*(1-factor)**2 * (np.log10(Dstar)-4.6) # shape effect
    R3 = (0.65-((factor/2.83) * np.tanh(np.log10(Dstar)-4.6)))**(1+(3.5-p)/2.5) # roundness effect
    Wstar = R3 * 10**(R1+R2) 
    sinking = ((Wstar*(pdensity-wdensity) * g*kinematic_viscosity) / wdensity)**(1/3) # introduced NaN for those particles that have no value for ESD because the dimensions were not measured
    return sinking

def zhiyao(pdensity,wdensity,g,kinematic_viscosity,pdiameter):
    # ... according to Zhiyao et al. 2008 (Water Science and Engineering)
    # args :
    #   pdensity : grain density [kg/m3]
    #   wdensity : seawater density [kg/m3]
    #   g : 9.81 #  - acceleration due to gravity [m.s-2]
    #   kinematic_viscosity : kinematic viscosity of water  [m2 s-1 ]
    #   pdiameter : grain diameter [m]
    delta = pdensity/wdensity-1 
    Dstar = ((delta*g)/kinematic_viscosity**2)**(1/3)*pdiameter # formula (5)
    sinking = (kinematic_viscosity/pdiameter)*Dstar**3 * (38.1+0.93*Dstar**(12/7))**(-7/8) # formula (11)
    return sinking

def ahrens(pdensity,wdensity,g,kinematic_viscosity,pdiameter):
    # args :
    #   pdensity : grain density [kg/m3]
    #   wdensity : seawater density [kg/m3]
    #   g : 9.81 #  - acceleration due to gravity [m.s-2]
    #   kinematic_viscosity : kinematic viscosity of water  [m2 s-1 ]
    #   pdiameter : grain diameter [m]

    Delta = (pdensity - wdensity) / wdensity
    particle_diameter_cm = pdiameter * 100
    kinematic_viscosity_cm2_s = kinematic_viscosity * 10000
    gravity_cm_s2 = g * 100
               
    A = Delta * gravity_cm_s2 * particle_diameter_cm**3 / kinematic_viscosity_cm2_s**2
    C1 = 0.055 * np.tanh(12*A**-0.59 * np.exp(-0.0004*A))
    Ct = 1.06 * np.tanh(0.016*A**0.50 * np.exp(-120/A))
        
    sinking_velocity_cm_sec = C1 * Delta * gravity_cm_s2 * particle_diameter_cm**2\
                             / kinematic_viscosity_cm2_s + Ct *\
                             np.sqrt(Delta * gravity_cm_s2 * particle_diameter_cm) # term associated with turbulent flow
    return sinking_velocity_cm_sec / 100

def komar(dynamic_viscosity,pdensity,wdensity,g,pdiameter,factor):
    # ... velocity for ellipsoid particles according to Komar 1980 (equation 2 in abstract)
    # args :
    #   pdensity : grain density [kg/m3]
    #   wdensity : seawater density [kg/m3]
    #   g : 9.81 #  - acceleration due to gravity [m.s-2]
    #   kinematic_viscosity : kinematic viscosity of water  [m2 s-1 ]
    #   pdiameter : grain diameter [m]
    #   factor : used in equation, 1 by default

    sinking = (1/18) * (1/dynamic_viscosity) * \
                  (pdensity - wdensity) * \
              g * pdiameter**2 * factor**0.380 
    return sinking

def soulsby1997(pdensity,wdensity,g,kinematic_viscosity,pdiameter):
    # args :
    #   pdensity : grain density [kg/m3]
    #   wdensity : seawater density [kg/m3]
    #   g : 9.81 #  - acceleration due to gravity [m.s-2]
    #   kinematic_viscosity : kinematic viscosity of water  [m2 s-1 ]
    #   pdiameter : grain diameter [m]

    v = kinematic_viscosity # m2 s-1 - kinematic viscosity of water
    ps = pdensity # kg m-3 - grain density 
    p = wdensity # kg m-3 - water density
    s = ps/p # ratio of densities of grain and water
     # dimensionless grain size calculation - Eq 98, Soulsby 1997
    d = pdiameter
    D_star = (((g*(s-1))/(v**2))**(1/3))*d 
    v_d = v/d 
    #settling velocity calculation - Eq 103, Soulsby 1997
    ws = v_d*(((10.36**2 + 1.049*D_star**3)**(1/2)) - 10.36) 
    return ws

def van_rijn1984(pdensity,wdensity,g,kinematic_viscosity,pdiameter):
    # settling velocity calculation for Van Rijn 1984 - equations 101a,b,c Soulsby 1997
    # 
    # args :
    #   pdensity : grain density [kg/m3]
    #   wdensity : seawater density [kg/m3]
    #   g : 9.81 #  - acceleration due to gravity [m.s-2]
    #   kinematic_viscosity : kinematic viscosity of water  [m2 s-1 ]
    #   pdiameter : grain diameter [m]
    # 
    v = kinematic_viscosity # - kinematic viscosity of water
    ps = pdensity # kg m-3 - grain density
    p = wdensity # kg m-3 - water density
    s = ps/p # ratio of densities of grain and water

    # dimensionless grain size calculation - Eq 98, Soulsby 1997
    d = pdiameter
    D_star = (((g*(s-1))/(v**2))**(1/3))*d 
    D_cube = D_star**3
    ws = np.nan*np.ones(d.shape[0])
    # if D_cube <= 16.187:
    #     ws = (v*D_cube)/(18*d)
    # elif D_cube > 16.187 and D_cube <= 16187:
    #     ws = ((10.*v)/d)*(((1+0.01*D_cube)**(1/2))-1)
    # elif D_cube > 16187:
    #     ws = (1.1*v*(D_star**(1.5)))/d
    id1 = np.where(D_cube <= 16.187)
    ws[id1] = (v[id1]*D_cube[id1])/(18*d[id1])

    id2 = np.where((D_cube > 16.187) & (D_cube <= 16187))
    ws[id2] = ((10.*v[id2])/d[id2])*(((1+0.01*D_cube[id2])**(1/2))-1)

    id3 = np.where(D_cube > 16187)
    ws[id3] = (1.1*v[id3]*(D_star[id3]**(1.5)))/d[id3]
    
    if (ws == np.nan).any():
        import pdb;pdb.set_trace()
        print('still some nans in van_rijn1984 terminal_velocity computation - check')
    return ws

#######################################################################################################################
# UNESCO equations for seawater density and viscosity
# 
# References
# ----------
# .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
#    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
#    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
#    http://unesdoc.UNESCO.org/images/0005/000598/059832eb.pdf
# 
# Consider updating to Thermodynamic Equation of Seawater 2010 (TEOS-10) :
# 
# https://github.com/TEOS-10/GSW-Python
#######################################################################################################################

def d2p(depth,lat):

    ## Checking for lat positive
    lat = np.abs(lat)
    plat = np.abs(lat*np.pi)/180 
    d = np.sin(plat)
    c1 = 5.92e-3 + d**2 * 5.25e-3

    pressure = ((1-c1)-np.sqrt(((1-c1)**2)-(8.84e-6*depth))) / 4.42e-6

    return pressure

def sw_smow(t):
    """
    Calculate the density of Standard Mean Ocean Water (pure water).
    Parameters
    ----------
    t : ndarray
        Temperature (1D array) in degrees Celsius.
    Returns
    -------
    rho : ndarray
        Density in kg m^{-3}.
    """
    c68 = 1.00024   # conversion constant to T68 temperature scale.

    # Coefficients
    a0 = 999.842594
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9

    T68 = t * c68

    dens = a0 + (a1 * T68) + (a2 * T68**2) + (a3 * T68**3) \
            + (a4 * T68**4) + (a5 * T68**5)

    return dens

def sw_dens0(t, s):
    """
    Calculate sea water density at atmospheric surface pressure.
    Parameters
    ----------
    t : ndarray
        Temperature (1D array) in degrees Celsius.
    s: ndarray
        Salinity (PSU). Must be the same size as t.
    Returns
    -------
    dens : ndarray
        Seawater density at atmospheric surface pressure (kg m^{-1}).
    """
    c68 = 1.00024   # conversion constant to T68 temperature scale.
    b0 = 8.24493e-1
    b1 = -4.0899e-3
    b2 = 7.6438e-5
    b3 = -8.2467e-7
    b4 = 5.3875e-9

    c0 = -5.72466e-3
    c1 = 1.0227e-4
    c2 = -1.6546e-6

    d0 = 4.8314e-4

    t68 = t * c68

    dens = s * (b0 + (b1 * t68) + (b2 * t68**2) + (b3 * t68**3) + (b4 * t68**4)) + \
            s**1.5 * (c0 + (c1 * t68) + (c2 * t68**2)) + (d0 * s**2)

    dens = dens + sw_smow(t68)

    return dens


def sw_seck(t, s, p):
    """
    Calculate Secant Bulk Modulus (K) of seawater.
    Parameters
    ----------
    t : ndarray
        Temperature (1D array) in degrees Celsius.
    s : ndarray
        Salinity (1D array) in practical salinity units (unitless). Must be the
        same shape as t.
    p : ndarray
        Pressure (1D array) in decibars. Must be the same shape as t.
    Returns
    -------
    k : ndarray
        Secant Bulk Modulus of seawater.
    """

    # Compression terms
    c68 = 1.00024   # conversion constant to T68 temperature scale.
    T68 = t * c68
    Patm = p / 10.0  # convert to bar

    h3 = -5.77905e-7
    h2 = 1.16092e-4
    h1 = 1.43713e-3
    h0 = 3.239908

    AW = h0 + (h1 * T68) + (h2 * T68**2) + (h3 * T68**3)

    k2 = 5.2787e-8
    k1 = -6.12293e-6
    k0 = 8.50935e-5

    BW = k0 + (k1 + k2 * T68) * T68

    e4 = -5.155288e-5
    e3 = 1.360477e-2
    e2 = -2.327105
    e1 = 148.4206
    e0 = 19652.21

    KW = e0 + (e1 + (e2 + (e3 + e4 * T68) * T68) * T68) * T68

    # K at atmospheric pressure

    j0 = 1.91075e-4

    i2 = -1.6078e-6
    i1 = -1.0981e-5
    i0 = 2.2838e-3

    A = AW + s * (i0 + (i1 * T68) + (i2 * T68**2)) + (j0 * s**1.5)

    m2 = 9.1697e-10
    m1 = 2.0816e-8
    m0 = -9.9348e-7

    # Equation 18
    B = BW + (m0 + (m1 * T68) + (m2 * T68**2)) * s

    f3 = -6.1670e-5
    f2 = 1.09987e-2
    f1 = -0.603459
    f0 = 54.6746

    g2 = -5.3009e-4
    g1 = 1.6483e-2
    g0 = 7.944e-2

    # Equation 16
    K0 = KW + s * (f0 + (f1 * T68) + (f2 * T68**2) + (f3 * T68**3)) + \
            s**1.5 * (g0 + (g1 * T68) + (g2 * T68**2))

    # K at t, s, p
    K = K0 + (A * Patm) + (B * Patm**2)  # Equation 15

    return K

def sw_dens(t, s, p):
    """
    Convert temperature, salinity and pressure to density.
    Parameters
    ----------
    t : ndarray
        Temperature (1D array) in degrees Celsius.
    s : ndarray
        Salinity (1D array) in practical salinity units (unitless). Must be the
        same shape as t.
    p : ndarray
        Pressure (1D array) in decibars. Must be the same shape as t.
    Returns
    -------
    rho : ndarray
        Density in kg m^{-3}.
    Notes
    -----
    Valid temperature range is -2 to 40C, salinity is 0-42 and pressure is
    0-10000 decibars. Warnings are issued if the data fall outside these
    ranges.
    """

    # Check for values outside the valid ranges.
    if t.min() < -2:
        n = np.sum(t < -2)
        print('WARNING: {} values below minimum value temperature (-2C)'.format(n))

    if t.max() > 40:
        n = np.sum(t > 40)
        print('WARNING: {} values above maximum value temperature (40C)'.format(n))

    if s.min() < 0:
        n = np.sum(s < 0)
        print('WARNING: {} values below minimum salinity value (0 PSU)'.format(n))

    if s.max() > 42:
        n = np.sum(s > 42)
        print('WARNING: {} values above maximum salinity value (42C)'.format(n))

    if p.min() < 0:
        n = np.sum(p < 0)
        print('WARNING: {} values below minimum pressure value (0 decibar)'.format(n))

    if p.max() > 10000:
        n = np.sum(p > 10000)
        print('WARNING: {} values above maximum pressure value (10000 decibar)'.format(n))

    dens0 = sw_dens0(t, s)
    k = sw_seck(t, s, p)
    Patm = p / 10.0  # pressure in bars
    rho = dens0 / (1 - Patm / k)

    return rho

def gravity(lat,method="UNESCO"):
    X = np.sin(lat * np.pi / 180.)
    X = X * X
    X2 = (np.sin(2 * lat * np.pi / 180.))**2
    if (method == "UNESCO"):
        grav =  9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * X) *X)
    else:
        grav =  9.780327*(1.0 + 0.0053024 * X -0.0000058 * X2)
  
    return grav

def viscosity(S,t,P):
#viscosity <- function (S = 35, t = 25, P = 1.013253) {

    vi=1.7910 - t*(6.144e-02 - t*(1.4510e-03 - t*1.6826e-05))+\
      - 1.5290e-04*P + 8.3885e-08*P*P + 2.4727e-03*S+\
      + (6.0574e-06*P - 2.6760e-09*P*P)*t + (t*(4.8429e-05+\
      - t*(4.7172e-06 - t*7.5986e-08)))*S

    return vi