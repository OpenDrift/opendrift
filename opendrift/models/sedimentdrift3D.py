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
import logging
from opendrift.models.opendrift3D import OpenDrift3DSimulation
from opendrift.models.oceandrift import OceanDrift
# from opendrift.elements.passivetracer import PassiveTracer
from opendrift.elements.buoyanttracer import BuoyantTracer
# from opendrift.elements.sedimenttracer import SedimentTracer


class SedimentDrift3D(OpenDrift3DSimulation): # based on OpenDrift3DSimulation base class
    """Trajectory model based on the OpenDrift framework using the OpenDrift3DSimulation baseclass

    Sediment 3D motion 
    Propagation with horizontal and vertical ocean currents, horizontal and 
    vertical diffusions (additional wind drag inherited from base class if needed).
    Suitable for sediment tracers, e.g. for tracking sediment particles.
    Adapted from OpenDrift3DSimulation by Simon Weppe - MetOcean Solutions.

    """
    ElementType = BuoyantTracer # simply use BuoyantTracer for now - will eventually move to SedimentTracer
    required_variables = [
        'x_sea_water_velocity',
        'y_sea_water_velocity',
        'x_wind', 'y_wind',
        'upward_sea_water_velocity',
        'ocean_vertical_diffusivity',
        'ocean_horizontal_diffusivity',
        'sea_surface_wave_significant_height',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'sea_floor_depth_below_sea_level'
        ]

    required_variables.append('land_binary_mask')

    required_profiles = ['ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    fallback_values = {
        'x_sea_water_velocity': 0,
        'y_sea_water_velocity': 0,
        'sea_surface_wave_significant_height': 0,
        'sea_surface_wave_stokes_drift_x_velocity': 0,
        'sea_surface_wave_stokes_drift_y_velocity': 0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
        'x_wind': 0,
        'y_wind': 0,
        'upward_sea_water_velocity': 0.0,
        'ocean_vertical_diffusivity': 0.0, # in m2/s
        'ocean_horizontal_diffusivity' : 0.0, # in m2/s2
        'sea_floor_depth_below_sea_level': 10000
        }
    
    # Adding some specs - inspired from basereader.py
    #
    # Default plotting colors of trajectory endpoints

    status_colors_default = {'initial': 'green',
                             'active': 'blue',
                             'missing_data': 'gray',
                             'settled': 'red'}
    # adding processes switch to contol resuspension
    configspec_sedimentdrift3d = '''
        [processes]
            resuspension = boolean(default=False)
    '''

    def __init__(self, *args, **kwargs):
        # update configstring
        self._add_configstring(self.configspec_sedimentdrift3d)

        # Calling general constructor of parent class
        super(SedimentDrift3D, self).__init__(*args, **kwargs)

        # Turbulent mixing switched off by default
        self.set_config('processes:turbulentmixing', False)
        # resuspension switched off by default
        self.set_config('processes:resuspension', False)

        self.max_speed = 1.5

    def update_terminal_velocity(self, *args, **kwargs):
        '''
        Terminal velocity due to buoyancy or sedimentation rate,
        to be used in turbulent mixing module.
        stored as an array in self.elements.terminal_velocity
        '''
        # conserve user-input terminal_velocity
        #
        # could include some  more sophisticated calculations of terminal_velocity/settling velocity based on temp/salinity etc..
        #
        # self.elements.terminal_velocity = 0.
        pass

    def horizontal_diffusion(self, *args, **kwargs):
        '''
        Horizontal diffusion based on random walk technique
        using diffusion coefficients K_xy (in m2/s - constant or interpolaote from field)
        and uniformly distributed random numbers R(-1,1) with a zero mean.

        
        Constant coefficients 'ocean_horizontal_diffusivity' can be set using:
        o.fallback_values['ocean_horizontal_diffusivity'] = 0.1 

        Time-Space varying coefficients  'ocean_horizontal_diffusivity' 
        can also be interpolated from an environment object (as in vertical_mixing() )
        
        * this should not be used in combination with drift:current_uncertainty ~=0 
        * as this model the same process, only with a different approach 

        The random component added to the particle position to reproduce turbulent horizontal 
        is computed using a approach similar to PyGnome, CMS :

        -https://github.com/beatrixparis/connectivity-modeling-system/blob/master/cms-master/src/mod_turb.f90
        -Garcia-Martinez and Tovar, 1999 - Computer Modeling of Oil Spill Trajectories With a High Accuracy Method
        -Lonin, S.A., 1999. Lagrangian model for oil spill diffusion at sea. Spill Science and Technology Bulletin, 5(5): 331-336
        -https://github.com/NOAA-ORR-ERD/PyGnome/blob/master/lib_gnome/Random_c.cpp#L50 

        stored as an array in self.elements.horizontal_diffusion
 
        '''
        if not self.environment.ocean_horizontal_diffusivity.any():
            logging.debug('No horizontal diffusion applied - ocean_horizontal_diffusivity = 0.0')
            pass

        # check if some diffusion is not already accounted for using drift:current_uncertainty
        if self.get_config('drift:current_uncertainty') != 0:
            logging.debug('Warning = some horizontal diffusion already accounted for using drift:current_uncertainty')

        diff_fac = 6 # 6 is used in PyGnome as well, but can be = 2 in other formulations, such as in the CMS model - hard coded for now
        K_xy = self.environment.ocean_horizontal_diffusivity # horizontal diffusion coefficient in [m2/s]
        # max diffusion distance in meters is : max_diff_distance = np.sqrt(diff_fac * K_xy * self.time_step.seconds)
        # here we need velocity to fit with the update_position() subroutine
        # max_diff_velocity = max_diff_distance / self.time_step.seconds = np.sqrt(diff_fac * K_xy / self.time_step.seconds) 
        max_diff_velocity = np.sqrt(diff_fac * K_xy / self.time_step.seconds) # max diffusion velocity in [m/s]
        # add randomness - random number in [-1,1] using uniform distribution i.e. same probabilty for each value (note some other formulations may use "normal" distribution)
        diff_velocity = np.random.uniform(-1,1,size = len(max_diff_velocity)) * max_diff_velocity # in [m/s]
        # random directions
        theta_rand = np.random.uniform( 0, 2*np.pi, size = len(max_diff_velocity) )
        # split diff_velocity into (u,v) velocities using random directions
        x_vel = diff_velocity * np.cos(theta_rand)
        y_vel = diff_velocity * np.sin(theta_rand)

        logging.debug('Applying horizontal diffusion to particle positions')
        logging.debug('\t\t%s   <- horizontal diffusion distance [m] ->   %s' % (np.min(diff_velocity*self.time_step.seconds), np.max(diff_velocity*self.time_step.seconds)))

        # update positions with the diffusion velocities     
        self.update_positions(x_vel, y_vel)

    def update(self):
        """Update positions and properties of elements."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Simply move particles with ambient current
        self.advect_ocean_current() # from physics_methods.py

        # Horizontal diffusion (using diffusion coefficient approach )
        self.horizontal_diffusion()

        # Advect particles due to wind drag
        # (according to specified wind_drift_factor)
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Turbulent Mixing
        self.update_terminal_velocity() # routine to estimate settling velocity - simply keeps the user-input one for now 
        self.vertical_mixing() # using subroutine from opendrift3D.py - include buoyancy-related settling and mixing

        # Vertical advection
        self.vertical_advection()

        # Sediment resuspension checks , if switched on
        # self.sediment_resuspension() - ToDo!

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')
        # When no resuspension is required, deactivate that reached the seabed
        if self.get_config('processes:resuspension') is False:
            self.deactivate_elements(self.elements.z ==
                                     -1.*self.environment.sea_floor_depth_below_sea_level,
                                     reason='settled')

        # Note the interaction with shoreline in taken care of by interact_with_coastline in basemodel.py
        # when run() is called

    def sediment_resuspension(self):
        """
        Compute ambient bed shear stresses at positions of particles that are (or just touched) the bottom
        and determine if they can settle or be resuspended based on critical_shear_stress 
        """
        # 1-find particles on the bottom
        # 2-compute bed shear stresses
        # 3-compare to critical_shear_stress
        # 4-resuspend or stay on seabed depending on 3)
        #   > probably need to use a cut-off age after which particles are de-activated anyway
        #   to prevent excessive build-up of "active" particle in the simulations


# General physics functions---------------------------------------------------------------------------------
# Could be moved to physics_methods.py once cross-checked / accepted 

def bedshearstress_cw(self,p,time=None,imax=2):
    """Computation of bed shear stress due to current and waves
    current-related stress is computed following a drag-coefficient approach
    wave-related stress is computed following Van Rijn approach
    combined wave-current mean and max stresses are computed followin Soulsby(1995) approach
    Arguments:
      self : Material object, expected to include fields movers, reactors (if input)
      p: particle positions array (Nx3)
      time: 
    Returns:
      tau_cur : current-related bed shear stress tau_cur
      tau_cw : combined mean current-wave bed shear stress 
      tau_max: combined max current-wave bed shear stress tau_max
      topo : water depth at particle positions (of first mover) 
    """
    rhow=1027 # default volumic mass for seawater
    tau_cur=numpy.tile(0.0,numpy.size(p,0)) # allocate
    tau_cw=numpy.tile(0.0,numpy.size(p,0)) # allocate
    tau_max=numpy.tile(0.0,numpy.size(p,0)) # allocate
    # current-related bed shear stress (sum of all movers)
    for mover in self.movers[0:]:
      if mover.topo: # topo needed to define bedshear stress
        topo=mover.topo.interp(p,None,3)     
        if (not mover.is3d) and (mover.z0>0): # mover is a 2D-depth averaged current
          # temporarily set mover.z0 to 0.0 so that mover.interp yields the un-corrected depth-averaged current (direct interpolation, no log profile )
          # not super elegant, probably a better way to do this - anyway to access GriddedTide or GriddedData from here ?
          z0_tmp=copy.copy(mover.z0)
          mover.z0=0.0
          u2dhim=mover.interp(p,time,imax)   #u2dhim=mover.interp(self,p,time,imax)        
          mover.z0=z0_tmp
          u2dhim_mag=(u2dhim[:,0]**2+u2dhim[:,1]**2)**0.5
          # Drag coefficient for 2D case using water depth and z0 (see COHERENS manual eq.7.2, or Delft3d)
          Cdrag=( 0.4 /(numpy.log(abs(topo[:,0] /mover.z0))-1) )**2
          #Now compute the bed shear stress [N/m2] 
          tau_cur+=rhow*Cdrag*u2dhim_mag**2             
        elif (mover.is3d) and (mover.z0>0):   # mover is a 3D current field
          #import pdb;pdb.set_trace()
          # Assume the first grid point above the bed is assumed to be the top of the logarithmic boundary layer
          # the log profile extends from that last wet bin level, to the bottom
          # see COHERENS manual eq 7.1/7.2
  
          # find closest "wet" vertical levels at each particle locations
          bin_lev=numpy.zeros(len(p[:,0]))
          for lev in mover.lev:
            bin_lev[topo[:,0]<=lev]=lev
          #vertical height from last wet vertical bin to seabed
          zb=bin_lev-topo[:,0]
          #current computed at last wet vertical bin
          uub=mover.interp(numpy.vstack((p[:,0],p[:,1],bin_lev)).T,time,imax)
          uub_mag=(uub[:,0]**2+uub[:,1]**2)**0.5
          # Drag coefficient for 3D case using zb and z0 (see COHERENS manual eq.7.2, or Delft3d)
          Cdrag=( 0.4 /(numpy.log(abs(zb /mover.z0))-1) )**2 
         #Now compute the bed shear stress [N/m2]
          tau_cur+=rhow*Cdrag*uub_mag**2

    # wave-related bed shear stress
    if len(self.reactors)>0: # check if wave forcing is included 
      # for now assume that if reactors exist, they will be correctly input
      # computation of wave-related and combined bed shear stresses based on code from 
      #https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/phys_fun/bedshearstresses.m 
      #https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/phys_fun/sandandmudtransport.m
      hs=self.reactors['hs'].interp(p[:1,:],time)[:,0] #wave height
      tp=self.reactors['tp'].interp(p[:1,:],time)[:,0] #wave period
      # wave-related roughness
      # vanRijn 2007 suggests same equations than for current-related roughness where 20* d50 <ksw<150*d50
      # here we are using nikuradse for consistency with the use of z0 in the mover class for now
      ksw=30*self.movers[0].z0  
      topo= self.movers[0].topo.interp(p,None,3)
      w=2*numpy.pi/tp
      kh = qkhfs( w, topo[:,0] ) # dispersion relationship
      Adelta = hs/(2*numpy.sinh(kh)) # peak wave orbital excursion
      Udelta = (numpy.pi*hs)/(tp*numpy.sinh(kh))  # peak wave orbital velocity
      fw = numpy.exp(-5.977+5.213*(Adelta/ksw)**-0.194)  # wave-related friction coefficient (van Rijn)
      fw = numpy.min(fw,0.3)
      tau_wave = 0.25 * rhow * fw * (Udelta)**2 # wave-related bed shear stress
      #cycle mean bed shear stress according to Soulsby,1995
      tau_cw=tau_cur*[1+1.2*(tau_wave/(tau_cur+tau_wave))**3.2]
      # max bed shear stress during wave cycle
      tau_max=[tau_wave**2+tau_cur**2]**0.5
    else:
      tau_max=tau_cur
      tau_cw=tau_cur
    # if (tau_cur==0).any():
    #   import pdb;pdb.set_trace()
    return tau_cur,tau_cw,tau_max,topo

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
    y = numpy.sqrt(x) * (x<1.) + x *(x>=1.)
    # is this faster than a loop?
    t = numpy.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = numpy.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = numpy.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    kh = y
    return kh