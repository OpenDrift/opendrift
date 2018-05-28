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
        # do nothing 
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
        # Horizontal diffusion
        self.horizontal_diffusion()

        # Advect particles due to wind drag
        # (according to specified wind_drift_factor)
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Turbulent Mixing
        self.update_terminal_velocity()
        self.vertical_mixing() # using subroutine from opendrift3D.py

        # Vertical advection
        self.vertical_advection()

        # Sediment resuspension checks , if switched on
        # self.sediment_resuspension() - ToDo!
        # 
        # 1-find particles on the bottom
        # 2-compute bed shear stresses
        # 3-compare to critical_shear_stress
        # 4-resuspend or stay on seabed depending on 3)
        #   > probably need to use a cut-off age after which particles are de-activated anyway
        #   to prevent excessive build-up of "active" particle in the simulations

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
