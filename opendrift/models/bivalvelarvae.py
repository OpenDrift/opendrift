# This file has been adapted from pelagicegg.py but introduces a settlement competency period
# added config setting for 'min_settlement_age_seconds' which defines the minimum age after which  settlement can occur.
# Coastal settlement occurs after a competent particle interacts with the coastline i.e. after age>min_settlement_age_seconds
# Benthic settlement occurs when a competent particle makes contact with the bottom i.e. after age>min_settlement_age_seconds
# Habitat type is not considered in either coastal or benthic settlement
# 
# The model overload the seabed and coastline interaction subfunction from basemodel.py :
#    > lift_elements_to_seafloor
#    > interact_with_coastline
# 
#  Authors : Craig Norrie, Simon Weppe
# 
# 
#  Under development - more testing to do
# 
# 


import numpy as np
from opendrift.models.opendrift3D import \
    OpenDrift3DSimulation, Lagrangian3DArray
from opendrift.elements import LagrangianArray
import logging

# Defining the  element properties from Pelagicegg model
class BivalveLarvaeObj(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs/larvae
    """

    variables = LagrangianArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': '[]',
                                       'default': 31.25}),  # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 1028.}),
        ('hatched', {'dtype': np.float32,
                     'units': '',
                     'default': 0.}),
        ('terminal_velocity', {'dtype': np.float32,
                       'units': 'm/s',
                       'default': 0.})])


class BivalveLarvae(OpenDrift3DSimulation):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Under construction.
    """

    ElementType = BivalveLarvaeObj
    # ElementType = BuoyantTracer 

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature',
                          'sea_water_salinity',
                          'sea_surface_height',
                          'surface_downward_x_stress',
                          'surface_downward_y_stress',
                          'turbulent_kinetic_energy',
                          'turbulent_generic_length_scale',
                          'upward_sea_water_velocity'
                          ]

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented)

    # required_profiles = ['sea_water_temperature',
    #                      'sea_water_salinity',
    #                      'ocean_vertical_diffusivity']

    # removing salt/water temp profile requirement for now
    # > need to get correct profiles from SCHISM reader

    required_profiles = ['ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'ocean_vertical_diffusivity': 0.02,  # m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity': 34.,
                       'sea_surface_height': 0.0,
                       'surface_downward_x_stress': 0,
                       'surface_downward_y_stress': 0,
                       'turbulent_kinetic_energy': 0,
                       'turbulent_generic_length_scale': 0,
                       'upward_sea_water_velocity': 0
                       }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'settled_on_coast': 'red', 'died': 'yellow', 'settled_on_bottom': 'magenta'}

    configspec_bivalvelarvae = ''' 
        [drift]
            min_settlement_age_seconds = float(default=0)
            '''

    def __init__(self, *args, **kwargs):
        
        ##add config spec
        self._add_configstring(self.configspec_bivalvelarvae)

        # Calling general constructor of parent class
        super(BivalveLarvae, self).__init__(*args, **kwargs)

        # By default, larvae do not strand when reaching shoreline. 
        # They are recirculated back to previous position instead
        self.set_config('general:coastline_action', 'previous')

        # resuspend larvae that reach seabed by default 
        self.set_config('drift:lift_to_seafloor', 'True')
        # set the deafult min_settlement_age_seconds to 0.0
        # self.set_config('drfit:min_settlement_age_seconds', '0.0')

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        pass
#       self.elements.terminal_velocity = W

    def sea_surface_height(self):
        '''Sea surface height for presently active elements

           sea_surface_height > 0 above mean sea level
           sea_surface_height < 0 below mean sea level
        '''
        if hasattr(self, 'environment') and \
                hasattr(self.environment, 'sea_surface_height'):
            if len(self.environment.sea_surface_height) == \
                    self.num_elements_active():
                sea_surface_height = \
                    self.environment.sea_surface_height
        if 'sea_surface_height' not in locals():
            env, env_profiles, missing = \
                self.get_environment(['sea_surface_height'],
                                     time=self.time, lon=self.elements.lon,
                                     lat=self.elements.lat,
                                     z=0*self.elements.lon, profiles=None)
            sea_surface_height = \
                env['sea_surface_height'].astype('float32') 
        return sea_surface_height   

    def update(self):
        """Update positions and properties of buoyant particles."""

        # Turbulent Mixing
        self.update_terminal_velocity()
        self.vertical_mixing()

        # Horizontal advection
        self.advect_ocean_current()

        # Vertical advection
        if self.get_config('processes:verticaladvection') is True:
            self.vertical_advection()           

    def lift_elements_to_seafloor(self):  ###Initiate settlement if particles touch bottom during competence period
        '''Lift any elements which are below seafloor and check age
          (overloads the lift_elements_to_seafloor() from basemodel.py)

           The methods will check age of larvae that touched the seabed.
             if age_particle < min_settlement_age_seconds : resuspend larvae
             if age_particle > min_settlement_age_seconds : larvaes settle and will be deactivated.

        '''
            
        if 'sea_floor_depth_below_sea_level' not in self.priority_list:
            return
        
        sea_floor_depth = self.sea_floor_depth() # returns a positive down water depth
        sea_surface_height = self.sea_surface_height() # returns surface elevation at particle positions (>0 above msl, <0 below msl)

        below = self.elements.z < -sea_floor_depth
        
        if self.get_config('drift:lift_to_seafloor') is True:
            # self.elements.z[below] = -sea_floor_depth[below] - intial code

            self.elements.z[below] = np.minimum(-sea_floor_depth[below], sea_surface_height[below])
            # make sure particles dont get above water at this stage i.e. z>sea_surface_height
            # this can happen when reader has negative values
            # e.g. : sea_floor_depth() returns e.g. -2.0 (i.e. wetting-drying points)
        else:
            self.deactivate_elements(below, reason='seafloor')

        # Deactivate elements that touched seabed and have age>min_settlement_age_seconds
        if self.get_config('drift:min_settlement_age_seconds') != 0.0:
            older = self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds')
            if older.any():
                self.deactivate_elements(older & below ,reason='settled_on_bottom')

    def surface_stick(self):
        '''Keep particles just below the surface.
           (overloads the OpenDrift3DSimulation version to allow for possibly time-varying
           sea_surface_height)
        '''
        
        sea_surface_height = self.sea_surface_height() # returns surface elevation at particle positions (>0 above msl, <0 below msl)
        
        # keep particle just below sea_surface_height (self.elements.z depth are negative down)
        surface = np.where(self.elements.z >= sea_surface_height)
        if len(surface[0]) > 0:
            self.elements.z[surface] = sea_surface_height -0.01 # set particle z at 0.01m below sea_surface_height

    def interact_with_coastline(self,final = False): 
        """Coastline interaction according to configuration setting
           (overloads the interact_with_coastline() from basemodel.py)
           
           The method checks for age of particles that intersected coastlines:

             if age_particle < min_settlement_age_seconds : move larvaes back to previous wet position
             if age_particle > min_settlement_age_seconds : larvaes become stranded and will be deactivated.

        """
        i = self.get_config('general:coastline_action')

        if not hasattr(self.environment, 'land_binary_mask'):
            return

        if i == 'none':  # Do nothing
            return

        if final is True:  # Get land_binary_mask for final location
            en, en_prof, missing = \
                self.get_environment(['land_binary_mask'],
                                     self.time,
                                     self.elements.lon,
                                     self.elements.lat,
                                     self.elements.z,
                                     None)
            self.environment.land_binary_mask = en.land_binary_mask

        if i == 'stranding':  # Deactivate elements on land
            self.deactivate_elements(
                self.environment.land_binary_mask == 1, reason='stranded')

        elif i == 'previous':  # Go back to previous position (in water)
            # previous_position_if = self.previous_position_if()
            if self.newly_seeded_IDs is not None:
                self.deactivate_elements(
                    (self.environment.land_binary_mask == 1) &
                    (self.elements.ID >= self.newly_seeded_IDs[0]),
                    reason='seeded_on_land')
            on_land = np.where(self.environment.land_binary_mask == 1)[0]

                # if previous_position_if is not None:
                #     self.deactivate_elements((previous_position_if*1 == 1) & (
                #                      self.environment.land_binary_mask == 0),
                #                          reason='seeded_at_nodata_position')

            # if previous_position_if is None:
            #     on_land = np.where(self.environment.land_binary_mask == 1)[0]
            # else:
            #     on_land = np.where((self.environment.land_binary_mask == 1) |
            #                        (previous_position_if == 1))[0]
            if len(on_land) == 0:
                logging.debug('No elements hit coastline.')
            else:                
                if self.get_config('drift:min_settlement_age_seconds') == 0.0 :
                    # No minimum age input, set back to previous position
                    logging.debug('%s elements hit coastline, '
                              'moving back to water' % len(on_land))
                    on_land_ID = self.elements.ID[on_land]
                    self.elements.lon[on_land] = \
                        np.copy(self.previous_lon[on_land_ID - 1])
                    self.elements.lat[on_land] = \
                        np.copy(self.previous_lat[on_land_ID - 1])
                    self.environment.land_binary_mask[on_land] = 0
                else:
                    # Minimum age before settling was input; check age of particle versus min_settlement_age_seconds
                    # and strand or recirculate accordingly
                    on_land_and_younger = (self.environment.land_binary_mask == 1) & (self.elements.age_seconds < self.get_config('drift:min_settlement_age_seconds'))
                    on_land_and_older = (self.environment.land_binary_mask == 1) & (self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds'))

                    logging.debug('%s elements hit coastline' % len(on_land))
                    logging.debug('moving %s elements younger than min_settlement_age_seconds back to previous water position' % len(np.where(on_land_and_younger)[0]) )
                    logging.debug('%s elements older than min_settlement_age_seconds remain stranded on coast' % len(np.where(on_land_and_older)[0]) )
                    
                    # refloat elements younger than min_settlement_age back to previous position(s)
                    if on_land_and_younger.any() :
                        self.elements.lon[np.where(on_land_and_younger)] = np.copy(self.previous_lon[np.where(on_land_and_younger)])  
                        self.elements.lat[np.where(on_land_and_younger)] = np.copy(self.previous_lat[np.where(on_land_and_younger)])
                        self.environment.land_binary_mask[on_land_and_younger] = 0 
                    # deactivate elements older than min_settlement_age & save position
                    # ** function expects an array of size consistent with self.elements.lon
                    self.deactivate_elements(on_land_and_older, reason='settled_on_coast')


    def increase_age_and_retire(self):  ####So that if max_age_seconds is exceeded particle is flagged as died
            """Increase age of elements, and retire if older than config setting."""
            # Increase age of elements
            self.elements.age_seconds += self.time_step.total_seconds()

            # Deactivate elements that exceed a certain age
            if self.get_config('drift:max_age_seconds') is not None:
                self.deactivate_elements(self.elements.age_seconds >=
                                         self.get_config('drift:max_age_seconds'),
                                         reason='died')

            # Deacticate any elements outside validity domain set by user
            if self.validity_domain is not None:
                W, E, S, N = self.validity_domain
                if W is not None:
                    self.deactivate_elements(self.elements.lon < W, reason='outside')
                if E is not None:
                    self.deactivate_elements(self.elements.lon > E, reason='outside')
                if S is not None:
                    self.deactivate_elements(self.elements.lat < S, reason='outside')
                if N is not None:
                    self.deactivate_elements(self.elements.lat > N, reason='outside')