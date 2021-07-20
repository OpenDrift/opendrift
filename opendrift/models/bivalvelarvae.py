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
# Modified: 28/06/2021. RChaput: Add settlement restricted to suitable habitat only, vertical swimming behavior, and maximum depth of dispersal.
#			         Add Haliotis iris specific behavior: o.set_config('drift:Haliotis_iris', True) => sink for 12 hours, then swim toward the surface using the user specified vertical terminal velocity
#
# Lines to add to script to run the vertical swimming code:
# ################################
# # Type of settlement
# ###############################
#o.habitat('./habitat/Name_of_your_shapefile.shp') # Location of the shapefile with the habitat
#o.set_config('drift:settlement_in_habitat', True) # Settlement restricted to habitat only, default=False
#o.set_config('drift:max_age_seconds', 10*24*3600) # Maximum PLD
#o.set_config('drift:min_settlement_age_seconds', 5*24*3600) # Beginning of the Competency period
#o.set_config('general:seafloor_action', 'lift_to_seafloor')
#
# ###############################
# # Vertical swimming
# ###############################
# # Need a null terminal_velocity
#o.set_config('drift:active_vertical_swimming', True) # Correlated random walk across the water column when advected away from coastal habitats
#o.set_config('drift:vertical_velocity', 0.0025) # Vertical swimming speed of the larvae: in meter/seconds
#o.set_config('drift:maximum_depth', -50.0) # Maximum depth of larvae: negative in meters
#o.set_config('drift:persistence', 50) # Control the persistence (memory) of the vertical movement to create a correlated random walk and sample the water column (0= pure random walk)


import numpy as np
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
import logging; logger = logging.getLogger(__name__)
from shapely.geometry import Polygon, Point, MultiPolygon # added for settlement in polygon only
import fiona # to import habitat
import random
from sklearn.neighbors import BallTree


# Defining the  element properties from Pelagicegg model
class BivalveLarvaeObj(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs/larvae
    """

    variables = Lagrangian3DArray.add_variables([
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
        ('vertical_movement', {'dtype': np.float32,
                     'units': 'none',
                     'default': 1.0}),
        ('terminal_velocity', {'dtype': np.float32,
                       'units': 'm/s',
                       'default': 0.})])


class BivalveLarvae(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Under construction.
    """

    ElementType = BivalveLarvaeObj
    # ElementType = BuoyantTracer 

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_ice_area_fraction': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
        'sea_water_temperature': {'fallback': 15, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'sea_surface_height': {'fallback': 0.0},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
      }

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

    # required_profiles = ['ocean_vertical_diffusivity']

    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'settled_on_coast': 'red', 'died': 'yellow', 'settled_on_bottom': 'magenta'}

    def __init__(self, *args, **kwargs):
        
        # Calling general constructor of parent class
        super(BivalveLarvae, self).__init__(*args, **kwargs)

        # By default, larvae do not strand when reaching shoreline. 
        # They are recirculated back to previous position instead
        self.set_config('general:coastline_action', 'previous')

        # resuspend larvae that reach seabed by default 
        self.set_config('general:seafloor_action','lift_to_seafloor')

        ##add config spec
        self._add_config({ 'drift:min_settlement_age_seconds': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'minimum age in seconds at which larvae can start to settle on seabed or stick to shoreline)',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'drift:settlement_in_habitat': {'type': 'bool', 'default': False,
                           'description': 'settlement restricted to suitable habitat only',
                           'level': self.CONFIG_LEVEL_BASIC}}) 
        self._add_config({ 'drift:active_vertical_swimming': {'type': 'bool', 'default': False,
                           'description': 'vertical movements correlated with the direction of inshore currents',
                           'level': self.CONFIG_LEVEL_BASIC}}) 
        self._add_config({ 'drift:persistence': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'none',
                           'description': 'persistence of the vertical swimming when sampling the currents',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'drift:vertical_velocity': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1, 'units': 'm/s',
                           'description': 'vertical swimming of larvae when advected away from habitat',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'drift:maximum_depth': {'type': 'float', 'default': None,'min': -10000.0, 'max': -1.0, 'units': 'm',
                           'description': 'maximum depth of the larvae',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'drift:Haliotis_iris': {'type': 'bool', 'default': False,
                           'description': 'turns on a sinking behavior for the first 12 hours of dispersal',
                           'level': self.CONFIG_LEVEL_BASIC}})

    
    def habitat(self, shapefile_location):
        """Suitable habitat in a shapefile"""
        polyShp = fiona.open(shapefile_location) # import shapefile
        polyList = []
        centers = []
        rad_centers = []
        for poly in polyShp: # create individual polygons from shapefile
             polyGeom = Polygon(poly['geometry']['coordinates'][0])
             polyList.append(polyGeom) # Compile polygon in a list 
             centers.append(polyGeom.centroid.coords[0]) # Compute centroid and return a [lon, lat] list
        for poly in range(len(centers)):
            rad_centers.append([np.deg2rad(centers[poly][1]),np.deg2rad(centers[poly][0])])
        self.multiShp = MultiPolygon(polyList).buffer(0) # Aggregate polygons in a MultiPolygon object and buffer to fuse polygons and remove errors
        self.ball_centers = BallTree(rad_centers, metric='haversine') # Create a Ball Tree with the centroids for faster computation
        return self.multiShp, self.ball_centers

    
    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        pass
#       self.elements.terminal_velocity = W

    

    def sea_surface_height(self):
        '''fetches sea surface height for presently active elements

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


#    def lift_elements_to_seafloor(self):  ###Initiate settlement if particles touch bottom during competence period
#        '''Lift any elements which are below seafloor and check age
#          (overloads the lift_elements_to_seafloor() from basemodel.py)
#
#           The methods will check age of larvae that touched the seabed.
#             if age_particle < min_settlement_age_seconds : resuspend larvae
#             if age_particle > min_settlement_age_seconds : larvaes settle and will be deactivated.
#
#        '''
#            
#        if 'sea_floor_depth_below_sea_level' not in self.priority_list:
#            return
#        
#        sea_floor_depth = self.sea_floor_depth() # returns a positive down water depth
#        sea_surface_height = self.sea_surface_height() # returns surface elevation at particle positions (>0 above msl, <0 below msl)
#
#        below = self.elements.z < -sea_floor_depth
#        
#        if self.get_config('drift:lift_to_seafloor') is True:
#            # self.elements.z[below] = -sea_floor_depth[below] - intial code
#
#            self.elements.z[below] = np.minimum(-sea_floor_depth[below], sea_surface_height[below])
#            # make sure particles dont get above water at this stage i.e. z>sea_surface_height
#            # this can happen when reader has negative values for sea_floor_depth
#            # e.g. : sea_floor_depth() returns e.g. -2.0 (i.e. wetting-drying points)
#            # 
#            # if sea_surface_height is not available from reader, fallback value is used (=0 by default).
#        else:
#            self.deactivate_elements(below, reason='seafloor')
#
#        # Deactivate elements that touched seabed and have age>min_settlement_age_seconds
#        if self.get_config('drift:min_settlement_age_seconds') != 0.0:
#            older = self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds')
#            if older.any():
#                self.deactivate_elements(older & below ,reason='settled_on_bottom')

                
   def interact_with_seafloor(self):
        """Seafloor interaction according to configuration setting"""
        # 
        # This function will overloads the version in basemodel.py
        if self.num_elements_active() == 0:
            return
        if 'sea_floor_depth_below_sea_level' not in self.priority_list:
            return
        sea_floor_depth = self.sea_floor_depth()
        below = np.where(self.elements.z < -sea_floor_depth)[0]
        if len(below) == 0:
                logger.debug('No elements hit seafloor.')
                return

        below_and_older = np.logical_and(self.elements.z < -sea_floor_depth, 
            self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds'))
        below_and_younger = np.logical_and(self.elements.z < -sea_floor_depth, 
            self.elements.age_seconds < self.get_config('drift:min_settlement_age_seconds'))
        
        # Move all elements younger back to seafloor 
        # (could rather be moved back to previous if relevant? )
        self.elements.z[np.where(below_and_younger)] = -sea_floor_depth[np.where(below_and_younger)]

        # deactivate elements that were both below and older
        if self.get_config('drift:settlement_in_habitat') is False:
            self.deactivate_elements(below_and_older ,reason='settled_on_bottom')
        # if elements can only settle in habitat then they are moved back to seafloor
        else:
            self.elements.z[np.where(below_and_older)] = -sea_floor_depth[np.where(below_and_older)]

        logger.debug('%s elements hit seafloor, %s were older than %s sec. and deactivated, %s were lifted back to seafloor' \
            % (len(below),len(below_and_older),self.get_config('drift:min_settlement_age_seconds'),len(below_and_younger)))             
                
                
    def surface_stick(self):
        '''Keep particles just below the surface.
           (overloads the OpenDrift3DSimulation version to allow for possibly time-varying
           sea_surface_height)
        '''
        
        sea_surface_height = self.sea_surface_height() # returns surface elevation at particle positions (>0 above msl, <0 below msl)
        
        # keep particle just below sea_surface_height (self.elements.z depth are negative down)
        surface = np.where(self.elements.z >= sea_surface_height)
        if len(surface[0]) > 0:
            self.elements.z[surface] = sea_surface_height[surface] -0.01 # set particle z at 0.01m below sea_surface_height

            
    def interact_with_coastline(self,final = False): 
        """Coastline interaction according to configuration setting
           (overloads the interact_with_coastline() from basemodel.py)
           
           The method checks for age of particles that intersected coastlines:
             if age_particle < min_settlement_age_seconds : move larvaes back to previous wet position
             if age_particle > min_settlement_age_seconds : larvaes become stranded and will be deactivated.
        """
        i = self.get_config('general:coastline_action') # will always be 'previous'

        if not hasattr(self.environment, 'land_binary_mask'):
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

        # if i == 'previous':  # Go back to previous position (in water)
        # previous_position_if = self.previous_position_if()
        if self.newly_seeded_IDs is not None:
                self.deactivate_elements(
                    (self.environment.land_binary_mask == 1) &
                    (self.elements.age_seconds == self.time_step.total_seconds()),
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
            logger.debug('No elements hit coastline.')
        else:
            if self.get_config('drift:settlement_in_habitat') is True:
                    # Particle can only settle in habitat, set back to previous location
                    logger.debug('%s elements hit coastline, '
                              'moving back to water' % len(on_land))
                    on_land_ID = self.elements.ID[on_land]
                    self.elements.lon[on_land] = \
                        np.copy(self.previous_lon[on_land_ID - 1])
                    self.elements.lat[on_land] = \
                        np.copy(self.previous_lat[on_land_ID - 1])
                    self.environment.land_binary_mask[on_land] = 0  
            elif self.get_config('drift:min_settlement_age_seconds') == 0.0 :
                # No minimum age input, set back to previous position (same as in interact_with_coastline() from basemodel.py)
                logger.debug('%s elements hit coastline, '
                          'moving back to water' % len(on_land))
                on_land_ID = self.elements.ID[on_land]
                self.elements.lon[on_land] = \
                    np.copy(self.previous_lon[on_land_ID - 1])
                self.elements.lat[on_land] = \
                    np.copy(self.previous_lat[on_land_ID - 1])
                self.environment.land_binary_mask[on_land] = 0
            else:
                #################################
                # Minimum age before settling was input; check age of particle versus min_settlement_age_seconds
                # and strand or recirculate accordingly
                on_land_and_younger = np.where((self.environment.land_binary_mask == 1) & (self.elements.age_seconds < self.get_config('drift:min_settlement_age_seconds')))[0]
                on_land_and_older = np.where((self.environment.land_binary_mask == 1) & (self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds')))[0]

                # this step replicates what is done is original code, but accounting for particle age. It seems necessary 
                # to have an array of ID, rather than directly indexing using the "np.where-type" index (in dint64)
                on_land_and_younger_ID = self.elements.ID[on_land_and_younger] 
                on_land_and_older_ID = self.elements.ID[on_land_and_older]

                logger.debug('%s elements hit coastline' % len(on_land))
                logger.debug('moving %s elements younger than min_settlement_age_seconds back to previous water position' % len(on_land_and_younger))
                logger.debug('%s elements older than min_settlement_age_seconds remain stranded on coast' % len(on_land_and_younger))
                
                # refloat elements younger than min_settlement_age back to previous position(s)
                if len(on_land_and_younger) > 0 :
                    # self.elements.lon[np.where(on_land_and_younger)] = np.copy(self.previous_lon[np.where(on_land_and_younger)])  
                    # self.elements.lat[np.where(on_land_and_younger)] = np.copy(self.previous_lat[np.where(on_land_and_younger)])
                    # self.environment.land_binary_mask[on_land_and_younger] = 0 

                    self.elements.lon[on_land_and_younger] = np.copy(self.previous_lon[on_land_and_younger_ID - 1])
                    self.elements.lat[on_land_and_younger] = np.copy(self.previous_lat[on_land_and_younger_ID - 1])
                    self.environment.land_binary_mask[on_land_and_younger] = 0

                # deactivate elements older than min_settlement_age & save position
                # ** function expects an array of size consistent with self.elements.lon
                self.deactivate_elements((self.environment.land_binary_mask == 1) & \
                                         (self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds')),
                                         reason='settled_on_coast')
     
    
    def interact_with_habitat(self):
        """Habitat interaction according to configuration setting
        The method checks if a particle is within the limit of an habitat before to allow settlement
        """        
        # Get age of particle
        old_enough = np.where(self.elements.age_seconds >= self.get_config('drift:min_settlement_age_seconds'))[0]
        # Extract particles positions
        if len(old_enough) > 0 :
            pts_lon = self.elements.lon[old_enough]
            pts_lat = self.elements.lat[old_enough]
            for i in range(len(pts_lon)): # => faster version
                 pt = Point(pts_lon[i], pts_lat[i])
                 in_habitat = pt.within(self.multiShp)
                 if in_habitat == True:
                     self.environment.land_binary_mask[old_enough[i]] = 6
                     
        # Deactivate elements that are within a polygon and old enough to settle
        # ** function expects an array of size consistent with self.elements.lon                
        self.deactivate_elements((self.environment.land_binary_mask == 6), reason='settled_in_habitat')
    
    
    def vertical_swimming(self):
        """"Vertical movements of the larvae when caught in inshore currents: larvae try to stay within the currents that move
        them closer to shore, otherwise the larvae go up and down to sample the water column
        """
        # Check distance with nearest_habitat using the Ball tree algorithm
        dist_old = self.ball_centers.query(list(zip(np.deg2rad(self.previous_lat), np.deg2rad(self.previous_lon))), k=1)
        dist = self.ball_centers.query(list(zip(np.deg2rad(self.elements.lat), np.deg2rad(self.elements.lon))), k=1)
        # Prompt vertical movement for particles moving away from the habitat
        for i in range(len(self.elements.lat)):
            if dist[0][i] > dist_old[0][i]:
                rand = random.uniform(0,1)
                movementr = self.elements.vertical_movement[i] if rand > 1/(2+self.get_config('drift:persistence')) else -1 if random.uniform(0,1) < 0.5 else 1 # Correlated random walk
                self.elements.vertical_movement[i] = movementr      
                self.elements.z[i] = self.elements.z[i] + movementr*self.get_config('drift:vertical_velocity') * self.time_step.total_seconds()
                # Control min and max depth of particles
                sea_surface_height = self.sea_surface_height()[i] # returns surface elevation at particle positions (>0 above msl, <0 below msl)
                if self.elements.z[i] >= sea_surface_height:
                    self.elements.z[i] = sea_surface_height - 0.01 # set particle z at 0.01m below sea_surface_height
                    self.elements.vertical_movement[i] = - abs(movementr)
                if self.get_config('drift:maximum_depth') is not None:
                    if self.elements.z[i] <= self.get_config('drift:maximum_depth'):
                        self.elements.vertical_movement[i] = + abs(movementr)
            else:
                pass
    
    
    def haliotis_iris_sinking(self):
        ''''Haliotis iris larvae sink durign the first 12 hours of their life, 
        then swim toward the surface'''
        young_haliotis_iris = np.where(self.elements.age_seconds <= 12*3600)[0]
        if len(young_haliotis_iris) > 0:
            self.elements.z[young_haliotis_iris] = self.elements.z[young_haliotis_iris] - 2*abs(self.elements.terminal_velocity[young_haliotis_iris]) * self.time_step.total_seconds()
    
    
    def maximum_depth(self):
            '''Turn around larvae that are going too deep'''
            if self.get_config('drift:maximum_depth') is not None:
                too_deep = np.where(self.elements.z < self.get_config('drift:maximum_depth'))[0]
                if len(too_deep) > 0:
                    self.elements.z[too_deep] = self.get_config('drift:maximum_depth')
    
    
    def increase_age_and_retire(self):  ####So that if max_age_seconds is exceeded particle is flagged as died
            """Increase age of elements, and retire if older than config setting.

               >essentially same as increase_age_and_retire() from basemodel.py, 
               only using a diffrent reason for retiring particles ('died' instead of 'retired')
               .. could probably be removed...
            """
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
                    
    def update(self):
        """Update positions and properties of buoyant particles."""

        # Update element age
        # self.elements.age_seconds += self.time_step.total_seconds()
        # already taken care of in increase_age_and_retire() in basemodel.py
        
        # Horizontal advection
        self.advect_ocean_current()

        # Check for presence in habitat
        if self.get_config('drift:settlement_in_habitat') is True:
            self.interact_with_habitat()
            
        # Additional module for Haliotis iris
        if self.get_config('drift:Haliotis_iris') is True:
            self.haliotis_iris_sinking()

        # Turbulent Mixing or settling-only 
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()  #compute vertical velocities, two cases possible - constant, or same as pelagic egg
            self.vertical_mixing()
        else:  # Buoyancy
            self.update_terminal_velocity()
            self.vertical_buoyancy()

        if self.get_config('drift:active_vertical_swimming') is True:
            self.vertical_swimming()

        self.vertical_advection()
        
        self.maximum_depth()
