# This file is intended for the orientation of fish larvae and developed from the models bivalvelarvae.py and larvalfish.py 
# 
# This code introduces different orientation behaviors for fish larvae:
# 
# 1. Direct orientation: Competent ('old enough to orient') larvae swim in the direction of the nearest habitat 
#    (represented by polygons in a shapefile) if they are close enough to detect it (controlled by user input). 
#    Based on Staaterman et al., 2012 and Codling et al., 2004.
# 
# 2. Rheotaxis orientation: Countercurrent swimming. Competent larvae swim against the direction of the current. 
#    Larvae swim at lower or equal speed as the currents, depending on their swimming speed. 
#    Larvae cannot swim faster than the currents.
#
# 3. Cardinal orientation: Competent larvae swim in a pre-determine direction (controlled by user input).
# 
# 4. Continuous orientation: 
#    Type 1 => Rheotaxis orientation when larvae are offshore, then reef orientation when they are close 
#    enough to the reefs AND ready to settle. 
#    OR
#    Type 2=> Cardinal orientation then reef orientation. User input cardinal heading, distance of 
#    reef detection, and beginning of orientation.
#
# Common to all orientation behavior: horizontal swimming speed developing with age of the larvae. 
# Coded following the relationship exposed in Fisher and Bellwood, 2003.
# 
# This code also introduces Ontogenetic Vertical Migration (OVM) where larvae swim to different depth following 
# their development (see Irisson et al., 2010 for distribution examples).
# 
# Finally, the settlement of larvae can be restricted to a specific habitat (polygons in shapefile).
#
# Sources for the algorithm of orientation: 
# 
# Codling et al., 2004, Staaterman et al., 2012, and Fisher and Bellwood, 2003. See PhD dissertation R Chaput 2020
# 
# Chaput et al., 2021. Uncertainty propagation into numerical biophysical model of fish larval connectivity in the Florida Keys
# submitted 
# 
# See example_fish_larvae_orientation.py for an example script to run the code
#  
# Author : Romain Chaput
# 
#  Under development - first draft: 25/06/2021 - Last update: 17/08/2021
# 
#  Edits - Simon Weppe 11/2021 
# 

import numpy as np
from opendrift.models.oceandrift import OceanDrift,Lagrangian3DArray
from opendrift.models.bivalvelarvae import BivalveLarvae
import logging; logger = logging.getLogger(__name__)
from datetime import timezone
from shapely.geometry import Polygon, Point, MultiPolygon, asPolygon # added for settlement in polygon only
import shapely
import random
from sklearn.neighbors import BallTree


# Defining the  element properties from Pelagicegg model
class FishLarvaeObj(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs/larvae
    """

    variables = Lagrangian3DArray.add_variables([
                ('flexion_stage', {'dtype': np.int64,
                             'units': '',
                             'default': 0}),
                ('u_swim', {'dtype': np.float32,
                             'units': 'm/s',
                             'default': 0.}),
                ('v_swim', {'dtype': np.float32,
                             'units': 'm/s',
                             'default': 0.})]) # pre_flexion,flexion,post_flexion 0,1,2 - to implement
    
class FishLarvaeOrient(BivalveLarvae,OceanDrift):
    """particle trajectory model based on the OpenDrift framework.
    """

    ElementType = FishLarvaeObj # consider adding a specific element where life stages are saved such as hatched, swimming, etc...

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

    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-200, 0]

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'settled_on_coast': 'red', 'died': 'yellow', 
                     'settled_on_bottom': 'magenta', 
                     'settled_on_habitat': 'orange'}

    def __init__(self, *args, **kwargs):
        
        # Calling general constructor of parent class
        super(FishLarvaeOrient, self).__init__(*args, **kwargs)

        # By default, larvae do not strand when reaching shoreline. 
        # They are recirculated back to previous position instead
        self.set_config('general:coastline_action', 'previous')

        # resuspend larvae that reach seabed by default 
        self.set_config('general:seafloor_action', 'lift_to_seafloor')

        ##add config spec
        self._add_config({ 'biology:orientation': {'type': 'enum', 'default': 'none',
                           'enum': ['none', 'direct', 'rheotaxis', 'cardinal','continuous_1','continuous_2'],
                           'description':  '"direct": biased correlated random walk toward the nearest habitat; \
                                            "rheotaxis": swim against the currents; \
                                            "cardinal": swim in a pre-determined direction; \
                                            "continuous_1": mix of rheotaxis and direct orientation; \
                                            "continuous_2": mix of cardinal and direct orientation; \
                                            "none": no orientation behavior.',
                            'level': self.CONFIG_LEVEL_ADVANCED}})
        # self._add_config({ 'biology:min_settlement_age_seconds': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
        #                    'description': 'minimum age in seconds at which larvae can start to settle on habitat, or seabed or stick to shoreline',
        #                    'level': self.CONFIG_LEVEL_BASIC}})
        # self._add_config({ 'biology:settlement_in_habitat': {'type': 'bool', 'default': False,
        #                    'description': 'settlement restricted to suitable habitat only, specified by shapefile or txt file with nan-delimited polys',
        #                    'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:OVM': {'type': 'bool', 'default': False,
                           'description': 'Ontogenetic Vertical Migration',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:cardinal_heading': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 360.0, 'units': 'degrees North',
                           'description': 'cardinal orientation of the larvae during dispersal, nautical convention (0.0 degree = North, positive clockwise)',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:max_orient_distance': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'meters',
                           'description': 'maximum detection distance of the habitat for orientation',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:settle_swimming_speed': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'centimeters.seconds^-1',
                           'description': 'maximum swimming speed of the larvae at settlement stage',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:hatch_swimming_speed': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'centimeters.seconds^-1',
                           'description': 'minimum swimming speed of the larvae at hatching',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:vertical_migration_speed_constant': {'type': 'float', 'default': None,'min': 0.0, 'max': 100.0, 'units': 'm/s',
                           'description': 'Constant vertical migration rate (m/s), if None, uses values from update_terminal_velocity()',
                           'level': self.CONFIG_LEVEL_BASIC}})            
        self._add_config({ 'biology:maximum_larvae_depth': {'type': 'float', 'default': None,'min': -10000.0, 'max': -1.0, 'units': 'm',
                           'description': 'maximum depth of the larvae',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:beginning_orientation': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'Beginning of the orientation of the larvae in seconds (can correspond to the flexion stage)',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:pre_flexion': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'age of the pre-flexion stage transition in seconds',
                           'level': self.CONFIG_LEVEL_BASIC}})    
        self._add_config({ 'biology:flexion': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'Flexion stage of the larvae in seconds',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:post_flexion': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'age of the post-flexion stage transition in seconds',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:depth_early_stage': {'type': 'float', 'default': 0.0,'max': 0.0, 'min': -1.0e10, 'units': 'meters',
                           'description': 'mean depth of larvae during the first 24 hours',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:depth_pre_flexion': {'type': 'float', 'default': 0.0,'max': 0.0, 'min': -1.0e10, 'units': 'meters',
                           'description': 'mean depth of larvae during the pre-flexion stage',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:depth_flexion': {'type': 'float', 'default': 0.0,'max': 0.0, 'min': -1.0e10, 'units': 'meters',
                           'description': 'mean depth of larvae during the flexion stage',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:depth_post_flexion': {'type': 'float', 'default': 0.0,'max': 0.0, 'min': -1.0e10, 'units': 'meters',
                           'description': 'mean depth of larvae during the post-flexion stage',
                           'level': self.CONFIG_LEVEL_BASIC}})

    #####################################################################################################################
    # accesory functions
    ##############mid_mid_stage_phyllosomastage_phyllosoma#######################################################################################################
     
    # Haversine formula to compute angles during orientation
    def haversine_angle(self, lon1, lat1, lon2, lat2):
        rlat1 = np.deg2rad(lat1)
        rlat2 = np.deg2rad(lat2)
        rlon1 = np.deg2rad(lon1)
        rlon2 = np.deg2rad(lon2)
        X = np.cos(rlat2)*np.sin(rlon2-rlon1)
        Y = np.cos(rlat1)*np.sin(rlat2)-np.sin(rlat1)*np.cos(rlat2)*np.cos(rlon2-rlon1)
        return np.arctan2(Y,X)
    
    #####################################################################################################################
    # Definition of habitat ...>> inherits from bivalvelarvae.py
    #####################################################################################################################
    
    #####################################################################################################################
    # IBM: horizontal and vertical movements
    #####################################################################################################################
            
    def direct_orientation_habitat(self):
            """"
            Direct orientation, where the larvae swim toward the nearest reef. 
            Biased correlated random walk toward the nearest habitat. 
            Equations described in Codling et al., 2004 and Staaterman et al., 2012

            ** algorithm could be vectorized to optimize speed (instead of looping through each particle)

            """
            # Check if the particles are old enough to orient
            old_enough = np.where(self.elements.age_seconds >= self.get_config('biology:beginning_orientation'))[0]
            logger.debug('Larvae : direct orientation - %s particles ready for orientation ' % (len(old_enough)))
            if len(old_enough) > 0 :
                habitat_near, habitat_id = self.ball_centers.query(list(zip(np.deg2rad(self.elements.lat[old_enough]), np.deg2rad(self.elements.lon[old_enough]))), k=1)
                # old_close_enough = np.where( (self.elements.age_seconds >= self.get_config('biology:beginning_orientation')) & \
                #                              (habitat_near.ravel()*6371. < self.get_config('biology:max_orient_distance')) )[0]
                close_enough = (habitat_near.ravel()*6371. < self.get_config('biology:max_orient_distance')) 
                logger.debug('Larvae : direct orientation - moving %s particles towards nearest reef' % (close_enough.sum()))
                # looping through each particles - could be vectorized
                for i in range(len(self.elements.lat[old_enough])):
                    if habitat_near[i][0]*6371 > self.get_config('biology:max_orient_distance'):
                        pass
                    else:
                        pt_lon = self.elements.lon[old_enough][i]
                        pt_lat = self.elements.lat[old_enough][i]
                        pt_lon_old = self.previous_lon[old_enough][i]
                        pt_lat_old = self.previous_lat[old_enough][i]
                        # Case where particle close enough and old enough to orient
                        # Strength of orientation (depend on distance to the habitat)
                        d = 1 - (habitat_near[i][0]*6371/self.get_config('biology:max_orient_distance'))
                        # Compute direction of nearest habitat. See Staaterman et al., 2012
                        theta_pref = - self.haversine_angle(pt_lon, pt_lat, self.centers_habitat[habitat_id[i][0]][0], self.centers_habitat[habitat_id[i][0]][1]) 
                        # Compute direction from previous timestep
                        theta_current = self.haversine_angle(pt_lon_old, pt_lat_old, pt_lon, pt_lat)
                        # Mean turning angle
                        mu = -d * (theta_current - theta_pref)
                        # New direction randomly selected in Von Mises distribution
                        ti  = np.random.vonmises(0, 5) # First parameter: mu, second parameter: kappa (control the uncertainty of orientation) 
                        theta = ti - theta_current - mu
                    
                        # Compute u and v velocity
                        self.elements.u_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.cos(theta)
                        self.elements.v_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.sin(theta) 
                            
            self.update_positions(self.elements.u_swim , self.elements.v_swim)
    
    
    def cardinal_orientation(self):
            """"
            Cardinal orientation, where the larvae use a magnetic, internal compass to swim
            toward a pre-determined heading
            """
            # Check if the particles are old enough to orient
            old_enough = np.where(self.elements.age_seconds >= self.get_config('biology:beginning_orientation'))[0]
            logger.debug('Larvae : cardinal orientation - moving %s particles towards cardinal_heading' % (len(old_enough)))

            if len(old_enough) > 0 :
                # Compute heading
                cardinal_heading_cartesian = (90 - self.get_config('biology:cardinal_heading')) # convert cardinal_heading to cartesian
                thetaCard = np.deg2rad(cardinal_heading_cartesian) 
                for i in range(len(self.elements.lat[old_enough])):
                    # Computing preferred direction
                    ti  = np.random.vonmises(0, 5)
                    theta = thetaCard + ti                
                    # Compute u and v velocity
                    self.elements.u_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.cos(theta)
                    self.elements.v_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.sin(theta)
            
            self.update_positions(self.elements.u_swim , self.elements.v_swim)
        
    def rheotaxis_orientation(self):
            """"
            Rheotaxis orientation, where the larvae swim against the local current at a
            speed slower or equal to that of the currents. 

            ** algorithm could be vectorized to optimize speed

            """            
            # Check if the particles are old enough to orient
            old_enough = np.where(self.elements.age_seconds >= self.get_config('biology:beginning_orientation'))[0]
            logger.debug('Larvae : rheotaxis orientation - moving %s particles against current' % (len(old_enough)))

            if len(old_enough) > 0 :            
                for i in range(len(self.elements.lat[old_enough])):
                    # Compute randomness of direction
                    ti  = np.random.vonmises(0, 5) # Samples from a von Mises distribution on the interval [-pi, pi].
                    #Compute rheotaxis heading i.e. against the current speed at particle position (180deg difference)
                    current_dir = np.arctan2(self.environment.y_sea_water_velocity[old_enough[i]], self.environment.x_sea_water_velocity[old_enough[i]])  # Direction of current
                    thetaRheo = np.mod(current_dir+np.pi,2*np.pi) # against current direction
                    theta = thetaRheo + ti  # add some uncertainty                  
                    # Compute current speed absolute value
                    uv = np.sqrt(self.environment.x_sea_water_velocity[old_enough[i]]**2 + self.environment.y_sea_water_velocity[old_enough[i]]**2)
                    # Compute norm of swimming speed
                    norm_swim = np.abs(self.swimming_speed(self.elements.age_seconds[old_enough][i]))
                    if norm_swim < uv:
                        # Compute u and v velocity
                        self.elements.u_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.cos(theta)
                        self.elements.v_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.sin(theta)
                    else:
                        # Note that if particle swim against current at same speed they will become static
                        self.elements.u_swim[old_enough[i]] = uv * np.cos(theta)
                        self.elements.v_swim[old_enough[i]] = uv * np.sin(theta)
            self.update_positions(self.elements.u_swim , self.elements.v_swim)

    def mix_orientation(self):
            """
            Continuous orientation of the larvae using either rheotaxis orientation and direct orientation when getting closer to the reefs, 
            or cardinal orientation and direct orientation when getting closer to the reefs

            continuous_1 :  the larvae swim first against the currents and switch to Direct Orientation when they approach the reefs
            continuous_2 :  the larvae orient with Cardinal Orientation first and switch to Direct orientation when they approach the reefs

            ** algorithm could be vectorized to optimize speed

                """
            # Check if the particles are old enough to orient toward the reefs and settle
            old_enough = np.where(self.elements.age_seconds >= self.get_config('biology:beginning_orientation'))[0]
            logger.debug('Larvae : mix orientation - moving %s particles according to option: %s' % (len(old_enough),self.get_config('biology:orientation')))

            if len(old_enough) > 0 :
                for i in range(len(self.elements.lat[old_enough])):
                    if self.get_config('biology:orientation')=='continuous_1':
                        # Compute randomness of direction
                        ti  = np.random.vonmises(0, 5)
                        #Compute rheotaxis heading
                        current_dir = np.arctan2(self.environment.y_sea_water_velocity[old_enough[i]], self.environment.x_sea_water_velocity[old_enough[i]])  # Direction of current
                        thetaRheo = np.mod(current_dir+np.pi,2*np.pi) # against current direction
                        # thetaRheo = -np.arctan2(self.environment.y_sea_water_velocity[old_enough[i]], self.environment.x_sea_water_velocity[old_enough[i]])
                        theta = thetaRheo + ti
                
                        # Compute current speed absolute value
                        uv = np.sqrt(self.environment.x_sea_water_velocity[old_enough[i]]**2 + self.environment.y_sea_water_velocity[old_enough[i]]**2)
                        # Compute norm of swimming speed
                        norm_swim = np.abs(self.swimming_speed(self.elements.age_seconds[old_enough][i]))
                
                        if norm_swim < uv:
                                # Compute u and v velocity
                                self.elements.u_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.cos(theta)
                                self.elements.v_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.sin(theta)
                        else:
                                self.elements.u_swim[old_enough[i]] = uv * np.cos(theta)
                                self.elements.v_swim[old_enough[i]] = uv * np.sin(theta)
                                
                    if self.get_config('biology:orientation')=='continuous_2':
                        # Computing preferred direction
                        cardinal_heading_cartesian = (90 - self.get_config('biology:cardinal_heading')) # convert cardinal_heading to cartesian
                        thetaCard = np.deg2rad(cardinal_heading_cartesian)
                        # thetaCard = np.deg2rad(self.get_config('biology:cardinal_heading'))
                        ti  = np.random.vonmises(0, 5)
                        theta = thetaCard + ti  
                        # Compute u and v velocity
                        self.elements.u_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.cos(theta)
                        self.elements.v_swim[old_enough[i]] = self.swimming_speed(self.elements.age_seconds[old_enough][i])*np.sin(theta)
            # Check if larvae are old enough to go back to the reef to settle
            looking_for_reef = np.where(self.elements.age_seconds >= self.get_config('biology:min_settlement_age_seconds'))[0]
            if len(looking_for_reef) > 0:
                # Check if larvae are close enough to detect the reef
                habitat_near, habitat_id = self.ball_centers.query(list(zip(np.deg2rad(self.elements.lat[looking_for_reef]), np.deg2rad(self.elements.lon[looking_for_reef]))), k=1)
                for i in range(len(self.elements.lat[looking_for_reef])):
                        if habitat_near[i][0]*6371 < self.get_config('biology:max_orient_distance'):
                            pt_lon = self.elements.lon[looking_for_reef][i]
                            pt_lat = self.elements.lat[looking_for_reef][i]
                            pt_lon_old = self.previous_lon[looking_for_reef][i]
                            pt_lat_old = self.previous_lat[looking_for_reef][i]
                            # Case where particle close enough and old enough to orient
                            # Strength of orientation (depend on distance to the habitat)
                            d = 1 - (habitat_near[i][0]*6371/self.get_config('biology:max_orient_distance'))
                            # Compute direction of nearest habitat. See Staaterman et al., 2012
                            theta_pref = - self.haversine_angle(pt_lon, pt_lat, self.centers_habitat[habitat_id[i][0]][0], self.centers_habitat[habitat_id[i][0]][1]) 
                            # Compute direction from previous timestep
                            theta_current = self.haversine_angle(pt_lon_old, pt_lat_old, pt_lon, pt_lat)
                            # Mean turning angle
                            mu = -d * (theta_current - theta_pref)
                            # New direction randomly selected in Von Mises distribution
                            ti  = np.random.vonmises(0, 5) # First parameter: mu, second parameter: kappa (control the uncertainty of orientation) 
                            theta = ti - theta_current - mu
                        
                            # Compute u and v velocity
                            self.elements.u_swim[looking_for_reef[i]] = self.swimming_speed(self.elements.age_seconds[looking_for_reef][i])*np.cos(theta)
                            self.elements.v_swim[looking_for_reef[i]] = self.swimming_speed(self.elements.age_seconds[looking_for_reef][i])*np.sin(theta)
                        
            self.update_positions(self.elements.u_swim , self.elements.v_swim)
    
    
    def swimming_speed(self, age):
            ''' Compute horizontal swimming speed of the larvae [m/s]
            Presented in Fisher and Bellwood (2003) and used in Staaterman et al. (2012)

            age in seconds
            ''' 
            S_hatch = self.get_config('biology:hatch_swimming_speed') # m/s
            S_settle = self.get_config('biology:settle_swimming_speed') # m/s
            PLD = self.get_config('drift:max_age_seconds') # average pelagic larval duration
            hor_swimming_speed = ( S_hatch + (S_settle - S_hatch) ** (np.log10(age) / np.log10(PLD)) ) / 100 # note x^ab = (x^a)^b to reduce Staaterman equation
            # hor_swimming_speed = (self.get_config('biology:hatch_swimming_speed') + (self.get_config('biology:settle_swimming_speed') - self.get_config('biology:hatch_swimming_speed'))    ** (np.log(age)/np.log(self.get_config('drift:max_age_seconds'))) )    / 100
            logger.debug('Larvae - computing horizontal swimming speed :  %s m/s < swimming speed < %s m/s ' % (np.min(hor_swimming_speed),np.max(hor_swimming_speed)))
            return hor_swimming_speed
        
    def update_terminal_velocity(self):
            '''
            Ontogenetic Vertical Migration. 


            Modified from pelagicplankton_moana.py developed by Simon Weppe
            Modifies self.elements.terminal_velocity, only if self.get_config('biology:OVM') is True
            Particles go towards their preferred depth.

            '''

            if self.get_config('biology:OVM') is True:
                vertical_velocity = np.abs(self.get_config('biology:vertical_migration_speed_constant'))  # magnitude in m/s 
                early_stage = np.where(self.elements.age_seconds < self.get_config('biology:pre_flexion'))[0]
                if len(early_stage) > 0 :
                    self.elements.terminal_velocity[early_stage] = - np.sign(self.elements.z[early_stage] - self.get_config('biology:depth_early_stage')) * vertical_velocity
                
                pre_flexion = np.where((self.elements.age_seconds >= self.get_config('biology:pre_flexion')) & (self.elements.age_seconds < self.get_config('biology:flexion')))[0]
                if len(pre_flexion) > 0 :
                    self.elements.terminal_velocity[pre_flexion] = - np.sign(self.elements.z[pre_flexion] - self.get_config('biology:depth_pre_flexion')) * vertical_velocity
                
                flexion = np.where((self.elements.age_seconds >= self.get_config('biology:flexion')) & (self.elements.age_seconds < self.get_config('biology:post_flexion')))[0]
                if len(flexion) > 0 :
                    self.elements.terminal_velocity[flexion] = - np.sign(self.elements.z[flexion] - self.get_config('biology:depth_flexion')) * vertical_velocity
                    
                post_flexion = np.where(self.elements.age_seconds >= self.get_config('biology:post_flexion'))[0]
                if len(post_flexion) > 0 :
                    self.elements.terminal_velocity[post_flexion] = - np.sign(self.elements.z[post_flexion] - self.get_config('biology:depth_post_flexion')) * vertical_velocity
                
                logger.debug('Larvae - computing vertical swimming speed :  %s cm/s < swimming speed < %s cm/s ' % (np.min(self.elements.terminal_velocity),np.max(self.elements.terminal_velocity)))
            else :
                pass
    def maximum_depth(self):
            '''Turn around larvae that are going too deep'''
            if self.get_config('biology:maximum_larvae_depth') is not None:
                too_deep = np.where(self.elements.z < self.get_config('biology:maximum_larvae_depth'))[0]
                if len(too_deep) > 0:
                    self.elements.z[too_deep] = self.get_config('biology:maximum_larvae_depth')    
                

###################################################################################################################
# Pelagic larval duration and mortality
###################################################################################################################

    def larval_mortality(self):
        ''' Mortality depending on environmental parameters.
        To implement in the future
        '''
        pass
                       
###################################################################################################################
# Update position of the larvae
###################################################################################################################    

    def update(self):
        """Update positions and properties of buoyant particles."""

        ## Horizontal advection
        self.advect_ocean_current() # Independent from age

        if False:    
            # Advect particles due to surface wind drag,
            # according to element property wind_drift_factor
            self.advect_wind()

            # Stokes drift
            self.stokes_drift()

        # Orientation behaviors
        if self.get_config('biology:orientation')=='none':
            pass
        else:
            self.reset_horizontal_swimming()
            if self.get_config('biology:orientation')=='direct':
                self.direct_orientation_habitat() # orientation toward the nearest reef
            if self.get_config('biology:orientation')=='rheotaxis':
                self.rheotaxis_orientation() # orientation against the currents
            if self.get_config('biology:orientation')=='cardinal':
                self.cardinal_orientation() # orientation toward pre-determined direction
            if self.get_config('biology:orientation')=='continuous_1':
                self.mix_orientation() # continuous orientation using rheotaxis and direct orientation
            if self.get_config('biology:orientation')=='continuous_2':
                self.mix_orientation() # continuous orientation using cardinal and direct orientation
        
        # Vertical advection
        self.vertical_advection()   

        # Turbulent Mixing or settling-only 
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity() # includes the OVM if specified
            self.vertical_mixing()
        else:  # Buoyancy
            self.update_terminal_velocity() # includes the OVM if specified
            self.vertical_buoyancy()

        # check for max depth 
        self.maximum_depth()
        
        #** settlement_in_habitat handled internallyt in interact_with_seafloor() 
            
        # Mortality
        self.larval_mortality()
            
