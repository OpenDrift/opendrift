#   
#   This file is intended for the rock lobster larave Jasus edwardsii and has been adapted from bivalvelarvae.py and larvalfish.py 
# 
#   It introduces a multi-step dispersal phase based on the development of the lobster larvae: nauplosioma, phyllosoma, and puerulus
# 
# . Nauplosioma stage:  Lobster larvae hatch into this stage and rise to the surface where they metamorphose into phyllosoma. 
#                       Brief stage (<12hrs) that is not represented in this code, but can be implied by the release of larvae 
#                       slightly offshore and at the surface  
#
# . Phyllosoma stage:   Planktonic stage where larvae can control their vertical position in the water column via buoyancy. 
#                       Late stage is characterized by diel vertical migration and slow swimming. 
#                       Rarely found inshore after the stage 5, so we remove phyllosoma larvae found within 20 km 
#                       of the coast 4 to 12 mo after hatching.
#
#   Puerulus stage:     At this stage the lobster larvae have developed horizontal swimming capabilities and 
#                       if close enough to the coast will swim to their settlement habitat. The pueruli have poorly 
#                       developed mouth and rely on stored energy for the final stretch of their dispersal.
# 
#  Authors : Romain Chaput, University of Auckland, Simon Weppe MetOcean Solutions/MetService New Zealand
# 

import numpy as np
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.models.bivalvelarvae import BivalveLarvae
import logging; logger = logging.getLogger(__name__)
from datetime import timezone
from shapely.geometry import Polygon, Point, MultiPolygon,asPolygon # added for settlement in polygon only
import shapely
import random
from sklearn.neighbors import BallTree

# Defining the  element properties from Pelagicegg model
class LobsterLarvaeObj(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs/larvae
    """

    variables = Lagrangian3DArray.add_variables([
                ('light', {'dtype': np.float32,
                             'units': 'ugEm2',
                             'default': 0.}),
                ('u_swim', {'dtype': np.float32,
                             'units': 'm/s',
                             'default': 0.}),
                ('v_swim', {'dtype': np.float32,
                             'units': 'm/s',
                             'default': 0.})])

class LobsterLarvae(BivalveLarvae):
    """
    LobsterLarvae inherits from the bivavel larvae module which includes a minimum age before 
    settling on coast or seafloor can occur, as well as the option to define suitable habitat polygons
    """

    ElementType = LobsterLarvaeObj

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
    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    required_profiles_z_range = [-200, 0] # The depth range (in m) which
                                          # profiles shall cover

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'settled_on_coast': 'red', 'died': 'yellow', 
                     'settled_on_bottom': 'magenta',
                     'settled_on_habitat': 'orange',
                     'swam_too_close_to_shore' : 'black'}

    def __init__(self, *args, **kwargs):
        
        # Calling general constructor of parent class
        super(LobsterLarvae, self).__init__(*args, **kwargs)

        # By default, larvae do not strand when reaching shoreline
        # They are recirculated back to previous position instead
        self.set_config('general:coastline_action', 'previous') # same as bivalve
        # resuspend larvae that reach seabed by default 
        self.set_config('general:seafloor_action', 'lift_to_seafloor') # same as bivalve

        ##add config spec
        # habitat
        # self._add_config({ 'biology:min_settlement_age_seconds': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
        #                    'description': 'minimum age in seconds at which larvae can start to settle on seabed or stick to shoreline',
        #                    'level': self.CONFIG_LEVEL_BASIC}})
        # self._add_config({ 'biology:settlement_in_habitat': {'type': 'bool', 'default': False,
        #                    'description': 'settlement restricted to suitable habitat only',
        #                    'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:direct_orientation_habitat': {'type': 'bool', 'default': False,
                           'description': 'biased correlated random walk toward the nearest habitat',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:max_orient_distance': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'meters',
                           'description': 'maximum detection distance of the habitat for orientation',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:age_beginning_orientation': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'Beginning of the orientation of the larvae in seconds (can correspond to the flexion stage)',
                           'level': self.CONFIG_LEVEL_BASIC}})
        # life stages
        self._add_config({ 'biology:stage_puerulus': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'minimum age of puerulus larvae',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:late_stage_phyllosoma': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'minimum age age of late stage phyllosoma',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:mid_stage_phyllosoma': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'minimum age age of middle stage phyllosoma',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:max_swimming_speed_puerulus': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'centimeters.seconds^-1',
                           'description': 'maximum swimming speed of the puerulus when crusing toward the habitat',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:min_swimming_speed_puerulus': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'centimeters.seconds^-1',
                           'description': 'minimum swimming speed of the puerulus when crusing toward the habitat',
                           'level': self.CONFIG_LEVEL_BASIC}})
        # vertical motion, as in pelagicplankton_moana.py
        self._add_config({ 'biology:vertical_position_daytime': {'type': 'float', 'default': -5.00,'min': -100.0, 'max':0.0, 'units': 'meters negative down',
                           'description': 'the depth a species is expected to inhabit during the day time, in meters, negative down',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:vertical_position_nighttime': {'type': 'float', 'default': -1.00,'min': -100.0, 'max':0.0, 'units': 'meters negative down',
                           'description': 'the depth a species is expected to inhabit during the night time, in meters, negative down',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:vertical_migration_speed_constant': {'type': 'float', 'default': 1.0e-4,'min': 0.0, 'max': 1.0e-3, 'units': 'm/s',
                           'description': 'Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity()',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:maximum_larvae_depth': {'type': 'float', 'default': -100.0,'min': -10000.0, 'max': -1.0, 'units': 'm',
                           'description': 'maximum depth the larvae can swim down to, larvae will swim up if reached',
                           'level': self.CONFIG_LEVEL_BASIC}})                       

  #####################################################################################################################
  # Accessory functions
  #####################################################################################################################
    
    # Haversine formula to compute angles during orientation
    def haversine_angle(self, lon1, lat1, lon2, lat2):
        rlat1 = np.deg2rad(lat1)
        rlat2 = np.deg2rad(lat2)
        rlon1 = np.deg2rad(lon1)
        rlon2 = np.deg2rad(lon2)
        X = np.cos(rlat2)*np.sin(rlon2-rlon1)
        Y = np.cos(rlat1)*np.sin(rlat2)-np.sin(rlat1)*np.cos(rlat2)*np.cos(rlon2-rlon1)
        return np.arctan2(Y,X)
     
    def calculateMaxSunLight(self):
        # Calculates the max sun radiation at given positions and dates (and returns zero for night time)
        # 
        # The method is using the third party library PySolar : https://pysolar.readthedocs.io/en/latest/#
        # 
        # 
        # some other available options:
        # https://pypi.org/project/solarpy/
        # https://github.com/trondkr/pyibm/blob/master/light.py
        # use calclight from Kino Module here  : https://github.com/trondkr/KINO-ROMS/tree/master/Romagnoni-2019-OpenDrift/kino
        # ERcore : dawn and sunset times : https://github.com/metocean/ercore/blob/ercore_opensrc/ercore/lib/suncalc.py
        # https://nasa-develop.github.io/dnppy/modules/solar.html#examples
        # 
        from pysolar import solar
        date = self.time
        date = date.replace(tzinfo=timezone.utc) # make the datetime object aware of timezone, set to UTC
        logger.debug('Assuming UTC time for solar calculations')
        # longitude convention in pysolar, consistent with Opendrift : negative reckoning west from prime meridian in Greenwich, England
        # the particle longitude should be converted to the convention [-180,180] if that is not the case
        sun_altitude = solar.get_altitude(self.elements.lon, self.elements.lat, date) # get sun altitude in degrees
        sun_azimut = solar.get_azimuth(self.elements.lon, self.elements.lat, date) # get sun azimuth in degrees
        sun_radiation = np.zeros(len(sun_azimut))
        # not ideal get_radiation_direct doesnt accept arrays...
        for elem_i,alt in enumerate(sun_altitude):
            sun_radiation[elem_i] = solar.radiation.get_radiation_direct(date, alt)  # watts per square meter [W/m2] for that time of day
        # save compute light for each particle 
        self.elements.light = sun_radiation * 4.6 #Converted from W/m2 to umol/m2/s-1"" - 1 W/m2 ≈ 4.6 μmole.m2/s
        logger.debug('Solar radiation from %s to %s [W/m2]' % (sun_radiation.min(), sun_radiation.max() ) )
        # print(np.min(sun_radiation))
        # print(date)

    #####################################################################################################################
    # IBM: horizontal and vertical movements
    #####################################################################################################################

    def reset_horizontal_swimming(self):
            # Create a  vector for swimming movement
            self.u_swim = np.array([0.0]*len(self.elements.lat))
            self.v_swim = np.array([0.0]*len(self.elements.lon))
            return self.u_swim, self.v_swim

    def direct_orientation_habitat(self):
            """"
            Direct orientation, where the larvae swim toward the nearest reef. 
            Biased correlated random walk toward the nearest habitat. 
            Equations described in Codling et al., 2004 and Staaterman et al., 2012

            ** Algorithm could be vectorized to optimize speed (instead of looping through each particle)

            """
            # Check if the particles are old enough to orient
            old_enough = np.where(self.elements.age_seconds >= self.get_config('biology:age_beginning_orientation'))[0]
            logger.debug('Larvae : direct orientation - %s particles ready for orientation ' % (len(old_enough)))
            if len(old_enough) > 0 :
                # find closest habitat distance, and its index
                habitat_near, habitat_id = self.find_nearest_habitat(self.elements.lon[old_enough],self.elements.lat[old_enough])
                # old_close_enough = np.where( (self.elements.age_seconds >= self.get_config('biology:age_beginning_orientation')) & \
                #                              (habitat_near.ravel()*6371. < self.get_config('biology:max_orient_distance')) )[0]
                close_enough = (habitat_near.ravel()*6371. < self.get_config('biology:max_orient_distance')) # 6371km is earth radius 
                logger.debug('Larvae : direct orientation - moving %s particles towards nearest reef' % (close_enough.sum()))

                cur_dir_rad = self.get_current_direction()

                # looping through each particle - could be vectorized
                for i in range(len(self.elements.lat[old_enough])):
                    if habitat_near[i][0]*6371 > self.get_config('biology:max_orient_distance'):
                        pass
                    else:
                        pt_lon = self.elements.lon[old_enough][i]
                        pt_lat = self.elements.lat[old_enough][i]
                        pt_lon_old = self.previous_lon[old_enough][i]
                        pt_lat_old = self.previous_lat[old_enough][i]
                        # Case where particle is old enough and close enough from habitat to orient
                        # Strength of orientation (depend on distance to the habitat). Eq. 3 Staaterman et al., 2012
                        d = 1 - (habitat_near[i][0]*6371/self.get_config('biology:max_orient_distance')) # 
                        # Compute direction of nearest habitat. See Staaterman et al., 2012
                        theta_pref = - self.haversine_angle(pt_lon, pt_lat, self.centers_habitat[habitat_id[i][0]][0], self.centers_habitat[habitat_id[i][0]][1]) 
                        # Compute current direction from previous timestep
                        # ** note this will include the swimming-induced motions from past timestep as well
                        theta_current = self.haversine_angle(pt_lon_old, pt_lat_old, pt_lon, pt_lat) # from previous positions 
                        # theta_current = cur_dir_rad[old_enough[i]] # from ambient current direction
                        # Mean turning angle, Eq. 2 Staaterman et al., 2012
                        mu = -d * (theta_current - theta_pref)
                        # Define larvae swimming direction
                        # direction uncertainity drawn from Von Mises distribution : First parameter: mu, second parameter: kappa
                        ti  = np.random.vonmises(0, 5) #  (control the uncertainty of orientation)
                        # final larvae swimming direction
                        theta = ti - theta_current - mu
                        # Compute u and v swim velocity component
                        self.elements.u_swim[old_enough[i]] = self.swimming_speed()*np.cos(theta)
                        self.elements.v_swim[old_enough[i]] = self.swimming_speed()*np.sin(theta)

            logger.debug('    %.5f [m/s] =< u_swim <=  %.5f [m/s]' % (np.min(self.elements.u_swim),np.max(self.elements.u_swim))) 
            logger.debug('    %.5f [m/s] =< v_swim <=  %.5f [m/s]' % (np.min(self.elements.u_swim),np.max(self.elements.u_swim))) 
            self.update_positions(self.elements.u_swim , self.elements.v_swim)
    
    def find_nearest_habitat(self,lon,lat):
        # self.ball_centers.query(list(zip(np.deg2rad(self.elements.lat[old_enough]), np.deg2rad(self.elements.lon[old_enough]))), k=1)
        return self.ball_centers.query(list(zip(np.deg2rad(lat), np.deg2rad(lon))), k=1)

    def swimming_speed(self):
            # Compute swimming speed of larvae - no dependence on age 
            # beta distribution to reproduce the findings in Jeffs and Hollands 2000 (speeds from 13 and 22cm/s with outlier at 30.7cm/s) 
            # and Wilkin and Jeffs 2011 (model predictions btwn 13 and 16cm/s) (mean = 1/(1+B/A) => 16.8cm/s here with a=2 and b=8)
            # Not age dependent because puerulus do not develop further until they reach an habitat
            # swim_speed = (Vmin + random.beta*(Vmax-Vmin) ) (* dt ???)
            # 
            # Jeffs and Hollands,2000,Swimming behaviour of the puerulus of the spiny lobster, Jasus edwardsii.Crustaceana 73(7):847-856
            
            Vmin = self.get_config('biology:min_swimming_speed_puerulus')
            Vmax = self.get_config('biology:max_swimming_speed_puerulus')
            swimming_speed = ( Vmin + np.random.beta(2, 8) * (Vmax-Vmin) ) / 100
            return swimming_speed 

            # original code..not sure about multiplication by dt ? 
            # swimming_speed = (self.get_config('biology:min_swimming_speed_puerulus') \
            #    + np.random.beta(2, 8) \
            #    * (self.get_config('biology:max_swimming_speed_puerulus') - self.get_config('biology:min_swimming_speed_puerulus'))) \
            #    * self.time_step.total_seconds()   
   
    def get_current_direction(self):
        ''' returns current direction in the trigonometric convention, as used in direct_orientation_habitat'''
        uu = self.environment['x_sea_water_velocity']
        vv = self.environment['y_sea_water_velocity']
        return np.arctan2(vv,uu)


    
    def update_terminal_velocity(self,Tprofiles=None, Sprofiles=None,z_index=None): 
            ''' Diel vertical migration for late life stages (VIII to XI) only. Modified from pelagicplankton_moana.py developed by Simon Weppe
            Keep all stages within 100m of the surface '''
            
            # Modifies the same variable than update_terminal_velocity(), self.elements.terminal_velocity = W, but using a different algorithm.
            # Larvae are assumed to move to daytime or nighttime vertical positions in the water column, at a constant rate
            # the actual settling is taken care of in vertical_mixing() or vertical_buoyancy() (i.e. from OceanDrift methods)
            # it is expected that larve will go down during day time and up during night time but that is not fixed in the code. 
            # Particles will simply go towards the daytime or nighttime positions.
            # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L80 

            late_stage_phy = self.elements.age_seconds > self.get_config('biology:late_stage_phyllosoma') 
            logger.debug('Larvae : update_terminal_velocity - %s particles in late_stage_phyllosoma' % (late_stage_phy.sum()) )
            if late_stage_phy.any() > 0 :
                self.calculateMaxSunLight() # compute solar radiation at particle positions (using PySolar)
                vertical_velocity = np.abs(self.get_config('biology:vertical_migration_speed_constant'))  # magnitude in m/s 
                z_day = self.get_config('biology:vertical_position_daytime')    #  the depth a species is expected to inhabit during the day time, in meters, negative down') #
                z_night =self.get_config('biology:vertical_position_nighttime') # 'the depth a species is expected to inhabit during the night time, in meters, negative down') #
                ind_day =   np.where(late_stage_phy & (self.elements.light>0) )     #np.where(self.elements.light[late_stage_phy]>0)
                ind_night = np.where(late_stage_phy & (self.elements.light == 0) )  #np.where(self.elements.light[late_stage_phy]==0)
                logger.debug('Using constant migration rate (%s m/s) towards day and night time positions' % (vertical_velocity) )
                logger.debug('%s particles in day time' % (len(ind_day[0])))
                logger.debug('%s particles in night time' % (len(ind_night[0])))
                # for particles in daytime : particles below the daytime position need to go up while particles above the daytime position need to go down
                # (same for for particles in nightime)
                # Note : depth convention is negative down in Opendrift
                # 
                # e.g. z=-5, z_day = -3, below daytime position,  need to go up (terminal_velocity>0) 
                #      diff = (z - z_day) = -2, so w = - np.sign(diff) * vertical_velocity
                self.elements.terminal_velocity[ind_day] = - np.sign(self.elements.z[ind_day] - z_day) * vertical_velocity
                self.elements.terminal_velocity[ind_night] = - np.sign(self.elements.z[ind_night] - z_night) * vertical_velocity
                # print(self.elements.z)
                
                # take into account maximum depth for other larval stages
                early_life = np.where(self.elements.age_seconds < self.get_config('biology:late_stage_phyllosoma'))[0]
                if len(early_life) > 0:
                    logger.debug('%s particles deeper than max depth %s m, setting positive terminal_velocity' % ( len(early_life),self.get_config('biology:maximum_larvae_depth') ))
                    too_deep = np.where(self.elements.z[early_life] < self.get_config('biology:maximum_larvae_depth') )[0]
                    self.elements.terminal_velocity[too_deep] = np.abs(vertical_velocity) # set vertical swimming speed positive so particles go up, in m/s
                
###################################################################################################################
# Pelagic larval duration and mortality
###################################################################################################################

    def phyllosoma_mortality(self):
        ''' Phyllosoma are not found inshore between 4 and 12 months of dispersal (i.e between mid and late Phyllosoma stages)
        This could be due to an increased in predatory pressure closer to the coast,
        therefore, we remove all the phyllosoma larvae that are found within 20km of the coast
        '''
        mid_stage_phyllosoma =  np.where( (self.elements.age_seconds >= self.get_config('biology:mid_stage_phyllosoma')) & \
                                          (self.elements.age_seconds <= self.get_config('biology:late_stage_phyllosoma')) )[0]
        logger.debug('Larvae : checking phyllosoma distance to shore - %s particles in mid to late_stage_phyllosoma' % (len(mid_stage_phyllosoma)))
        if len(mid_stage_phyllosoma) > 0:
            for i in range(len(self.elements.lon[mid_stage_phyllosoma])):
                pt_lon = self.elements.lon[mid_stage_phyllosoma][i]
                pt_lat = self.elements.lat[mid_stage_phyllosoma][i] 
                # Remove phyllosoma found 20km inshore, or any other habitats defined by user (habitat_near[i]*6371 = distance in km)
                lon_circle ,lat_circle = self.get_circle(pt_lon, pt_lat, 20e3) 
                # check if any land points within the 20km radius
                on_land, prof, missing = self.get_environment(['land_binary_mask'], self.time, \
                                                            np.array(lon_circle), np.array(lat_circle), \
                                                            0.0*np.array(lat_circle), None)
                if np.array(on_land.view('f')).any():
                    self.environment.land_binary_mask[mid_stage_phyllosoma[i]] = 8
        if (self.environment.land_binary_mask == 8).any():
            logger.debug('Larvae : removing %s phyllosoma that swam too close to shore' % (self.environment.land_binary_mask == 8).sum() )
            self.deactivate_elements((self.environment.land_binary_mask == 8), reason='swam_too_close_to_shore')


    def get_circle(self,centerLon,centerLat,radius_meters):
        # https://stackoverflow.com/questions/15886846/python-elegant-way-of-finding-the-gps-coordinates-of-a-circle-around-a-certain
        radius_earth_meters = 6378137 
        angle = np.linspace(0,2*np.pi,36)
        dx = radius_meters * np.cos(angle)
        dy = radius_meters * np.sin(angle)
        lat_circle = centerLat + (180 / np.pi) * (dy / radius_earth_meters)
        lon_circle = centerLon + (180 / np.pi) * (dx / radius_earth_meters) / np.cos(centerLat * np.pi / 180)
        # import matplotlib.pyplot as plt
        # plt.ion()
        # plt.plot(centerLon,centerLat,'r.')
        # plt.plot(lon_circle,lat_circle,'g.')
        # import pdb;pdb.set_trace()
        return lon_circle ,lat_circle

        # see also..but slower
        # https://gis.stackexchange.com/questions/289044/creating-buffer-circle-x-kilometers-from-point-using-python/289923
        # 
        # from functools import partial
        # import pyproj
        # from shapely.ops import transform
        # from shapely.geometry import Point
        #  proj_wgs84 = pyproj.Proj('+proj=longlat +datum=WGS84')
        # # Azimuthal equidistant projection
        # local_azimuthal_projection = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
        # project = partial(pyproj.transform,
        #                   pyproj.Proj(local_azimuthal_projection.format(lat=centerLat, lon=centerLon)),
        #                   proj_wgs84)
        # buf = Point(0, 0).buffer(radius_meters)  # distance in metres
        # circle = transform(project, buf).exterior.xy
        # lon_circle = circle[0]
        # lat_circle = circle[1]   

###################################################################################################################
# Update position of the larvae
###################################################################################################################    

    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Horizontal swimming: Puerulus stage swimming to shore
        if self.get_config('biology:direct_orientation_habitat') is True:
            self.direct_orientation_habitat() # self.update_positions in the function direct_orientation_habitat

        if False:    
            # Advect particles due to surface wind drag,
            # according to element property wind_drift_factor
            self.advect_wind()

            # Stokes drift
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

        # Mortality due to shore proximity (<20km)
        # slow, may be more relevant to apply in post-processing 
        if True: 
            self.phyllosoma_mortality()
            
