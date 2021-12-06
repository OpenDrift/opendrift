# This file has been adapted from oceandrift.py and introduces a settlement competency period, after which the particle
# can settle on the seafloor, the shoreline or within user-defined habitat polygons
# 
#       Added config setting for 'min_settlement_age_seconds' which defines the minimum age after which  settlement can occur.
#       Added config setting for 'settlement_in_habitat'
# 
# Settlement behaviour : 
# 
# If 'settlement_in_habitat' is False:
#       - coast settlement occurs after a competent particle interacts with the coastline    i.e. comptent after age>min_settlement_age_seconds
#       - seafloor settlement occurs when a competent particle makes contact with the bottom i.e. comptent after age>min_settlement_age_seconds
# If 'settlement_in_habitat' is True:
#       'seafloor settlement occurs when a competent particle makes contact with user-defined habitat polygon(s) i.e. comptent after age>min_settlement_age_seconds          
# 
# The model overloads the seabed and coastline interaction subfunctions from basemodel.py :
#    > lift_elements_to_seafloor()
#    > interact_with_coastline()
# 
#  Authors : Craig Norrie, Simon Weppe, Romain Chaput
# 

import numpy as np
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from shapely.geometry import Polygon, Point, MultiPolygon,asPolygon # added for settlement in polygon only
import shapely
import shapely.vectorized
import random
from sklearn.neighbors import BallTree
import logging; logger = logging.getLogger(__name__)


class BivalveLarvae(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Under construction.
    """

    ElementType = Lagrangian3DArray

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
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
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
    required_profiles_z_range = [-200, 0]

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'settled_on_coast': 'red', 'died': 'yellow', 
                     'settled_on_bottom': 'magenta',
                     'settled_on_habitat': 'orange'}

    def __init__(self, *args, **kwargs):
        
        # Calling general constructor of parent class
        super(BivalveLarvae, self).__init__(*args, **kwargs)

        # By default, larvae do not strand when reaching shoreline. 
        # They are recirculated back to previous position instead
        self._set_config_default('general:coastline_action', 'previous')
        # resuspend larvae that reach seabed by default
        self._set_config_default('general:seafloor_action', 'lift_to_seafloor')  #['previous', 'deactivate', 'lift_to_seafloor']

        # add config spec
        self._add_config({ 'biology:min_settlement_age_seconds': {'type': 'float', 'default': 0.0,'min': 0.0, 'max': 1.0e10, 'units': 'seconds',
                           'description': 'minimum age in seconds at which larvae can start to settle on seabed or stick to shoreline)',
                           'level': self.CONFIG_LEVEL_BASIC}})
        self._add_config({ 'biology:settlement_in_habitat': {'type': 'bool', 'default': False,
                           'description': 'settlement restricted to suitable user-defined habitat only ',
                           'level': self.CONFIG_LEVEL_BASIC}})

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


    def surface_stick(self):
        '''Keep particles just below the surface.
           (overloads the OceanDrift version to allow for possibly time-varying
           sea_surface_height)
        '''
        
        sea_surface_height = self.sea_surface_height() # returns surface elevation at particle positions (>0 above msl, <0 below msl)
        
        # keep particle just below sea_surface_height (self.elements.z depth are negative down)
        surface = np.where(self.elements.z >= sea_surface_height)
        if len(surface[0]) > 0:
            self.elements.z[surface] = sea_surface_height[surface] -0.01 # set particle z at 0.01m below sea_surface_height

    #####################################################################################################################
    # Definition of habitat
    #####################################################################################################################
      
    def add_settlement_habitat(self, shapefile_location):
        """Suitable habitat polygons in a shapefile.shp , or nan-delimited polygons in textfile.txt """
        # 
        # remove dependency to Fiona, and use cartopy instead
        # allow for input of nan-delimited polys from text file as well
        polyList = []
        self.centers_habitat = []
        rad_centers = []
        if '.shp' in shapefile_location :
            from cartopy import io
            shp = io.shapereader.Reader(shapefile_location)
            for shape_i in shp.geometries() : 
                polyList.append(shape_i) # Compile polygon in a list
                self.centers_habitat.append(shape_i.centroid.coords[0]) # Compute centroid and return a [lon, lat] list
        elif '.txt' in  shapefile_location:            
            polys = np.loadtxt(shapefile_location)  # nan-delimited closed polygons for land and islands
            # make sure that start and end lines are [nan,nan] as well
            if not np.isnan(polys[0, 0]):
                polys = np.vstack(([np.nan, np.nan], polys))
            if not np.isnan(polys[-1, 0]):
                polys = np.vstack((polys, [np.nan, np.nan]))
            id_nans = np.where(np.isnan(polys[:, 0]))[0]
            for cnt, _id_i in enumerate(id_nans[:-1]):
                # The shapely.geometry.asShape() family of functions can be used to wrap Numpy coordinate arrays
                # https://shapely.readthedocs.io/en/latest/manual.html
                shape_i = asPolygon(polys[id_nans[cnt] + 1:id_nans[cnt + 1] - 1, :])
                polyList.append(shape_i)
                self.centers_habitat.append(shape_i.centroid.coords[0]) # Compute centroid and return a [lon, lat] list

        for poly in range(len(self.centers_habitat)):
            rad_centers.append([np.deg2rad(self.centers_habitat[poly][1]),np.deg2rad(self.centers_habitat[poly][0])])
        self.multiShp = MultiPolygon(polyList).buffer(0) # Aggregate polygons in a MultiPolygon object and buffer to fuse polygons and remove errors
        # save as shapely prepared geometry for efficient in-poly checks 
        self.habitat_mask = shapely.prepared.prep(MultiPolygon(polyList).buffer(0)) # same as in landmask custom
        self.ball_centers = BallTree(rad_centers, metric='haversine') # Create a Ball Tree with the centroids for faster computation

    #################################################################################################################################
    # SEAFLOOR AND SEABED INTERACTION
    #################################################################################################################################
    def interact_with_seafloor(self): 
        """
        Seafloor interaction according to configuration setting
            Method interaction_with_seafloor() called in vertical_buoyancy(), which is called in update() function       
        """
        # 
        # This function will overloads the version in basemodel.py

        if self.get_config('biology:settlement_in_habitat') is False : 
            # if no specific user-defined habitat, deactivate particles below seabed and older than min_settlement_age_seconds
            self.interact_with_seafloor_bivalve()
        else :
            # if specific user-defined habitat, checks if particles are within habitat polygons             
            self.interact_with_habitat()

        # original code
        # https://github.com/OpenDrift/opendrift/blob/4dbd9a607fe23e64dcbf9fd05905af8713dc74d1/opendrift/models/basemodel.py#L613 
    
    def interact_with_seafloor_bivalve(self):
            '''
            Upon touching the seafloor, the bivalve particle settles only if older than min_settlement_age_seconds
            '''

            if self.num_elements_active() == 0:
                return
            # if 'sea_floor_depth_below_sea_level' not in self.priority_list:
            #     return
            sea_floor_depth = self.sea_floor_depth()
            below = np.where(self.elements.z < -sea_floor_depth)[0]
            if len(below) == 0:
                    logger.debug('No elements hit seafloor.')
                    return 
            below_and_older = np.logical_and(self.elements.z < -sea_floor_depth, 
                self.elements.age_seconds >= self.get_config('biology:min_settlement_age_seconds'))
            below_and_younger = np.logical_and(self.elements.z < -sea_floor_depth, 
                self.elements.age_seconds < self.get_config('biology:min_settlement_age_seconds'))
            
            # Move all elements younger back just above seafloor 
            # (could rather be moved back to previous if relevant? )
            self.elements.z[np.where(below_and_younger)] = -sea_floor_depth[np.where(below_and_younger)] + 0.2
            
            if (below_and_older.any()) & (self.get_config('biology:settlement_in_habitat') is False) :
                # deactivate elements that were both below and older
                self.deactivate_elements(below_and_older ,reason='settled_on_bottom')
                logger.debug('%s elements hit seafloor, %s were older than %s sec. and were deactivated, %s were lifted back to seafloor' \
                    % (len(below),below_and_older.sum(),self.get_config('biology:min_settlement_age_seconds'),below_and_younger.sum()) )

    def interact_with_habitat(self):
            '''
            The bivalve particle settles only if older than min_settlement_age_seconds and within user-defined habitat polygon
 
            To be defined : do we need the particle to be touching seabed AND within habitat or within habitat only?
            '''
            
            if self.num_elements_active() == 0:
                return
            # if 'sea_floor_depth_below_sea_level' not in self.priority_list:
            #     return
            sea_floor_depth = self.sea_floor_depth()
            # below = np.where(self.elements.z < -sea_floor_depth)[0]

            # if len(below) == 0:
            #         logger.debug('No elements hit seafloor.')
            #         return 

            # check if there are any particles old enough to settle and within habitat
            #
            # ** NOTE ** vertical position of particle not taken into account (can be updated)
            old_and_in_habitat = np.logical_and(self.elements.age_seconds >= self.get_config('biology:min_settlement_age_seconds'), 
                                                shapely.vectorized.contains(self.habitat_mask, self.elements.lon, self.elements.lat))
            below = self.elements.z < -sea_floor_depth

            # Move all elements younger than min_settlement_age_seconds, that tocuhed seabed back to just above seafloor 
            self.elements.z[np.where( (~old_and_in_habitat) & (below) )] = -sea_floor_depth[np.where((~old_and_in_habitat) & (below))] + 0.2
            
            if old_and_in_habitat.any():
                self.deactivate_elements(old_and_in_habitat, reason='settled_on_habitat')
                logger.debug('%s elements reached habitat and were older than %s sec. and were deactivated' \
                             % (old_and_in_habitat.sum(),self.get_config('biology:min_settlement_age_seconds')))


    def interact_with_coastline(self,final = False): 
        """Coastline interaction according to configuration setting
           (overloads the interact_with_coastline() from basemodel.py)
           
           The method checks for age of particles that intersected coastlines:

            If no specified habitat :
                - if age_particle < min_settlement_age_seconds : move larvaes back to previous wet position
                - if age_particle > min_settlement_age_seconds : larvaes become stranded and will be deactivated.
            If specified habitat
                - move larvaes back to previous wet position (note one can specify coast or nearshore areas in habitat polygons if relevant)

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
            logger.debug('No elements hit coastline.')
        else: 
            if (self.get_config('biology:min_settlement_age_seconds') == 0.0) | (self.get_config('biology:settlement_in_habitat') is True) :
                # No minimum age input, or user specified some habitat polygon(s)
                # set back to previous "wet" positions (same as in interact_with_coastline() from basemodel.py)
                logger.debug('%s elements hit coastline, '
                          'moving back to water' % len(on_land))
                on_land_ID = self.elements.ID[on_land]
                self.elements.lon[on_land] = \
                    np.copy(self.previous_lon[on_land_ID - 1])
                self.elements.lat[on_land] = \
                    np.copy(self.previous_lat[on_land_ID - 1])
                self.environment.land_binary_mask[on_land] = 0
            else:
                # Minimum age before settling was input, and is not zero
                # check age of particle versus min_settlement_age_seconds and strand/recirculate accordingly
                on_land_and_younger = np.where((self.environment.land_binary_mask == 1) & (self.elements.age_seconds < self.get_config('biology:min_settlement_age_seconds')))[0]
                on_land_and_older = np.where((self.environment.land_binary_mask == 1) & (self.elements.age_seconds >= self.get_config('biology:min_settlement_age_seconds')))[0]

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
                                         (self.elements.age_seconds >= self.get_config('biology:min_settlement_age_seconds')),
                                         reason='settled_on_coast')