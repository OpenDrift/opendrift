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
# 
##################################################################################################
# 
# Module to simulate behaviour of Plankton, including phytoplankton (plants) and zooplankton (animals)
#
# Developed by Simon Weppe based on several sources to fit requirements for use in MetOceanTrack model
# developped by MetOcean/MetService NZ as part of the Moana project (https://www.moanaproject.org/)
# 
##################################################################################################
# 
# Sources:
# 
# This module is based on several sources:
# a) the work and codes of Kristiansen,Romagnoni,Kvile
#   >> https://github.com/trondkr/KINO-ROMS/tree/master/Romagnoni-2019-OpenDrift/kino
#   >> https://github.com/trondkr/KINO-ROMS/blob/master/Romagnoni-2019-OpenDrift/kino/pelagicplankton.py
#
# "This code simulates cod egg and larvae and was developed by Trond Kristiansen (me (at) trondkristiansen.com)
# and Kristina Kvile (kristokv (at) gmail.com)""
# 
# b) Opendrift's pelagicegg module 
# >> https://github.com/OpenDrift/opendrift/blob/master/opendrift/models/pelagicegg.py
# 
# c ) Some components of code are taken from 
# >> https://github.com/metocean/ercore/tree/ercore_opensrc/ercore/lib
# 
##################################################################################################

import os
import numpy as np
from datetime import datetime, timedelta,timezone

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.elements import LagrangianArray

# Necessary to import the Fortran module that calculates light
if False: # disabling for now
    from kino import calclight
    from kino import bioenergetics

# Defining the PelagicPlankton element properties
class PelagicPlankton(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic plankton

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
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),
        ('hatched', {'dtype': np.int64,
                     'units': '',
                     'default': 0}),
        ('stage_fraction', {'dtype': np.float32,   #KK: Added to track percentage of development time completed
                         'units': '',
                         'default': 0.}),               
        ('length', {'dtype': np.float32,
                         'units': 'mm',
                         'default': 4.0}),
        ('weight', {'dtype': np.float32,
                         'units': 'mg',
                         'default': 0.08}),
        ('Eb', {'dtype': np.float32,
                     'units': 'ugEm2',
                     'default': 0.}),
         ('light', {'dtype': np.float32,
                     'units': 'ugEm2',
                     'default': 0.}),
         ('growth_rate', {'dtype': np.float32,   
                          'units': '%',
                          'default': 0.}),
         ('ingestion_rate', {'dtype': np.float32,   
                          'units': '%',
                          'default': 0.}),
         ('stomach_fullness', {'dtype': np.float32, 
                          'units': '%',
                          'default': 0.}),
         ('stomach', {'dtype': np.float32, 
                          'units': '',
                          'default': 0.}),
         ('survival', {'dtype': np.float32, # generic load to track mortality
                          'units': '',
                          'default': 1.})])

class PelagicPlanktonDrift(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

        Under construction.
    """

    ElementType = PelagicPlankton

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          #'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature',
                          'sea_water_salinity',
                          'surface_downward_x_stress',
                          'surface_downward_y_stress',
                          'turbulent_kinetic_energy',
                          'turbulent_generic_length_scale',
                          'upward_sea_water_velocity'
                          # 'cloud_coverage' ?
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
    required_profiles_z_range = [-150, 0]  # The depth range (in m) which
                                          # profiles shall cover

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       #'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 1000,
                       'ocean_vertical_diffusivity': 0.02,  # m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity': 34.,
                       'surface_downward_x_stress': 0,
                       'surface_downward_y_stress': 0,
                       'turbulent_kinetic_energy': 0,
                       'turbulent_generic_length_scale': 0,
                       'upward_sea_water_velocity': 0
                       }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(PelagicPlanktonDrift, self).__init__(*args, **kwargs)

        # By default, eggs do not strand towards coastline
        self.set_config('general:coastline_action', 'previous')

        # Vertical mixing is enabled by default
        self.set_config('drift:vertical_mixing', True)

        #IBM-specific configuration options
        # specifications on pelagicplankton module by Trond Kristiansen
        # https://github.com/trondkr/KINO-ROMS/blob/master/Romagnoni-2019-OpenDrift/kino/pelagicplankton.py
        self._add_config('biology:constantIngestion', 'float(min=0.0, max=1.0, default=0.5)', comment='Ingestion constant')
        self._add_config('biology:activemetabOn', 'float(min=0.0, max=1.0, default=1.0)', comment='Active metabolism')
        self._add_config('biology:attenuationCoefficient', 'float(min=0.0, max=1.0, default=0.18)', comment='Attenuation coefficient')
        self._add_config('biology:fractionOfTimestepSwimming', 'float(min=0.0, max=1.0, default=0.15)', comment='Fraction of timestep swimming')
        self._add_config('biology:lowerStomachLim', 'float(min=0.0, max=1.0, default=0.3)', comment='Limit of stomach fullness for larvae to go down if light increases')
        self._add_config('biology:haddock', 'boolean(default=False)', comment='Species=haddock')
        self._add_config('biology:cod', 'boolean(default=True)', comment='Species=cod')

         # ERcore specifications:
         # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L25 
        self._add_config('biology:mortality_daily_rate', 'float(min=0.0, max=100.0, default=0.05)', comment='Mortality rate (percentage of biomass dying per day)') 
        self._add_config('biology:min_settlement_age_seconds', 'float(min=0.0, max=7.0e6, default=0.0)', comment='Minimum age before beaching can occur, in seconds')
        self._add_config('biology:vertical_position_daytime', 'float(min=-1000.0, max=0.0, default=-5.0)',   comment='the depth a species is expected to inhabit during the day time, in meters, negative down') #
        self._add_config('biology:vertical_position_nighttime', 'float(min=-1000.0, max=0.0, default=-1.0)', comment='the depth a species is expected to inhabit during the night time, in meters, negative down') #
        # vertical_migration_speed_constant is the speed at which larvae will move to vertical_position_daytime or vertical_position_nighttime (depeding if it is day or nighttime)
        # if set to None the model will use the vertical velocity defined in update_terminal_velocity() (and so will not take into account vertical_position_daytime/vertical_position_nighttime)
        self._add_config('biology:vertical_migration_speed_constant', 'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity()') #
        self._add_config('biology:temperature_min', 'float(min=0.0, max=100.0, default=None)', comment=' lower threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
        self._add_config('biology:temperature_max', 'float(min=0.0, max=100.0, default=None)', comment=' upper threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
        self._add_config('biology:temperature_tolerance', 'float(min=0.0, max=1.0, default=1.0)', comment=' temperature tolerance before dying in degrees Celsius') #
        self._add_config('biology:salinity_min', 'float(min=0.0, max=100.0, default=None)', comment=' lower threshold salinity where a species population quickly declines to extinction in ppt') #
        self._add_config('biology:salinity_max', 'float(min=0.0, max=100.0, default=None)', comment=' upper threshold salinity where a species population quickly declines to extinction in ppt') #
        self._add_config('biology:salinity_tolerance', 'float(min=0.0, max=1.0, default=1.0)', comment=' salinity tolerance before dying in ppt') #
        # below not implemented yet
        self._add_config('biology:thermotaxis', 'float(min=0.0, max=1.0, default=None)', comment='movement of an organism towards or away from a source of heat, in (m/s per C/m)') #
        self._add_config('biology:halotaxis', 'float(min=0.0, max=1.0, default=None)', comment='movement of an organism towards or away from a source of salt, in (m/s per PSU/m)') #

    
    def update_terminal_velocity(self): 
        # function to update vertical velocities (negative = settling down)
        # called in update() method
        if self.get_config('biology:vertical_migration_speed_constant') is None:
            self.update_terminal_velocity_pelagicegg() # same as pelagicegg.py
        else:
            self.update_terminal_velocity_constant() # constant migration rate towards day or night time positions

    def update_terminal_velocity_pelagicegg(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for Pelagic Egg

            according to
            S. Sundby (1983): A one-dimensional model for the vertical
            distribution of pelagic fish eggs in the mixed layer
            Deep Sea Research (30) pp. 645-661

            Method copied from ibm.f90 module of LADIM:
            Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
            Fish. Oceanogr. (16) pp. 216-228

            same as Opendrift's PelagicEggDrift model. This function is called in update()
        """

        g = 9.81  # ms-2

        # Pelagic Egg properties that determine buoyancy
        eggsize = self.elements.diameter  # 0.0014 for NEA Cod
        eggsalinity = self.elements.neutral_buoyancy_salinity
        # 31.25 for NEA Cod

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

        # The density difference bettwen a pelagic egg and the ambient water
        # is regulated by their salinity difference through the
        # equation of state for sea water.
        # The Egg has the same temperature as the ambient water and its
        # salinity is regulated by osmosis through the egg shell.
        DENSw = self.sea_water_density(T=T0, S=S0)
        DENSegg = self.sea_water_density(T=T0, S=eggsalinity)
        dr = DENSw-DENSegg  # density difference

        # water viscosity
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0)
        # ~0.0014 kg m-1 s-1

        # terminal velocity for low Reynolds numbers
        W = (1.0/my_w)*(1.0/18.0)*g*eggsize**2 * dr

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W*1000*eggsize/my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * T0)  # in cm2/s
        d0 = (eggsize * 100) - 0.4 * \
            (9.0 * my_w**2 / (100 * g) * DENSw / dr)**(1.0 / 3.0)  # cm
        W2 = 19.0*d0*(0.001*dr)**(2.0/3.0)*(my_w*0.001*DENSw)**(-1.0/3.0)
        # cm/s
        W2 = W2/100.  # back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    ##################################################################################################
    # IBM-specific routines
    ##################################################################################################
    # >> reproduce behaviours included in Plankton class
    #    in https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py
    # 
    # 
    ##################################################################################################
    
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
        self.logger.debug('Assuming UTC time for solar calculations')
        # longitude convention in pysolar, consistent with Opendrift : negative reckoning west from prime meridian in Greenwich, England
        # the particle longitude should be converted to the convention [-180,180] if that is not the case
        sun_altitude = solar.get_altitude(self.elements.lon, self.elements.lat, date) # get sun altitude in degrees
        sun_azimut = solar.get_azimuth(self.elements.lon, self.elements.lat, date) # get sun azimuth in degrees
        sun_radiation = np.zeros(len(sun_azimut))
        # not ideal get_radiation_direct doesnt accept arrays...
        for elem_i,alt in enumerate(sun_altitude):
            sun_radiation[elem_i] = solar.radiation.get_radiation_direct(date, alt)  # watts per square meter [W/m2] for that time of day
        self.elements.light = sun_radiation * 4.6 #Converted from W/m2 to umol/m2/s-1"" - 1 W/m2 ≈ 4.6 μmole.m2/s
        self.logger.debug('Solar radiation from %s to %s [W/m2]' % (sun_radiation.min(), sun_radiation.max() ) )
        # print(np.min(sun_radiation))
        # print(date)

        # Opendrift has a built-in function to work out solar elevation at current particle positions and date (from physics_methods.py)
        # but it returns results that are quite different from Pysolar.
        # print(sun_altitude)
        # print(self.solar_elevation()) 

        # Quick test on night time vs daytime
        # lon = -3 #(France)
        # lat = 47
        # date=  datetime(2020, 10, 6, 12, 0)
        # date = date.replace(tzinfo=timezone.utc)
        # date_vec = date + np.arange(48) * timedelta(hours=1)
        # sun_radiation = np.zeros(len(date_vec)) 
        # sun_radiation = np.zeros(len(date_vec))       
        # for ii,date_i in enumerate(date_vec):
        #     sun_altitude[ii] = solar.get_altitude(lon,lat, date_i)
        #     sun_radiation[ii] = solar.radiation.get_radiation_direct(date_i, sun_altitude[ii])
        # import matplotlib.pyplot as plt
        # plt.ion()
        # plt.plot(date_vec,sun_radiation)
        # import pdb;pdb.set_trace()
    
    def plankton_development(self):

        self.update_survival() # mortality based on user-input mortality rate
        self.update_weight_temperature() # # mortality based on "liveable" temperature range
        self.update_weight_salinity() # # mortality based on "liveable" salinity range

    def update_survival(self):
        # update survval fraction based on mortality rate in [day-1]
        mortality_daily_rate = self.get_config('biology:mortality_daily_rate')
        # update survival fraction, accounting for mortality rate over that timestep
        fraction_died = self.elements.survival * mortality_daily_rate * (self.time_step.total_seconds()/(3600*24))
        self.elements.survival -=  fraction_died
        # import pdb;pdb.set_trace()

    def update_weight_temperature(self):
        #update particle survival fraction based on temperature, if a temperature range was input
        temp_tol = self.get_config('biology:temperature_tolerance')
        temp_min = self.get_config('biology:temperature_min')
        temp_max = self.get_config('biology:temperature_max')
        temp_xy = self.environment.sea_water_temperature # at particle positions

        if (temp_min is not None) or (temp_max is not None) :
            if temp_max is not None :
                m=(temp_xy-temp_max+temp_tol)/temp_tol # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L60
                if (m>0).any():
                  self.logger.debug('Maximum temperature reached for %s particles' % np.sum(m>0))
                self.elements.survival -= np.maximum(np.minimum(m,1),0)*self.elements.survival # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L62
            if temp_min is not None :
                m=(temp_min+temp_tol-temp_xy)/temp_tol # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L64
                if (m>0).any():
                  self.logger.debug('Minimum temperature reached for %s particles' % np.sum(m>0))
                self.elements.survival -= np.maximum(np.minimum(m,1),0)*self.elements.survival # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L65
        # print('TEMP')
        # print(self.elements.survival)
        # import pdb;pdb.set_trace()

    def update_weight_salinity(self):
        #update particle survival fraction based on salinity, if a salinity range was input
        salt_tol = self.get_config('biology:salinity_tolerance')
        salt_min = self.get_config('biology:salinity_min')
        salt_max = self.get_config('biology:salinity_max')
        salt_xy = self.environment.sea_water_salinity # at particle positions

        if (salt_min is not None) or (salt_max is not None) :
            if salt_max is not None :
                m=(salt_xy-salt_max+salt_tol)/salt_tol # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L60
                if (m>0).any():
                  self.logger.debug('Maximum salinity reached for %s particles' % np.sum(m>0))
                self.elements.survival -= np.maximum(np.minimum(m,1),0)*self.elements.survival # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L73
            if salt_min is not None :
                m=(salt_min+salt_tol-salt_xy)/salt_tol # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L64
                if (m>0).any():
                  self.logger.debug('Minimum salinity reached for %s particles' % np.sum(m>0))
                self.elements.survival -= np.maximum(np.minimum(m,1),0)*self.elements.survival # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L77
        # print('SALT')
        # print(self.elements.survival)
        # import pdb;pdb.set_trace()

    def update_terminal_velocity_constant(self):
        # modifies the same variable than update_terminal_velocity(), self.elements.terminal_velocity = W, but using a different algorithm.
        # Larvae are assumed to move to daytime or nighttime vertical positions in the water column, at a constant rate
        #
        # the actual settling is taken care of in vertical_mixing() or vertical_buoyancy() (i.e. from OceanDrift methods)
        self.calculateMaxSunLight() # compute solar radiation at particle positions (using PySolar)
        # it is expected that larve will go down during day time and up during night time but that is not fixed in the code. 
        # Particles will simply go towards the daytime or nighttime poistions.
        # https://github.com/metocean/ercore/blob/ercore_nc/ercore/materials/biota.py#L80
        vertical_velocity = np.abs(self.get_config('biology:vertical_migration_speed_constant'))  # magnitude in m/s 
        z_day = self.get_config('biology:vertical_position_daytime')    #  the depth a species is expected to inhabit during the day time, in meters, negative down') #
        z_night =self.get_config('biology:vertical_position_nighttime') # 'the depth a species is expected to inhabit during the night time, in meters, negative down') #
        ind_day = np.where(self.elements.light>0)
        ind_night = np.where(self.elements.light==0)
        self.logger.debug('Using constant migration rate (%s m/s) towards day and night time positions' % (vertical_velocity) )
        self.logger.debug('%s particles in day time' % (len(ind_day[0])))
        self.logger.debug('%s particles in night time' % (len(ind_night[0])))
        # for particles in daytime : particles below the daytime position need to go up while particles above the daytime position need to go down
        # (same for for particles in nightime)
        # Note : depth convention is neagtive down in Opendrift
        # 
        # e.g. z=-5, z_day = -3, below daytime position,  need to go up (terminal_velocity>0) 
        #      diff = (z - z_day) = -2, so w = - np.sign(diff) * vertical_velocity
        self.elements.terminal_velocity[ind_day] = - np.sign(self.elements.z[ind_day] - z_day) * vertical_velocity
        self.elements.terminal_velocity[ind_night] = - np.sign(self.elements.z[ind_night] - z_night) * vertical_velocity
        # print(self.elements.z)

    def spawn(self):
        pass
        # particles change types after a certain time
        # need more info to implement
        # this could call a new seed_elements() method (that may or may not need to be changed relative to the one in basemodel.py) or 
        # simply change some parameters of existing particles at a certan time...
        # 

    ##################################################################################################
    # Not used for now
    ##################################################################################################
    def calculateMaximumDailyLight(self):
        """LIGHT == Get the maximum light at the current latitude, and the current surface light at the current time.
        These values are used in calculateGrowth and  predation.FishPredAndStarvation, but need only be calcuated once per time step. """

        tt = self.time.timetuple()
      
        num=np.shape(self.elements.lat)[0]
        dayOfYear = float(tt.tm_yday); month=float(tt.tm_mon); hourOfDay = float(tt.tm_hour); daysInYear=365.0
        radfl0=np.zeros((num)); maxLight=np.zeros((num)); 
        cawdir=np.zeros((num)); clouds=0.0; 
        sunHeight=np.zeros((num)); surfaceLight=np.zeros((num))

        """Calculate the maximum and average light values for a given geographic position
        for a given time of year. Notice we use radfl0 instead of maximum values maxLight
        in light caclualtions. Seemed better to use average values than extreme values.
        NOTE: Convert from W/m2 to umol/m2/s-1"""
        radfl0,maxLight,cawdir = calclight.calclight.qsw(radfl0,maxLight,cawdir,clouds,self.elements.lat*np.pi/180.0,dayOfYear,daysInYear,num)
        maxLight = radfl0/0.217
        """Calculate the daily variation of irradiance based on hour of day - store in self.elements.light"""
        for ind in range(len(self.elements.lat)):
        
          sunHeight, surfaceLight = calclight.calclight.surlig(hourOfDay,float(maxLight[ind]),dayOfYear,float(self.elements.lat[ind]),sunHeight,surfaceLight)
          self.elements.light[ind]=surfaceLight

    def updateSurvival(self):
        # Update the size dependent mortality
        k = 0.06
        x = 0.4
        mortality_scale = 1 # To increase or decrease the mortality rate

        for ind in range(len(self.elements.lat)):
            if (self.elements.hatched[ind]>=1):
                #Larval mortality (per day):
                 mortality = k*(self.elements.weight[ind]**-x)*(self.time_step.total_seconds()/(24*3600.))
                 mortality = mortality_scale*mortality
            else:
                #Egg mortality:
                mortality = 0.2*(self.time_step.total_seconds()/(24*3600.))
                mortality = mortality_scale*mortality
                self.elements.survival[ind] = self.elements.survival[ind]*(np.exp(-mortality))

    def updateEggDevelopment(self):
        # Update percentage of egg stage completed
        amb_duration = np.exp(3.65 - 0.145*self.environment.sea_water_temperature) #Total egg development time (days) according to ambient temperature (Ellertsen et al. 1988)
        days_in_timestep = self.time_step.total_seconds()/(60*60*24)  #The fraction of a day completed in one time step
        amb_fraction = days_in_timestep/(amb_duration) #Fraction of development time completed during present time step 
        self.elements.stage_fraction += amb_fraction #Add fraction completed during present timestep to cumulative fraction completed
        self.elements.hatched[self.elements.stage_fraction>=1] = 1 #Eggs with total development time completed are hatched (1)

    def updateVerticalPosition(self,length,oldLight,currentLight,currentDepth,stomach_fullness,dt):
        # Update the vertical position of the current larva
        swimSpeed=0.261*(length**(1.552*length**(0.920-1.0)))-(5.289/length)
        fractionOfTimestepSwimming = self.get_config('biology:fractionOfTimestepSwimming')
        maxHourlyMove = swimSpeed*fractionOfTimestepSwimming*dt
        maxHourlyMove =  round(maxHourlyMove/1000.,1) # depth values are negative
        lowerStomachLim = self.get_config('biology:lowerStomachLim')

        if (oldLight <= currentLight and stomach_fullness >= lowerStomachLim): #If light increases and stomach is sufficiently full, go down
          depth = min(0.0,currentDepth - maxHourlyMove)
        else: #If light decreases or stomach is not sufficiently full, go up
          depth = min(0,currentDepth + maxHourlyMove)
        #print "current depth %s new depth %s light %s lastLight %s"%(currentDepth,depth,currentLight,oldLight)
        #print "maximum mm to move %s new depth %s old depth %s"%(maxHourlyMove,depth,currentDepth)
        return depth

    def updatePlanktonDevelopment(self):
   
        constantIngestion = self.get_config('biology:constantIngestion')
        activemetabOn = self.get_config('biology:activemetabOn')
        attCoeff = self.get_config('biology:attenuationCoefficient')
        haddock = self.get_config('biology:haddock')
        fractionOfTimestepSwimming = self.get_config('biology:fractionOfTimestepSwimming')

        self.updateEggDevelopment()
        
        # Save the light from previous timestep to use for vertical behavior
        lastLight=np.zeros(np.shape(self.elements.light))
        lastLight[:]=self.elements.light[:]

        self.calculateMaximumDailyLight()

        self.elements.Eb=self.elements.light*np.exp(attCoeff*(self.elements.z))
     
        dt=self.time_step.total_seconds()

        for ind in range(len(self.elements.lat)):
          if (self.elements.hatched[ind]>=1):
          # Calculate biological properties of the individual larva                                 
            larvamm,larvawgt,stomach,ingrate,growthrate,sfullness=bioenergetics.bioenergetics.growth(self.elements.length[ind],
              self.elements.weight[ind],
              self.elements.stomach[ind],
              self.elements.ingestion_rate[ind],
              self.elements.growth_rate[ind],
              self.elements.stomach_fullness[ind],
              haddock,
              activemetabOn,
              constantIngestion,
              fractionOfTimestepSwimming,
              self.elements.Eb[ind],
              dt,
              self.environment.sea_water_temperature[ind])
            # Update the element

            self.elements.length[ind]=larvamm
            self.elements.weight[ind]=larvawgt
            self.elements.stomach[ind]=stomach
            self.elements.ingestion_rate[ind]=ingrate
            self.elements.growth_rate[ind]=growthrate
            self.elements.stomach_fullness[ind]=sfullness

            self.elements.z[ind] = self.updateVertialPosition(self.elements.length[ind],
              lastLight[ind],
              self.elements.light[ind],
              self.elements.z[ind],
              self.elements.stomach_fullness[ind],
              dt)
    ##################################################################################################
    # above not used for now
    ##################################################################################################

    def update(self):
        """Update positions and properties of buoyant particles."""

        # Update element age
        self.elements.age_seconds += self.time_step.total_seconds()

        # move particles with ambient current
        self.advect_ocean_current()

        # Disabled for now
        if False:
            # Advect particles due to surface wind drag,
            # according to element property wind_drift_factor
            self.advect_wind()

            # Stokes drift
            self.stokes_drift()

        # Turbulent Mixing or settling-only 
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()  #compute vertical velocities, two cases possible - constant, or same as pelagic egg
            self.vertical_mixing()
        else:  # Buoyancy
            self.update_terminal_velocity()
            self.vertical_buoyancy()

        # Vertical advection
        self.vertical_advection()
        
        # Plankton specific
        if True:
          self.plankton_development()
        
        if False:
          # Plankton development
          self.updatePlanktonDevelopment()
          self.updateSurvival()