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
# Copyright 2021, Trond Kristiansen, Niva
# Jan 2021 Simplified by Knut-Frode Dagestad, MET Norway, and adapted to to Kvile et al. (2018)

import datetime
import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift


class SeaLiceElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for larval and
    juvenile stages of sea lice into super individuals
    """

    variables = Lagrangian3DArray.add_variables([
        ('particle_biomass',{'dtype': np.float,
                            'units': 'kg',
                            'default': 1000.}),
        ('hatched',{'dtype': np.float,
                            'units': '',
                            'default': 0.}),
        ('nauplii', {'dtype': np.float,
                            'units': '',
                            'default': 0.}),
        ('copepodid', {'dtype': np.float,
                     'units': '',
                     'default': 0.}),
        ('dead', {'dtype': np.bool,
                     'units': '',
                     'default': False}),
        ('degree_days', {'dtype': np.float,
                     'units': '',
                    'default': 0}), #range 40-170
        ('salinity_above', {'dtype': np.bool,
                     'units': '',
                    'default': False}),
        ('temperature_above', {'dtype': np.float,
                     'units': '',
                    'default': 10}),
        ('temperature_below', {'dtype': np.float,
                     'units': '',
                    'default': 10}),
        ('light', {'dtype': np.float,
                     'units': 'µmol photon  s−1 m−2',
                     'default': 0.})
                      ])


class SeaLice(OceanDrift):
    """
    Particle trajectory model based on the OpenDrift framework.
    Developed by Julien Moreau (Plastics@Bay CIC)
    Generic module for particles that are subject to vertical turbulent
    mixing with the possibility for positive or negative buoyancy
    Particles are sea-lice (Lepeophtheirus salmonis).
    """

    ElementType = SeaLiceElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 50},
        'surface_net_downward_radiative_flux':{'fallback': None},
        'ocean_vertical_diffusivity': {'fallback': 0.01, 'profiles': True},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True}
    }

    # required_profiles_z_range = [0, -50]  # The depth range (in m) which profiles should cover

    def __init__(self,*args, **kwargs):

        # Calling general constructor of parent class
        super(SeaLice, self).__init__(*args, **kwargs)

        # Configuration options
        self._add_config({
            'general:LicePerFish':{'type':np.float, 'default':0.5, # this could become an element property
                                'min': 0.0, 'max': None, 'units': '',
                                'description': 'Number of lice per fish',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:AvFishW8':{'type':np.float, 'default':4.5, # this could become an element property
                                'min': 0.0, 'max': 1.0, 'units': 'fraction',
                                'description': 'Average fish weight',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:seeding_time_step':{'type':np.float, 'default':None,
                                'min': None, 'max': None, 'units': 'seconds',
                                'description': 'Time between particle release',
                                'level': self.CONFIG_LEVEL_ESSENTIAL},
            'general:death_rate':{'type':np.float, 'default':0.01/3600,
                                'min': 0., 'max': None, 'units': 's-1',
                                'description': 'Rate of Larvae death per seconds',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:maturation_rate':{'type':np.float, 'default':0.1/3600,
                                'min': 0., 'max': None, 'units': 's-1',
                                'description': 'Rate of Nauplii maturation in Copepodids',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:maturity_date':{'type':np.float, 'default':3.63,
                                'min': 0., 'max': None, 'units': 'days',
                                'description': 'Days to start maturing into Copepodids',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:sinking_velocity':{'type':np.float, 'default':0.00025,
                                'min': 0., 'max': 0.01, 'units': 'm.s-1',
                                'description': 'Larvae sinking velocity',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:vertical_migration_speed':{'type':np.float, 'default':0.00075,
                                'min': 0., 'max': 0.01, 'units': 'm.s-1',
                                'description': 'Larvae vertical speed',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:freezing_salinity':{'type':np.float, 'default':27,
                                'min': 0., 'max': 35., 'units': 'PSU',
                                'description': 'Salinity immobilising larvae',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:avoided_salinity':{'type':np.float, 'default':32,
                                'min': 0., 'max': 50., 'units': 'PSU',
                                'description': 'Salinity actively avoided',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:nu':{'type':np.float, 'default':500.,
                                'min': 100, 'max': 1000, 'units': 'nm',
                                'description': 'Wavelength used to calculate irradiance',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'general:k_water':{'type':np.float, 'default':-0.2,
                                'min': -10., 'max': 0., 'units': '',
                                'description': 'coefficient of exponential decay of light in water',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            'general:Nauplii_light_trigger':{'type':np.float, 'default':2.E-5,
                                'min': 0., 'max': 1., 'units': 'µmol photon  s−1 m−2',
                                'description': 'light detection threshold of Nauplii',
                                'level': self.CONFIG_LEVEL_BASIC},
            'general:Copepodid_light_trigger':{'type':np.float, 'default':0.392,
                                'min': 0.0, 'max': 1.0, 'units': 'µmol photon  s−1 m−2',
                                'description': 'light detection threshold of copepodids',
                                'level': self.CONFIG_LEVEL_BASIC},
             'general:twilight':{'type':np.float, 'default':15,
                                'min': 0., 'max': 90., 'units': 'degrees',
                                'description': 'angle below the horizon for twilight',
                                'level': self.CONFIG_LEVEL_ADVANCED},
            })



        # self.LicePerFish=parameters['LicePerFish'] # need to be added in module parameters
        # self.AvFishW8 = parameters['AvFishW8'] # need to be added in module parameters
        # # self.licing_time_step= parameters['licing_time_step']# need to be added in module parameters kee it in seconds
        # self.particle_biomass= parameters['particle_biomass'] #number of particle per farm per emission per hour
        # self.death_rate = parameters['Killing_rate'] #kill per second
        # self.Mp = parameters['Maturation_rate'] #maturity per second
        # # change that into degree days in next generation
        #
        # self.sinking_velocity = parameters['sinking_velocity'] # m.s-1
        # self.vertical_migration_speed=parameters['vertical_migration_speed'] #m.s-1
        #
        # self.freezing_salinity = parameters['freezing_salinity'] #PSU
        # self.avoided_salinity=parameters['avoided_salinity'] # PSU
        # self.nu=parameters['nu'] #mainwavelength of irradiance in nanometres
        # self.k_water=parameters['k_water'] # parameter decay of light in m-1
        # self.Nauplii_light_trigger= parameters['Nauplii_light_trigger']# µmol photon  s−1 m−2
        # self.Copepodid_light_trigger=parameters['Copepodid_light_trigger'] # µmol photon  s−1 m−2
        # self.twilight = np.deg2rad(parameters['twilight']) #Nautical horizon


    def prepare_run(self):
        self.sensing_distance=2*self.vertical_migration_speed*self.time_step.total_seconds()
        self.Mat = int(np.ceil(self.maturity_date*24*3600/self.time_step.total_seconds())) # maturity age in timestep
        self.tau= np.pi+2*self.twilight #width of solar irridation distribution
        self.central_lon, self.central_lat=self.xy2lonlat((self.xmax-self.xmin)/2,
        (self.ymax-self.ymin)/2)
        logger.info("center of area for solar angle is lat: {}, lon: {}"\
        .format(self.central_lat,self.central_lon))
        logger.info("Making super individual population evolution...")
        t=np.arange(0,(self.end_time-self.start_time).total_seconds(),\
            self.time_step.total_seconds())
        self.juv, self.adult, self.dead=population(t, self.Mat, self.maturation_rate, self.death_rate)
        logger.info("Loading NASA irradiance data...")
        from opendrift.readers.reader_netCDF_CF_generic import Reader
        r = Reader('https://opendap.larc.nasa.gov/opendap/hyrax/POWER/daily/\
        power_801_daily_allsky_sfc_sw_dwn_lst.nc', standard_name_mapping=\
        {'ALLSKY_SFC_SW_DWN': 'surface_net_downward_radiative_flux'})
        self.add_reader(r)
        logger.info("NASA irradiance data reader ready.")

    def new_born(self, Particle_biomass):
        """
        Approach by Rittenhouse et al., (2016) [Rittenhouse, M.A., C.W. Revie
        and A. Hurford, 2016. A model for sea lice (Lepeophtheirus salmonis)
        dynamics in a seasonally changing environment. Energetics, 16, 8–16.,
        whereby the numbers of lice nauplii released every day (NP) are
        estimated according to the following equation:
        NP(t)=ηενA(t)
        where η = 592 is the number of eggs per clutch (Heuch et al., 2000;
        Rittenhouse et al., 2016),
        ε = 0.0476 is the egg string production rate per day (Heuch et al.,
        2000; Rittenhouse et al., 2016),
        ν = 0.6 is the hatching success i.e. the proportion of eggs which
        produce viable nauplii (Johnson and Albright, 1991),
        and A(t) is the total number of adult female lice on each farm,
        derived here from lice target levels and number of fish on site.
        """
        Ne, nu = 592 , 0.6
        Eps=0.0476/timedelta(days=1).total_seconds()*self.seeding_time_step
        return Ne* Eps* nu*Particle_biomass*self.LicePerFish/self.AvFishW8

    def population(t, Mat, Mp, Kp):
        """
        Simulate the population evolution in each farm according to its Biomass
        Parameters
        ----------
        t: numpy array (int)
        Time abscissa in timesteps
        Mat: int
        Age in timestep when the nauplii start to become adult copepodids
        Mp: float
        Probability of maturation by time-step (t>= Mat)
        Kp: float
        Probability of death by time-step
        Lice_parts: float
        Number of lice per particle when generated
        Exp: int
        Number of timesteps in the experiment
        Returns
        -------
        juv, adults, dead: np.arrays(int)
        Juveniles(Nauplii), adults(Copepodid), dead lices arrays through
        the duration of the experiment
        """
        juv=np.exp(-Kp*t)
        dead=Lice_parts-juv
        adult=np.zeros_like(t)
        decayjuv=juv[Mat-1]*np.exp(-(Mp+Kp)*(t[t>=Mat]-Mat))
        juv[t>=Mat] =decayjuv
        adult[t>=Mat]-=dead[t>=Mat]-juv[t>=Mat]
        return juv, adult, dead

    # def super_indiv_aging(self):
    #     """
    #     Create a dictionary of population evolution to read through the model
    #     """
    #     t = np.arange(0,(self.end_time-self.start_time).total_seconds(),\
    #         self.time_step.total_seconds())
    #     self.Juvs, self.Adults, self.Deads = np.zeros((len(self.particle_biomass), \
    #         len(t))), np.zeros((len(self.particle_biomass), len(t))), \
    #         np.zeros((len(self.particle_biomass), len(t)))
    #     self.desactivation_time=np.zeros(len(self.particle_biomass))
    #     for k in range(len(self.particle_biomass)):
    #         N0 = new_born(self, self.particle_biomass[k])
    #         self.Juvs[k], self.Adults[k], self.Deads[k] = population(t, self.Mat,
    #                                                     self.Mp, self.death_rate,N0)
    #         self.desactivation_time[k]= t[np.where(Deads[k]>N0-1)[0][0]] # all deads

    def SI_pop(self):
        """
        distribute the age fractions in the particles
        """
        # First desactivate dead particles
        dead=self.elements.dead
        alive= np.logical_not(dead)
        Dying=self.elements.hatched[alive]<self.elements.nauplii[alive]+self.elements.copepodid[alive]+1
        self.elements.dead[Dying]=True
        self.deactivate_elements(self.elements.dead, reason="All dead")
        # Take care of new borns
        alive= np.logical_not(dead)
        New_release = self.elements.age_seconds[alive]<self.time_step
        self.elements.hatched[New_release]=new_born(self, self.elements.Particle_biomass[New_release])
        # Make others mature
        Free = np.logical_not(New_release & dead)
        self.elements.nauplii[Free]=self.elements.hatched[Free]*self.juv[self.elements.age_seconds[Free]]
        self.elements.copepodid[Free]=self.elements.hatched[Free]*self.adult[self.elements.age_seconds[Free]]

    def sensing(self, Normal_salt):
        """
        Lice sensing if above or bellow the conditions are better for them
        within a specified distance.
        moving the particles up and then down then back in their initial position
        """
        Backup= self.elements.z[Normal_salt]
        self.elements.z[Normal_salt] +=self.sensing_distance
        self.elements.salinity_above[Normal_salt]=self.environment.sea_water_salinity>self.avoided_salinity
        self.elements.temperature_above[Normal_salt]= self.environment.sea_water_temperature
        self.elements.z[Normal_salt] -= 2* self.sensing_distance
        self.elements.temperature_below[Normal_salt]= self.environment.sea_water_temperature
        self.elements.z[Normal_salt] = Backup

    def degree_days(self):
        """
        Calculate the degree days of a particles
        >>> under development <<<
        """
        self.elements.degree_days+=self.environment.sea_water_temperature* \
                self.time_step.total_seconds()/timedelta(days=1).total_seconds()

    def solar_energy_distribution(self):
        """
        Distribute the daily energy from irradiance with a gaussian distribution.
        We use the twilight
        """
        self.solar_angle=np.deg2rad(self.solar_elevation(self.time,
                                    self.central_lon, self.central_lat))
        self.solar_coeff= np.sqrt(self.tau/(2*np.pi)*exp(-self.tau/2*self.solar_angle**2))

    def irradiance(self,Irradiance_watt, z):
        """
        Convert irradiance from W.m-2 to micromol photon.s-1.m-2
        https://www.berthold.com/en/bioanalytic/knowledge/faq/irradiance-to-photon-flux/
        calculate the photon flux in the water according to exponential decay from the sea surface
        """
        Irr0=self.solar_coeff*Irradiance_watt*self.nu*0.00836
        return Irr0*np.exp(self.k_water*z)

    def Lice_vertical_migration(self):
        """
        Make larval sea lice to migrate vertically according to salinity, light
        and temperature triggers.
        """
        ### all lice sink at the same speed
        self.elements.z -= self.sinking_velocity*self.time_step.total_seconds()

        ### filter actively avoiding salt
        Frozen=self.environment.sea_water_salinity<self.freezing_salinity
        Avoiding = np.logical_and((self.freezing_salinity< self.environment.sea_water_salinity),
                (self.environment.sea_water_salinity<self.avoided_salinity))
        Normal_salt = self.environment.sea_water_salinity>=self.avoided_salinity

        self.sensing(self, Normal_salt)
        ### Irradiance computation
        self.elements.light = irradiance(self, self.environment.surface_net_downward_radiative_flux,
                            self.elements.z)
        ### identify the elements involved in the different scenarios
        Filter_N= self.elements.copepodid < self.elements.nauplii
        Filter_C=np.logical_not(Filter_N)
        light_mig_N= self.environment.light[Filter_N&self.elements.salinity_above[Normal_salt]] > self.Nauplii_light_trigger
        light_mig_C= self.environment.light[Filter_C&self.elements.salinity_above[Normal_salt]] > self.Copepodid_light_trigger
        up_temp_mig= self.elements.warmer_above[np.logical_not(light_mig_N&light_mig_C) & self.elements.salinity_above[Normal_salt]] > self.environment.sea_water_temperature
        down_temp_mig=self.elements.warmer_below[np.logical_not(light_mig_N & light_mig_C & up_temp_mig)]>self.environment.sea_water_temperature

        ### Active migration
        self.elements.z[Avoiding&down_temp_mig] -= self.vertical_migration_speed*self.time_step.total_seconds()
        self.elements.z[light_mig_N&light_mig_C&up_temp_mig] += self.vertical_migration_speed*self.time_step.total_seconds()
        self.elements.z=depth_test(self.elements.z)

    def update(self):
        self.solar_energy_distribution()
        self.advect_ocean_current()
        self.SI_pop()
        self.Lice_vertical_migration()
