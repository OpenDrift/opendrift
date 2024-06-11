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
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# Copyright 2021, Gaute Hope, MET Norway
"""
OpenOil is a 3D oil drift module bundled within the OpenDrift framework.

The oil weathering calculations is based on the NOAA-ERR-ERD OilLibrary
package, which is installed as a dependency. The code for evaporation and
emulsification in OpenOil is borrowed from the NOAA PyGnome code, and adapted
to the OpenDrift architecture.

Example of ship leaking oil along the coast of Northern Norway
##############################################################

.. image:: https://dl.dropboxusercontent.com/s/ty6dmqf0oohewky/oilspill_tromsoe.gif?dl=0

Simulation of Deepwater Horizon (Macondo) accident, initiated from satellite images
###################################################################################

.. image:: https://dl.dropboxusercontent.com/s/ghi7crtmwpyjgto/macondo_simulation.gif?dl=0

Satellite images provided by Prof. Chuanmin Hu, and ocean model output provided by Prof. Vassiliki Kourafalou

Example oil budget for a simulation
###################################

.. image:: https://dl.dropboxusercontent.com/s/pb0h6tlev9pnoh3/oil_budget_draugen.png?dl=0

Oil properties affecting the drift
***********************************
The vertical (and thus indirectly also the horizontal) motion of oil (droplets) is affected by oil density and droplet diameters.

When using the NOAA oil weathering model (``o = OpenOil(weathering_model='noaa')``), which is the default, the density is obtained from the NOAA database according to the oiltype selected when seeding. This value can not be overridden by the user, and it will also change during the simulation due to oil weathering processes (evaporation and emulsification).

The droplet diameter may be given explicitly when seeding, e.g.:

.. testcode::


    o = OpenOil()
    o.set_config('environment:constant:x_wind', 0)
    o.set_config('environment:constant:y_wind', 0)
    o.set_config('environment:constant:x_sea_water_velocity', 0)
    o.set_config('environment:constant:y_sea_water_velocity', 0)
    o.seed_elements(4, 60, number=100, time=datetime.now(), diameter=1e-5)

In this case, the diameter will not change during the simulation, which is useful e.g. for sensitivity tests. The same diameter will be used for all elements for this example, but an array of the same length as the number of elements may also be provided.

If a constant droplet diameter is not given by the user, it will be chosen randomly within given config limits for a subsea spill ('blowout'), and modified after any later wave breaking event. Oil droplets seeded under sea surface (z<0) will be assigned initial diameters between the following limits, typical for a subsea blowout (Johansen, 2000)::

.. code::

    o.set_config('seed:droplet_diameter_min_subsea', 0.0005)  # 0.5 mm
    o.set_config('seed:droplet_diameter_max_subsea', 0.005)   # 5 mm

Alternatively, the user can specify normal or lognormal initial subsea droplet size distributions, which are later modified by wave breaking events. In these cases the user must specify the mean and standard deviation of the distribution::

.. code::

    o.set_config('seed:droplet_size_distribution','lognormal')
    o.set_config('seed:droplet_diameter_mu',0.001)  # 1 mm
    o.set_config('seed:droplet_diameter_sigma',0.0008) # 0.8 mm

Note that these config settings must be adjusted before the seeding call.
After each wave breaking event, a new droplet diameter will be chosen based on the config setting for droplet size distribution.
"""

from io import open
import numpy as np
from datetime import datetime
import pyproj
import matplotlib.pyplot as plt
import logging
import json
from importlib import resources

logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from . import noaa_oil_weathering as noaa
from . import adios
from adios_db.computation.physical_properties import KinematicViscosity, Density
from adios_db.computation import gnome_oil

from opendrift.models.physics_methods import oil_wave_entrainment_rate_li2017
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


# Defining the oil element properties
class Oil(Lagrangian3DArray):
    """Extending LagrangianArray with variables relevant for oil particles."""

    variables = Lagrangian3DArray.add_variables([
        ('mass_oil', {
            'dtype': np.float32,
            'units': 'kg',
            'seed': False,
            'default': 1
        }),
        (
            'viscosity',
            {
                'dtype': np.float32,
                'units': 'm2/s',
                'seed': False,  # Taken from NOAA database
                'description': 'Kinematic viscosity of oil (emulsion)',
                'default': 0.005
            }),
        (
            'density',
            {
                'dtype': np.float32,
                'units': 'kg/m^3',
                'seed': False,  # Taken from NOAA database
                'default': 880
            }),
        (
            'wind_drift_factor',
            {
                'dtype':
                np.float32,  # TODO: inherit from
                'units':
                '%',  # OceanDrift
                'description':
                'Elements at the ocean surface are moved by '
                'this fraction of the wind vector, in addition to '
                'currents and Stokes drift',
                'default':
                0.03
            }),
        ('bulltime', {
            'dtype': np.float32,
            'units': 's',
            'seed': False,
            'default': 0
        }),
        ('interfacial_area', {
            'dtype': np.float32,
            'units': 'm2',
            'seed': False,
            'default': 0
        }),
        ('mass_dispersed', {
            'dtype': np.float32,
            'units': 'kg',
            'seed': False,
            'default': 0
        }),
        ('mass_evaporated', {
            'dtype': np.float32,
            'units': 'kg',
            'seed': False,
            'default': 0
        }),
        ('mass_biodegraded', {
            'dtype': np.float32,
            'units': 'kg',
            'seed': False,
            'default': 0
        }),
        ('biodegradation_half_time_droplet',  # Note: if unit is "days", Xarray tries to decode as datetime with error
                {'dtype': np.float32, 'units': 'Days', 'seed': False, 'default': 1,
                 'description': 'Biodegradation half time in days for submerged oil droplets'}),
        ('biodegradation_half_time_slick',
                {'dtype': np.float32, 'units': 'Days', 'seed': False, 'default': 3,
                 'description': 'Biodegradation half time in days for surface oil slick'}),
        ('fraction_evaporated', {
            'dtype': np.float32,
            'units': '%',  # TODO: should be fraction and not percent
            'seed': False,
            'default': 0
        }),
        ('water_fraction', {
            'dtype': np.float32,
            'units': '%',
            'seed': False,
            'default': 0
        }),
        ('oil_film_thickness', {
            'dtype': np.float32,
            'units': 'm',
            'default': 0.001
        }),
        (
            'diameter',
            {
                'dtype': np.float32,  # Particle diameter
                'units': 'm',
                'seed': False,
                'default': 0.
            })
    ])


class OpenOil(OceanDrift):
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil

    required_variables = {
        'x_sea_water_velocity': {
            'fallback': None
        },
        'y_sea_water_velocity': {
            'fallback': None
        },
        'x_wind': {
            'fallback': None
        },
        'y_wind': {
            'fallback': None
        },
        'sea_surface_height': {'fallback': 0},
        'upward_sea_water_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_significant_height': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_stokes_drift_x_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_stokes_drift_y_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_period_at_variance_spectral_density_maximum': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment':
        {
            'fallback': 0,
            'important': False
        },
        'sea_ice_area_fraction': {
            'fallback': 0,
            'important': False
        },
        'sea_ice_x_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_ice_y_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_water_temperature': {
            'fallback': 10,
            'profiles': True
        },
        'sea_water_salinity': {
            'fallback': 34,
            'profiles': True
        },
        'sea_floor_depth_below_sea_level': {
            'fallback': 10000
        },
        'ocean_vertical_diffusivity': {
            'fallback': 0.02,
            'important': False,
            'profiles': True
        },
        'land_binary_mask': {
            'fallback': None
        },
        'ocean_mixed_layer_thickness': {
            'fallback': 50,
            'important': False
        },
    }


    # Default colors for plotting
    status_colors = {
        'initial': 'green',
        'active': 'blue',
        'missing_data': 'gray',
        'stranded': 'red',
        'evaporated': 'yellow',
        'dispersed': 'magenta'
    }

    duplicate_oils = ['ALVHEIM BLEND, STATOIL', 'DRAUGEN, STATOIL',
                  'EKOFISK BLEND 2000', 'EKOFISK BLEND, STATOIL',
                  'EKOFISK, CITGO', 'EKOFISK, EXXON', 'EKOFISK, PHILLIPS',
                  'EKOFISK, STATOIL', 'ELDFISK', 'ELDFISK B',
                  'GLITNE, STATOIL', 'GOLIAT BLEND, STATOIL',
                  'GRANE BLEND, STATOIL', 'GUDRUN BLEND, STATOIL',
                  'GULLFAKS A, STATOIL', 'GULLFAKS C, STATOIL',
                  'GULLFAKS, SHELL OIL', 'GULLFAKS SOR',
                  'GULLFAKS, STATOIL', 'HEIDRUN, STATOIL',
                  'NJORD, STATOIL', 'NORNE, STATOIL',
                  'OSEBERG BLEND, STATOIL', 'OSEBERG EXXON',
                  'OSEBERG, PHILLIPS', 'OSEBERG, SHELL OIL',
                  'SLEIPNER CONDENSATE, STATOIL',
                  'STATFJORD BLEND, STATOIL', 'VARG, STATOIL']


    def __init__(self, weathering_model='noaa', *args, **kwargs):
        self.oil_weathering_model = weathering_model

        if self.oil_weathering_model == 'noaa':  # Currently the only option
            self.oiltypes = adios.get_oil_names(
                location=kwargs.get('location', None))

            # Update config with oiltypes
            self.oiltypes.extend(adios.oil_name_alias.keys())

            # Sort alphabetically, but put GENERIC oils first
            generic_oiltypes = [o for o in self.oiltypes if o[0:7] == 'GENERIC']
            other_oiltypes = [o for o in self.oiltypes if o[0:7] != 'GENERIC']
            self.oiltypes = sorted([o for o in generic_oiltypes]) + sorted([o for o in other_oiltypes])
            self.oiltypes = [ot for ot in self.oiltypes if ot not in self.duplicate_oils]

            # For Norwegian oils, max water fraction from Sintef is overriding NOAA value
            self.max_water_fraction = None
        else:
            raise ValueError('Weathering model unknown: ' + weathering_model)

        kwargs.pop('location', None)

        # Calling general constructor of parent class
        super(OpenOil, self).__init__(*args, **kwargs)

        self._add_config({
            'seed:m3_per_hour': {
                'type': 'float',
                'default': 1,
                'min': 0,
                'max': 1e10,
                'units': 'm3 per hour',
                'description':
                'The amount (volume) of oil released per hour (or total amount if release is instantaneous)',
                'level': CONFIG_LEVEL_ESSENTIAL
            },
            'seed:droplet_size_distribution': {
                'type':
                'enum',
                'enum': ['uniform', 'normal', 'lognormal'],
                'default':
                'uniform',
                'level':
                CONFIG_LEVEL_ADVANCED,
                'description':
                'Droplet size distribution used for subsea release.'
            },
            'seed:droplet_diameter_mu': {
                'type': 'float',
                'default': 0.001,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description':
                'The mean diameter of oil droplet for a subsea release, used in normal/lognormal distributions.',
                'level': CONFIG_LEVEL_BASIC
            },
            'seed:droplet_diameter_sigma': {
                'type': 'float',
                'default': 0.0005,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description':
                'The standard deviation in diameter of oil droplet for a subsea release, used in normal/lognormal distributions.',
                'level': CONFIG_LEVEL_BASIC
            },
            'seed:droplet_diameter_min_subsea': {
                'type': 'float',
                'default': 0.0005,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description':
                'The minimum diameter of oil droplet for a subsea release, used in unifrom distribution.',
                'level': CONFIG_LEVEL_BASIC
            },
            'seed:droplet_diameter_max_subsea': {
                'type': 'float',
                'default': 0.005,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description':
                'The maximum diameter of oil droplet for a subsea release, used in uniform distribution.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:dispersion': {
                'type': 'bool',
                'default': True,
                'description':
                'Oil is removed from simulation (dispersed), if entrained as very small droplets.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:evaporation': {
                'type': 'bool',
                'default': True,
                'description': 'Surface oil is evaporated.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:emulsification': {
                'type': 'bool',
                'default': True,
                'description':
                'Surface oil is emulsified, i.e. water droplets are mixed into oil due to wave mixing, with resulting increas of viscosity.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:biodegradation': {
                'type': 'bool',
                'default': False,
                'description': 'Oil mass is biodegraded (eaten by bacteria).',
                'level': CONFIG_LEVEL_BASIC
            },
            'biodegradation:method': {
                'type': 'enum',
                'enum': ['Adcroft', 'half_time'],
                'default': 'Adcroft', 'level': CONFIG_LEVEL_ADVANCED,
                'description': 'Alogorithm to be used for biodegradation of oil'
            },
            'processes:update_oilfilm_thickness': {
                'type': 'bool',
                'default': False,
                'description':
                'Oil film thickness is calculated at each time step. The alternative is that oil film thickness is kept constant with value provided at seeding.',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'wave_entrainment:droplet_size_distribution': {
                'type':
                'enum',
                'enum': ['Johansen et al. (2015)', 'Li et al. (2017)'],
                'default':
                'Johansen et al. (2015)',
                'level':
                CONFIG_LEVEL_ADVANCED,
                'description':
                'Algorithm to be used for calculating oil droplet size spectrum after entrainment by breaking waves.'
            },
            'wave_entrainment:entrainment_rate': {
                'type':
                'enum',
                'enum': ['Li et al. (2017)'],
                'default':
                'Li et al. (2017)',
                'level':
                CONFIG_LEVEL_ADVANCED,
                'description':
                'Algorithm to be used for calculating the entrainment rate of oil due to wave breaking.'
            },
            'seed:oil_type': {
                'type':
                'enum',
                'enum':
                self.oiltypes,
                'default':
                self.oiltypes[0],
                'level':
                CONFIG_LEVEL_ESSENTIAL,
                'description':
                'Oil type to be used for the simulation, from the NOAA ADIOS database.'
            },
        })

        self._set_config_default('drift:vertical_advection', False)
        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:current_uncertainty', 0.05)
        self._set_config_default('drift:wind_uncertainty', 0.5)
        self._set_config_default('drift:max_speed', 1.3)

    def update_surface_oilfilm_thickness(self):
        '''The mass of oil is summed within a grid of 100x100
        cells covering the oil at a given time. Each oil particle
        within each cell is given a film thickness as the amount of
        oil divided by the cell area.
        '''
        from scipy.stats import binned_statistic_2d
        surface = np.where(self.elements.z == 0)[0]
        if len(surface) == 0:
            logger.debug('No oil at surface, no film thickness to update')
            return
        logger.debug(
            'Updating oil film thickness for %s of %s elements at surface' %
            (len(surface), self.num_elements_active()))
        meanlon = self.elements.lon[surface].mean()
        meanlat = self.elements.lat[surface].mean()
        # Using stereographic coordinates to get regular X and Y
        psproj = pyproj.Proj('+proj=stere +lat_0=%s +lat_ts=%s +lon_0=%s' %
                             (meanlat, meanlat, meanlon))
        X, Y = psproj(self.elements.lon[surface], self.elements.lat[surface])
        mass_bin, x_edge, y_edge, binnumber = binned_statistic_2d(
            X,
            Y,
            self.elements.mass_oil[surface],
            expand_binnumbers=True,
            statistic='sum',
            bins=100)
        bin_area = (x_edge[1] - x_edge[0]) * (y_edge[1] - y_edge[0])
        oil_density = 1000  # ok approximation here
        film_thickness = (mass_bin / oil_density) / bin_area
        # Postulating min and max film thickness
        max_thickness = 0.01  # 1 cm
        min_thickness = 1e-9  # 1 nanometer
        if film_thickness.max() > max_thickness:
            logger.debug(
                'Warning: decreasing thickness to %sm for %s of %s bins' %
                (max_thickness, np.sum(film_thickness > max_thickness),
                 film_thickness.size))
            film_thickness[film_thickness > max_thickness] = max_thickness
        num_too_thin = np.sum((film_thickness < min_thickness)
                              & (film_thickness > 0))
        if num_too_thin > 0:
            logger.debug(
                'Warning: increasing thickness to %sm for %s of %s bins' %
                (min_thickness, num_too_thin, film_thickness.size))
            film_thickness[film_thickness < min_thickness] = min_thickness

        # https://github.com/scipy/scipy/issues/7010
        binnumber = binnumber - 1

        bx = binnumber[0, :]
        by = binnumber[1, :]
        # Update thickness
        self.elements.oil_film_thickness[
            surface] = self.elements.oil_film_thickness[surface] * np.nan
        self.elements.oil_film_thickness[surface] = \
            film_thickness[bx, by]

    def biodegradation(self):
        if self.get_config('processes:biodegradation') is True:
            method = self.get_config('biodegradation:method')
            logger.debug(f'Calculating: biodegradation ({method})')
            if method == 'Adcroft':
                self.biodegradation_adcroft()
            elif method == 'half_time':
                self.biodegradation_half_time()

    def biodegradation_half_time(self):
        '''Oil biodegradation with exponential decay'''

        surface = np.where(self.elements.z == 0)[0]  # of active elements
        age0 = self.time_step.total_seconds()/(3600*24)  # days

        half_time = self.elements.biodegradation_half_time_droplet.copy()
        half_time[surface] = self.elements.biodegradation_half_time_slick[surface]

        fraction_biodegraded = (1 - np.exp(-age0 / half_time))
        biodegraded_now = self.elements.mass_oil * fraction_biodegraded

        self.elements.mass_biodegraded = \
            self.elements.mass_biodegraded + biodegraded_now
        self.elements.mass_oil = \
            self.elements.mass_oil - biodegraded_now
        if self.oil_weathering_model == 'noaa':
            self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :] = \
            self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :]*(1-fraction_biodegraded[:, np.newaxis])

    def biodegradation_adcroft(self):
        '''
        Oil biodegradation function based on the article:
        Adcroft et al. (2010), Simulations of underwater plumes of
        dissolved oil in the Gulf of Mexico.
        '''

        swt = self.environment.sea_water_temperature.copy()
        swt[swt > 100] -= 273.15  # K to C
        age0 = self.time_step.total_seconds() / (3600 * 24)

        # Decay rate in days (temperature in Celsius)
        tau = (12) * (3**((20 - swt) / 10))

        fraction_biodegraded = (1 - np.exp(-age0 / tau))
        biodegraded_now = self.elements.mass_oil * fraction_biodegraded

        self.elements.mass_biodegraded = \
            self.elements.mass_biodegraded + biodegraded_now
        self.elements.mass_oil = \
            self.elements.mass_oil - biodegraded_now
        if self.oil_weathering_model == 'noaa':
            self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :] = \
            self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :]*(1-fraction_biodegraded[:, np.newaxis])

    def disperse(self):
        if self.get_config('processes:dispersion') is True:

            logger.debug('   Calculating: dispersion')
            windspeed = np.sqrt(self.environment.x_wind**2 +
                                self.environment.y_wind**2)
            # From NOAA PyGnome model:
            # https://github.com/NOAA-ORR-ERD/PyGnome/
            v_entrain = 3.9E-8
            sea_water_density = 1028
            fraction_breaking_waves = 0.02
            wave_significant_height = \
                self.significant_wave_height()
            wave_significant_height[wave_significant_height == 0] = \
                0.0246*windspeed[wave_significant_height == 0]**2
            dissipation_wave_energy = \
                (0.0034*sea_water_density*9.81*wave_significant_height**2)
            c_disp = np.power(dissipation_wave_energy, 0.57) * \
                fraction_breaking_waves
            # Roy's constant
            C_Roy = 2400.0 * np.exp(-73.682 * np.sqrt(
                self.elements.viscosity / self.elements.density))

            q_disp = C_Roy * c_disp * v_entrain / self.elements.density

            oil_mass_loss = (q_disp * self.time_step.total_seconds() *
                             self.elements.density) * self.elements.mass_oil

            self.elements.mass_oil -= oil_mass_loss
            self.elements.mass_dispersed += oil_mass_loss

            #self.elements.z = \
            #    self.environment.sea_surface_wave_significant_height

            ## Marks R. (1987), Marine aerosols and whitecaps in the
            ## North Atlantic and Greenland sea regions
            ## Deutsche Hydrografische Zeitschrift, Vol 40, Issue 2 , pp 71-79
            #whitecap_coverage = (2.54E-4)*np.power(windspeed, 3.58)

            ## Martinsen et al. (1994), The operational
            ## oil drift system at DNMI
            ## DNMI Technical report No 125, 51
            #wave_period = 3.85*np.sqrt(
            #    self.environment.sea_surface_wave_significant_height)  # sec

            #time_between_breaking_events = 3.85/whitecap_coverage  # TBC

            #rho_w = 1025  # kg/m3
            #wave_period[wave_period == 0] = 5  # NB: temporal fix for no waves
            #dissipation_energy = 0.0034*rho_w*9.81*(wave_period**2)
            #dsize_coeff = 2100  # Wettre et al., Appendix A
            #c_oil = dsize_coeff*np.power(self.elements.viscosity, -0.4)
            ## Random numbers between 0 and 1:
            #p = np.random.rand(self.num_elements_active(), )
            #oil_per_unit_surface = 1
            #droplet_size = np.power((p*oil_per_unit_surface)/(c_oil *
            #                        np.power(dissipation_energy, 0.57)),
            #                        1/1.17)
            #self.deactivate_elements(droplet_size < 5E-7, reason='dispersed')

    def oil_weathering(self):
        if self.time_step.days < 0:
            logger.debug('Skipping oil weathering for backwards run')
            return
        self.timer_start('main loop:updating elements:oil weathering')
        if self.oil_weathering_model == 'noaa':
            self.oil_weathering_noaa()
        self.timer_end('main loop:updating elements:oil weathering')

    def prepare_run(self):
        if self.oil_weathering_model == 'noaa':
            self.noaa_mass_balance = {}
            # Populate with seeded mass spread on oiltype.mass_fraction
            mass_oil = np.atleast_1d(self.elements_scheduled.mass_oil)
            if len(mass_oil) == 1:
                mass_oil = mass_oil * np.ones(self.num_elements_total())
            self.noaa_mass_balance['mass_components'] = \
                np.asarray(self.oiltype.mass_fraction)*(mass_oil.reshape(
                    (self.num_elements_total(), -1)))
            self.noaa_mass_balance['mass_evaporated'] = \
                self.noaa_mass_balance['mass_components']*0

            self.oil_water_interfacial_tension = \
                self.oiltype.oil_water_surface_tension()
            logger.info('Oil-water surface tension is %f Nm' %
                        self.oil_water_interfacial_tension)
        try:
            max_water_fractions = json.loads(
                    resources.read_text('opendrift.models.openoil.adios', 'max_water_fraction.json'))
            if self.oil_name in max_water_fractions:
                self.max_water_fraction = max_water_fractions[self.oil_name]
                T = self.max_water_fraction['temperatures']
                wf = self.max_water_fraction['max_water_fraction']
                logger.info(f'Using max water fractions {wf} for temperatures {T} for oiltype {self.oil_name}')
                logger.info('Corresponding max water fraction from GNOME is '
                            f'{self.oiltype.gnome_oil["emulsion_water_fraction_max"]}')
            else:
                logger.info(f'Max water fraction not available for {self.oil_name}, using default')
        except Exception as e:
            logger.warning('Could not load max water content file')
            print(e)

        super(OpenOil, self).prepare_run()

    def oil_weathering_noaa(self):
        '''Oil weathering scheme adopted from NOAA PyGNOME model:
        https://github.com/NOAA-ORR-ERD/PyGnome
        '''
        logger.debug('NOAA oil weathering')
        # C to K
        self.environment.sea_water_temperature[
            self.environment.sea_water_temperature < 100] += 273.15

        #########################################################
        # Update density and viscosity according to temperature
        #########################################################
        self.timer_start(
            'main loop:updating elements:oil weathering:updating viscosities')
        oil_viscosity = self.KinematicViscosity.at_temp(
                self.environment.sea_water_temperature)
        self.timer_end(
            'main loop:updating elements:oil weathering:updating viscosities')
        self.timer_start(
            'main loop:updating elements:oil weathering:updating densities')
        oil_density = self.Density.at_temp(self.environment.sea_water_temperature)
        self.timer_end(
            'main loop:updating elements:oil weathering:updating densities')

        # Calculate emulsion density
        self.elements.density = (
            self.elements.water_fraction * self.sea_water_density() +
            (1 - self.elements.water_fraction) * oil_density)

        # Calculate emulsion viscosity
        visc_f_ref = 0.84  # From PyGNOME
        visc_curvfit_param = 1.5e3  # units are sec^0.5 / m
        fw_d_fref = self.elements.water_fraction / visc_f_ref
        kv1 = np.sqrt(oil_viscosity) * visc_curvfit_param
        kv1[kv1 < 1] = 1
        kv1[kv1 > 10] = 10
        # NB: neglecting dispersed and biodegraded in calculation of fraction_evaporated
        self.elements.fraction_evaporated = self.elements.mass_evaporated / (
            self.elements.mass_oil + self.elements.mass_evaporated)

        self.elements.viscosity = (
            oil_viscosity * np.exp(kv1 * self.elements.fraction_evaporated) * 
                (1 + (fw_d_fref / (1.187 - fw_d_fref)))**2.49)

        if self.get_config('processes:evaporation') is True:
            self.timer_start(
                'main loop:updating elements:oil weathering:evaporation')
            self.evaporation_noaa()
            self.timer_end(
                'main loop:updating elements:oil weathering:evaporation')

        if self.get_config('processes:emulsification') is True:
            self.timer_start(
                'main loop:updating elements:oil weathering:emulsification')
            self.emulsification_noaa()
            self.timer_end(
                'main loop:updating elements:oil weathering:emulsification')

        if self.get_config('processes:dispersion') is True:
            self.timer_start(
                'main loop:updating elements:oil weathering:dispersion')
            self.disperse_noaa()
            self.timer_end(
                'main loop:updating elements:oil weathering:dispersion')

        if self.get_config('processes:biodegradation') is True:
            self.timer_start(
                'main loop:updating elements:oil weathering:biodegradation')
            self.biodegradation()
            self.timer_end(
                'main loop:updating elements:oil weathering:biodegradation')

        if self.elements.mass_oil.min() < 0:  # Should not happen
            logger.warning('NEGATIVE OIL MASS!')

    def disperse_noaa(self):
        logger.debug('    Calculating: dispersion - NOAA')
        # From NOAA PyGnome model:
        # https://github.com/NOAA-ORR-ERD/PyGnome/
        c_disp = np.power(self.wave_energy_dissipation(), 0.57) * \
            self.sea_surface_wave_breaking_fraction()
        # Roy's constant
        C_Roy = 2400.0 * np.exp(
            -73.682 * np.sqrt(self.elements.viscosity / self.elements.density))
        v_entrain = 3.9E-8
        q_disp = C_Roy * c_disp * v_entrain / self.elements.density
        fraction_dispersed = (q_disp * self.time_step.total_seconds() *
                              self.elements.density)
        if fraction_dispersed.max() >= 1:
            logger.warning(
                'Dispersion fraction larger than 1 -> setting to 0.99')
            fraction_dispersed[fraction_dispersed >= 1] = .99
        oil_mass_loss = fraction_dispersed * self.elements.mass_oil

        self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :] = \
            self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :]*(1-fraction_dispersed[:, np.newaxis])

        self.elements.mass_oil -= oil_mass_loss
        self.elements.mass_dispersed += oil_mass_loss

    def plot_droplet_spectrum(self):
        '''Plotting distribution of droplet radii, for debugging'''
        plt.hist(self.elements.diameter / 2.0)
        plt.show()

    def evaporation_noaa(self):
        #############################################
        # Evaporation, for elements at surface only
        #############################################
        logger.debug('    Calculating evaporation - NOAA')
        surface = np.where(self.elements.z == 0)[0]  # of active elements
        if len(surface) == 0:
            logger.debug('All elements submerged, no evaporation')
            return
        if self.elements.age_seconds[surface].min() > 3600 * 24:
            logger.debug('All surface oil elements older than 24 hours, ' +
                         'skipping further evaporation.')
            return
        surfaceID = self.elements.ID[surface] - 1  # of any elements
        # Area for each element, repeated for each component
        volume = (self.elements.mass_oil[surface] /
                  self.elements.density[surface])
        area = volume / self.elements.oil_film_thickness[surface]
        evap_decay_constant = noaa.evap_decay_constant(
            self.oiltype,
            self.wind_speed()[surface],
            self.environment.sea_water_temperature[surface], area,
            self.noaa_mass_balance['mass_components'][surfaceID, :])
        mass_remain = \
            (self.noaa_mass_balance['mass_components'][surfaceID, :] *
             np.exp(evap_decay_constant*self.time_step.total_seconds()))
        self.elements.mass_evaporated[surface] += \
            np.sum(self.noaa_mass_balance[
                    'mass_components'][surfaceID, :] - mass_remain, 1)
        self.noaa_mass_balance['mass_components'][surfaceID, :] = \
            mass_remain
        self.elements.mass_oil[surface] = np.sum(mass_remain, 1)

    def emulsification_noaa(self):
        #############################################
        # Emulsification (surface only)
        #############################################
        logger.debug('    Calculating emulsification - NOAA')
        emul_time = self.oiltype.bulltime
        emul_constant = self.oiltype.bullwinkle

        # Emulsify...
        # f ((le_age >= emul_time && emul_time >= 0.) || frac_evap[i] >= emul_C && emul_C > 0.)

        start_emulsion = np.where((
            (self.elements.age_seconds >= emul_time) & (emul_time >= 0))
                                  | ((self.elements.fraction_evaporated >= emul_constant)
                                     & (emul_constant > 0)))[0]
        if len(start_emulsion) == 0:
            logger.debug('        Emulsification not yet started')
            return

        # max water content fraction - get from database
        Y_max = np.atleast_1d(self.oiltype.emulsion_water_fraction_max)
        if self.max_water_fraction is not None:
            wf = self.max_water_fraction['max_water_fraction']
            wft = self.max_water_fraction['temperatures']
            if len(wf) == 1:
                wf = [wf, wf]
                wft = [wft, wft]
            swt = self.environment.sea_water_temperature[start_emulsion] - 273.15  # to Celcius
            weights = (wft[1] - swt) / (wft[1] - wft[0])
            weights[swt>wft[1]] = 0
            weights[swt<=wft[0]] = 1
            max_water_fraction_sintef = weights*wf[0] + (1-weights)*wf[1]

            if (Y_max - max_water_fraction_sintef).min() > 0:
                logger.debug(
                    f'Overriding max water fraction {Y_max} with linear fit to SINTEF max values:'
                    f' T: {wft}, Fraction: {wf}')
                Y_max = np.array(np.minimum(Y_max, max_water_fraction_sintef))
        # emulsion
        if Y_max.max() <= 0:
            logger.debug('Oil does not emulsify, returning.')
            return
        # Constants for droplets
        drop_min = 1.0e-6
        drop_max = 1.0e-5
        S_max = (6. / drop_min) * (Y_max / (1.0 - Y_max))
        S_min = (6. / drop_max) * (Y_max / (1.0 - Y_max))
        if self.oiltype.bulltime > 0:  # User has set value
            start_time = self.oiltype.bulltime * np.ones(len(start_emulsion))
        else:
            start_time = self.elements.age_seconds[start_emulsion]
            # TODO: it should be possible to specify oil age at seeding
            start_time[self.elements.age_seconds[start_emulsion] >=
                       0] = self.elements.bulltime[start_emulsion]
        # Update droplet interfacial area
        k_emul = noaa.water_uptake_coefficient(self.oiltype, self.wind_speed()[start_emulsion])
        self.elements.interfacial_area[start_emulsion] = \
            self.elements.interfacial_area[start_emulsion] + \
            (k_emul*self.time_step.total_seconds()* np.exp((-k_emul/S_max)*(
                self.elements.age_seconds[start_emulsion] - start_time)))
        self.elements.interfacial_area[start_emulsion] = np.minimum(self.elements.interfacial_area[start_emulsion], S_max)
        # Update water fraction
        self.elements.water_fraction[start_emulsion] = (
            self.elements.interfacial_area[start_emulsion] * drop_max / (6.0 +
             (self.elements.interfacial_area[start_emulsion] * drop_max)))
        self.elements.water_fraction[start_emulsion] = np.minimum(self.elements.water_fraction[start_emulsion], Y_max)

    def update_terminal_velocity(self,
                                 Tprofiles=None,
                                 Sprofiles=None,
                                 z_index=None):
        """Calculate terminal velocity for oil droplets

        according to:

        * Tkalich et al. (2002): Vertical mixing of oil droplets by breaking waves

        * Marine Pollution Bulletin 44, 1219-1229

        If profiles of temperature and salt are passed into this function,
        they will be interpolated from the profiles.
        if not, T,S will be fetched from reader.
        """
        g = 9.81  # ms-2

        r = self.elements.diameter  # NB: r is diameter, not radius

        # Prepare interpolation of temp, salt

        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i,
                                   bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.uint8), 0)
            lower = np.minimum(upper + 1, Tprofiles.shape[0] - 1)
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

        T0 = T0 - 273.15 # convert to Celcius - needed for the calcs in this method

        rho_oil = self.elements.density
        rho_water = self.sea_water_density(T=T0, S=S0)

        # dynamic water viscosity
        my_w = 0.001 * (1.7915 - 0.0538 * T0 + 0.007 *
                        (T0**(2.0)) - 0.0023 * S0)
        # ~0.0014 kg m-1 s-1
        # kinemativ water viscosity
        ny_w = my_w / rho_water
        rhopr = rho_oil / rho_water

        # terminal velocity for low Reynolds numbers
        kw = 2 * g * (1 - rhopr) / (9 * ny_w)
        W = kw * (r/2)**2 # r is diameter so divide by 2 for radius

        # check if we are in a high Reynolds number regime
        Re = r * W / ny_w # r is diameter so no need to multiply by 2
        highRe = np.where(Re > 50)

        # Terminal velocity in high Reynolds numbers
        kw = (16 * g * (1 - rhopr) / 3)**0.5
        W2 = kw * (r/2)**0.5

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    def oil_wave_entrainment_rate(self):
        er = self.get_config('wave_entrainment:entrainment_rate')
        if er == 'Li et al. (2017)':
            entrainment_rate = oil_wave_entrainment_rate_li2017(
                dynamic_viscosity=self.elements.viscosity *
                self.elements.density,
                oil_density=self.elements.density,
                interfacial_tension=self.oil_water_interfacial_tension,
                significant_wave_height=self.significant_wave_height(),
                wave_breaking_fraction=self.sea_surface_wave_breaking_fraction(
                ),
                sea_water_density=self.sea_water_density())
            return entrainment_rate
        else:
            logger.error("no entrainment rate mechanism configured")
            raise Exception("no entrainment rate mechanism configured")

    def prepare_vertical_mixing(self):
        '''Calculate entrainment probability before main loop'''
        self.oil_entrainment_probability = \
            1 - np.exp(-self.oil_wave_entrainment_rate()*\
                       self.get_config('vertical_mixing:timestep'))
        # Calculate a random droplet diameter for each particle,
        # to be used if this particle gets entrained
        self.droplet_diameter_if_entrained = \
            self.get_wave_breaking_droplet_diameter()
        # Uncomment lines below to plot droplet size distribution at each step
        #import matplotlib.pyplot as plt
        #plt.hist(self.droplet_diameter_if_entrained, 200)
        #plt.gca().set_xscale("log")
        #plt.gca().set_yscale("log")
        #plt.show()

    def surface_wave_mixing(self, time_step_seconds):
        """Mix surface oil into water column."""
        # Entrain oil into uppermost layer (whitecapping from waves)
        # TODO: optimise this by only calculate for surface elements
        surface = self.elements.z >= 0
        random_number = np.random.uniform(0, 1, len(self.elements.z))
        entrained = np.logical_and(
            surface, random_number < self.oil_entrainment_probability)

        # Intrusion depth for wave entrainment from
        # Delvigne and Sweeney (1988), Li et al. (2017):
        if entrained.sum() > 0:
            logger.debug('Entraining %i of %i surface elements' %
                         (entrained.sum(), surface.sum()))
            zb = 1.5 * self.significant_wave_height()  # between 0 and zb
            intrusion_depth = np.random.uniform(0, np.mean(zb),
                                                entrained.sum())
            self.elements.z[entrained] = -intrusion_depth
            if self.keep_droplet_diameter is False:
                # Give entrained elements a random diameter
                self.elements.diameter[entrained] = \
                    self.droplet_diameter_if_entrained[entrained]

    def surface_stick(self):
        """set surfaced particles to exactly zero depth to let them form a slick """

        surface = np.where(self.elements.z >= 0)
        if len(surface[0]) > 0:
            self.elements.z[surface] = 0.

    def get_wave_breaking_droplet_diameter(self):
        dm = self.get_config('wave_entrainment:droplet_size_distribution')
        if dm == 'Johansen et al. (2015)':
            return self.get_wave_breaking_droplet_diameter_johansen2015()
        elif dm == 'Li et al. (2017)':
            return self.get_wave_breaking_droplet_diameter_liz2017()
        else:
            raise Exception("no wave entrainment droplet size distribution specified")

    def get_wave_breaking_droplet_diameter_liz2017(self):
        # Li,Zhengkai, M. Spaulding, D. French-McCay, D. Crowley, J.R. Payne: "Development of a unified oil droplet size distribution model
        # with application to surface breaking waves and subsea blowout releases considering dispersant effects" Mar. Pol. Bul.
        # DOI: 10.1016/j.marpolbul.2016.09.008
        # Should be prefered when the oil film thickness is unknown.
        if not hasattr(self, 'droplet_spectrum_pdf'):
            # Generate droplet spectrum as in Li (Zhengkai) et al. (2017)
            # Bounds are hardcoded to 1 micron and 3mm
            logger.debug('Generating wave breaking droplet size spectrum')
            self.droplet_spectrum_diameter = np.linspace(1e-6, 3e-3, 1000000)
            g = 9.81
            interfacial_tension = self.oil_water_interfacial_tension
            delta_rho = self.sea_water_density() - self.elements.density
            d_o = 4 * (interfacial_tension / (delta_rho * g))**0.5
            we = (self.sea_water_density() * g *
                  self.significant_wave_height() * d_o) / interfacial_tension
            oh = self.elements.viscosity * self.elements.density * (
                self.elements.density * interfacial_tension *
                d_o)**-0.5  # From kin. to dyn. viscosity by * density
            r = 1.791
            p = 0.460
            q = -0.518
            dV_50 = d_o * r * (
                1 + 10 * oh
            )**p * we**q  # median droplet diameter in volume distribution
            sd = 0.4  # log standard deviation in log10 units
            Sd = np.log(10) * sd  # log standard deviation in natural log units
            # TODO: calculation below with scalars, but we have arrays, with varying oil properties
            # treat all particle in one go:
            dV_50 = np.mean(dV_50)  # mean log diameter
            dN_50 = np.exp(
                np.log(dV_50) - 3 *
                Sd**2)  # convert number distribution to volume distribution
            logger.debug(
                'Droplet distribution median diameter dV_50: %f, dN_50: %f ' %
                (dV_50, np.mean(dN_50)))
            spectrum = (np.exp(
                -(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 /
                (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd *
                                 np.sqrt(2 * np.pi))
            self.droplet_spectrum_pdf = spectrum / np.sum(spectrum)
        if ~np.isfinite(np.sum(self.droplet_spectrum_pdf)) or \
                np.abs(np.sum(self.droplet_spectrum_pdf) - 1) > 1e-6:
            logger.warning('Could not update droplet diameters.')
            return self.elements.diameter
        else:
            return np.random.choice(self.droplet_spectrum_diameter,
                                    size=self.num_elements_active(),
                                    p=self.droplet_spectrum_pdf)

    def get_wave_breaking_droplet_diameter_johansen2015(self):
        # Johansen O, Reed M, Bodsberg NR, Natural dispersion revisited
        # DOI: 10.1016/j.marpolbul.2015.02.026
        # requires oil film thickness
        if not hasattr(self, 'droplet_spectrum_pdf') or self.get_config(
                'processes:update_oilfilm_thickness') is True:
            # Generate droplet spectrum as in Johansen et al. (2015)
            # Bounds are hardcoded to 1micron and 3mm
            logger.debug('Generating wave breaking droplet size spectrum')
            self.droplet_spectrum_diameter = np.linspace(1e-6, 3e-3, 1000000)
            g = 9.81
            interfacial_tension = self.oil_water_interfacial_tension
            H = self.significant_wave_height(
            )  # fall height = 2 * wave amplitude
            # Reyolds number (Eq. 7a from Johansen et al. 2015)
            re = (self.elements.density * self.elements.oil_film_thickness *
                  (g * H)**0.5) / (self.elements.viscosity *
                                   self.elements.density)
            # Weber number (Eq. 7b from Johansen et al.2015)
            we = (self.elements.density * self.elements.oil_film_thickness *
                  g * H) / interfacial_tension  # Weber number
            A = 2.251  # parameters from Johansen et al. 2015
            Bp = 0.027
            B = A * Bp
            dN_50 = (A * self.elements.oil_film_thickness * we**-0.6) + (
                B * self.elements.oil_film_thickness * re**-0.6)
            # median droplet diameter in number distribution
            sd = 0.4  # log standard deviation in log10 units
            Sd = np.log(10) * sd  # log standard deviation in natural log units
            # Convert number distribution to volume distribution
            dV_50 = np.exp(np.log(dN_50) + 3 * Sd**2)
            # TODO: calculation below with scalars, but we have
            # arrays, with varying oil properties
            # treat all particle in one go:
            dV_50 = np.mean(dV_50)  # mean log diameter
            logger.debug(
                'Droplet distribution median diameter dV_50: %f, dN_50: %f ' %
                (dV_50, np.mean(dN_50)))
            spectrum = (np.exp(
                -(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 /
                (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd *
                                 np.sqrt(2 * np.pi))
            self.droplet_spectrum_pdf = spectrum / np.sum(spectrum)
        if ~np.isfinite(np.sum(self.droplet_spectrum_pdf)) or \
                np.abs(np.sum(self.droplet_spectrum_pdf) - 1) > 1e-6:
            logger.warning('Could not update droplet diameters.')
            return self.elements.diameter
        else:
            return np.random.choice(self.droplet_spectrum_diameter,
                                    size=self.num_elements_active(),
                                    p=self.droplet_spectrum_pdf)

    def resurface_elements(self, minimum_depth=None):
        """Oil elements reaching surface (or above) form slick, not droplet"""
        surface = np.where(self.elements.z >= 0)[0]
        self.elements.z[surface] = 0

    def advect_oil(self):

        # Calculating various drift factors according to ice concentration
        if hasattr(self.environment, 'sea_ice_area_fraction'):
            A = self.environment.sea_ice_area_fraction
            # According to
            # Nordam T, Beegle-Krause CJ, Skancke J, Nepstad R, Reed M.
            # Improving oil spill trajectory modelling in the Arctic.
            # Mar Pollut Bull. 2019;140:65-74.
            # doi:10.1016/j.marpolbul.2019.01.019
            k_ice = (A - 0.3) / (0.8 - 0.3)
            k_ice[A < 0.3] = 0
            k_ice[A > 0.8] = 1
            if k_ice.max() > 0:
                logger.info(
                    'Ice concentration above 30%, using Nordam scheme for advection in ice'
                )
            # Using decreased Stokes drift according to
            # Arneborg, L. (2017). Oil drift modellling in pack ice
            # - Sensitivity of oil-in-ice parameters.
            # Ocean Engineering 144 (2017) 340-350
            factor_stokes = (0.7 - A) / 0.7
            factor_stokes[A > 0.7] = 0
        else:
            k_ice = 0
            factor_stokes = 1

        # Simply move particles with ambient current
        self.advect_ocean_current(factor=1 - k_ice)

        # Wind drag for elements at ocean surface
        self.advect_wind(factor=1 - k_ice)

        # Stokes drift
        self.stokes_drift(factor_stokes)

        # Advect with ice
        self.advect_with_sea_ice(factor=k_ice)

    def update(self):
        """Update positions and properties of oil particles."""

        # TODO: move all config-checking inside respective methods
        if self.get_config('processes:update_oilfilm_thickness') is True:
            self.update_surface_oilfilm_thickness()

        # Oil weathering (inherited from OpenOil)
        self.oil_weathering()

        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
            del self.droplet_spectrum_pdf

        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()

        # Horizontal advection (inherited from OpenOil)
        self.advect_oil()

    def get_oil_budget(self):
        """Get oil budget for the current simulation

        The oil budget consists of the following categories:

        * surface: the sum of variable mass_oil for all active elements where z = 0
        * submerged: the sum of variable mass_oil for all active elements where z < 0
        * stranded: the sum of variable mass_oil for all elements which are stranded
        * evaporated: the sum of variable mass_evaporated for all elements
        * dispersed: the sum of variable mass_dispersed for all elements

        The sum (total mass) shall equal the mass released. Note that the mass of oil
        is conserved, whereas the volume may change with changes in density and
        water uptake (emulsification). Therefore mass should be used for budgets,
        eventually converted to volume (by dividing on density) in the final step
        before presentation.

        Note that mass_oil is the mass of pure oil. The mass of oil emulsion
        (oil containing entrained water droplets) can be calculated as:

        .. code::

            mass_emulsion = mass_oil / (1 - water_fraction)

        I.e. water_fraction = 0 means pure oil, water_fraction = 0.5 means mixture of
        50% oil and 50% water, and water_fraction = 0.9 (which is maximum)
        means 10% oil and 90% water.
        """

        if self.time_step.days < 0:  # Backwards simulation
            return None

        z, dummy = self.get_property('z')
        mass_oil, status = self.get_property('mass_oil')
        density = self.get_property('density')[0][0, 0]

        if 'stranded' not in self.status_categories:
            self.status_categories.append('stranded')
        mass_submerged = np.ma.masked_where(
            ((status == self.status_categories.index('stranded')) |
             (z == 0.0)), mass_oil)
        mass_submerged = np.ma.sum(mass_submerged, axis=1).filled(0)

        mass_surface = np.ma.masked_where(
            ((status == self.status_categories.index('stranded')) | (z < 0.0)),
            mass_oil)
        mass_surface = np.ma.sum(mass_surface, axis=1).filled(0)

        mass_stranded = np.ma.sum(np.ma.masked_where(
            status != self.status_categories.index('stranded'), mass_oil),
                                  axis=1).filled(0)
        mass_evaporated, status = self.get_property('mass_evaporated')
        mass_evaporated = np.sum(mass_evaporated, axis=1).filled(0)
        mass_dispersed, status = self.get_property('mass_dispersed')
        mass_dispersed = np.sum(mass_dispersed, axis=1).filled(0)
        mass_biodegraded, status = self.get_property('mass_biodegraded')
        mass_biodegraded = np.sum(mass_biodegraded, axis=1).filled(0)

        oil_budget = {
            'oil_density':
            density,
            'mass_dispersed':
            mass_dispersed,
            'mass_submerged':
            mass_submerged,
            'mass_surface':
            mass_surface,
            'mass_stranded':
            mass_stranded,
            'mass_evaporated':
            mass_evaporated,
            'mass_biodegraded':
            mass_biodegraded,
            'mass_total': (mass_dispersed + mass_submerged + mass_surface +
                           mass_stranded + mass_evaporated + mass_biodegraded)
        }

        return oil_budget

    def plot_oil_budget(self,
                        filename=None,
                        ax=None,
                        show_watercontent_and_viscosity=True,
                        show_wind_and_current=True):

        if self.time_step.days < 0:  # Backwards simulation
            fig = plt.figure(figsize=(10, 6.))
            plt.text(0.1, 0.5, 'Oil weathering deactivated for '
                     'backwards simulations')
            plt.axis('off')
            if filename is not None:
                plt.savefig(filename)
                plt.close()
            else:
                plt.show()
            return

        b = self.get_oil_budget()

        oil_budget = np.row_stack(
            (b['mass_dispersed'], b['mass_submerged'], b['mass_surface'],
             b['mass_stranded'], b['mass_evaporated'], b['mass_biodegraded']))
        oil_density = b['oil_density']

        budget = np.cumsum(oil_budget, axis=0)

        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds() / 3600. for t in time_relative])

        if ax is None:
            # Left axis showing oil mass
            nrows = 1
            if show_watercontent_and_viscosity is True:
                nrows = nrows + 1
            if show_wind_and_current is True:
                nrows = nrows + 1
            fig, axs = plt.subplots(
                nrows=nrows, ncols=1,
                figsize=(10, 6. + (nrows - 1) * 3))  # Suitable aspect ratio
            #ax1 = fig.add_subplot(nrows=nrows, 1, 1)
            if nrows == 1:
                ax1 = axs
            elif nrows >= 2:
                ax1 = axs[0]
                if show_watercontent_and_viscosity is True:
                    self.plot_oil_watercontent_and_viscosity(ax=axs[1], show=False)
                if show_wind_and_current is True:
                    self.plot_environment(ax=axs[nrows - 1], show=False)
        else:
            ax1 = ax

        # Hack: make some emply plots since fill_between does not support label
        if np.sum(b['mass_dispersed']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='darkslategrey',
                              label='dispersed'))
            ax1.fill_between(time, 0, budget[0, :], facecolor='darkslategrey')
        if np.sum(b['mass_submerged']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='darkblue',
                              label='submerged'))
            ax1.fill_between(time,
                             budget[0, :],
                             budget[1, :],
                             facecolor='darkblue')
        if np.sum(b['mass_surface']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0), 0, 0, color='royalblue',
                              label='surface'))
            ax1.fill_between(time,
                             budget[1, :],
                             budget[2, :],
                             facecolor='royalblue')
        if np.sum(b['mass_stranded']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0), 0, 0, color='black', label='stranded'))
            ax1.fill_between(time,
                             budget[2, :],
                             budget[3, :],
                             facecolor='black')
        if np.sum(b['mass_evaporated']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='skyblue',
                              label='evaporated'))
            ax1.fill_between(time,
                             budget[3, :],
                             budget[4, :],
                             facecolor='skyblue')
        if np.sum(b['mass_biodegraded']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='indigo',
                              label='biodegraded'))
            ax1.fill_between(time,
                             budget[4, :],
                             budget[5, :],
                             facecolor='indigo')

        ax1.set_ylim([0, budget.max()])
        ax1.set_xlim([0, time.max()])
        ax1.set_ylabel('Mass oil  [%s]' %
                       self.elements.variables['mass_oil']['units'])
        ax1.set_xlabel('Time  [hours]')
        # Right axis showing volume
        ax2 = ax1.twinx()
        mass_total = b['mass_total'][-1]
        ax2.set_ylim([0, mass_total / oil_density])
        ax2.set_ylabel('Volume oil [m3]')
        plt.title('%s (%.1f kg/m3) - %s to %s' %
                  (self.get_oil_name(), oil_density,
                   self.start_time.strftime('%Y-%m-%d %H:%M'),
                   self.time.strftime('%Y-%m-%d %H:%M')))
        # Shrink current axis's height by 10% on the bottom
        #box = ax1.get_position()
        #ax1.set_position(
        #    [box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
        #ax2.set_position(
        #    [box.x0, box.y0 + box.height * 0.5, box.width, box.height * 0.6])
        ax1.legend(bbox_to_anchor=(0., -0.1, 1., -0.04),
                   loc=1,
                   ncol=6,
                   mode="expand",
                   borderaxespad=0.,
                   fontsize=10)
        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()

    def get_oil_name(self):
        try:
            return self.get_config('seed:oil_type')
        except:  # fallback if importing old files
            return 'unknown oiltype'

    def cumulative_oil_entrainment_fraction(self):
        '''Returns the fraction of oil elements which has been entrained vs time'''
        z = self.get_property('z')[0].copy()
        z = np.ma.masked_where(z == 0, z)
        me = np.ma.notmasked_edges(z, axis=0)
        maskfirst = me[0][0]
        maskrow = me[0][1]
        z = z * 0
        for mf, mr in zip(maskfirst, maskrow):
            z[mf:z.shape[0], mr] = 1  # has been entrained
        totentrained = np.sum(z, 1)
        cumulative_fraction_entrained = np.sum(z, 1) / z.shape[1]
        return cumulative_fraction_entrained

    def plot_oil_watercontent_and_viscosity(self, ax=None, show=True):
        if ax is None:
            fig, ax = plt.subplots()
        import matplotlib.dates as mdates

        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds() / 3600. for t in time_relative])
        kin_viscosity = self.history['viscosity']
        dyn_viscosity = kin_viscosity * self.history['density'] * 1000  # unit of mPas
        dyn_viscosity_mean = dyn_viscosity.mean(axis=0)
        dyn_viscosity_std = dyn_viscosity.std(axis=0)
        density = self.history['density'].mean(axis=0)
        density_std = self.history['density'].std(axis=0)
        watercontent = self.history['water_fraction'].mean(axis=0)*100
        watercontent_std = self.history['water_fraction'].std(axis=0)*100

        ax.plot(time,
                dyn_viscosity_mean,
                'g',
                lw=2,
                label='Emulsion viscosity')
        ax.fill_between(time,
                        dyn_viscosity_mean - dyn_viscosity_std,
                        dyn_viscosity_mean + dyn_viscosity_std,
                        color='g',
                        alpha=0.5)
        ax.set_ylim([0, max(dyn_viscosity_mean + dyn_viscosity_std)])
        ax.set_ylabel(r'Emulsion viscosity  [cPoise] / [mPas]', color='g')
        ax.tick_params(axis='y', colors='g')

        axb = ax.twinx()
        axb.plot(time, watercontent, 'b', lw=2, label='Water content')
        axb.fill_between(time,
                         watercontent - watercontent_std,
                         watercontent + watercontent_std,
                         color='b', alpha=0.5)
        ax.set_xlim([0, time.max()])
        ax.set_xlabel('Time [hours]')
        axb.set_ylim([0, 100])
        axb.set_ylabel(r'Water content  [%]', color='b')
        axb.tick_params(axis='y', colors='b')

        ax.legend(loc='upper left')
        axb.legend(loc='upper center')
        if show is True:
            plt.show()

    def set_oiltype(self, oiltype):
        """
        Sets the oil type by specifying the name, the first match will be chosen. See the `ADIOS database <https://adios.orr.noaa.gov/oils>`_ for a list. OpenDrift provides a small set of extra oils.
        """

        if self.get_config('seed:oil_type') != oiltype:
            self.__set_seed_config__('seed:oil_type', oiltype)
            logger.info(f'setting oil_type to: {oiltype}')

        oiltype = adios.oil_name_alias.get(oiltype, oiltype)
        self.oil_name = oiltype

        if self.oil_weathering_model == 'noaa':
            self.oiltype = adios.find_full_oil_from_name(self.oil_name)
            if not self.oiltype.valid():
                logger.error(
                    f"{self.oiltype} is not a valid oil for Opendrift simulations"
                )
                raise ValueError()
        else:
            raise ValueError("unsupported oil weathering model")

    def set_oiltype_by_id(self, oiltypeid):
        """
        Sets the oil type by specifying the ADIOS ID. See the `ADIOS database <https://adios.orr.noaa.gov/oils>`_ for a list. OpenDrift provides a small set of extra oils.
        """
        if self.oil_weathering_model == 'noaa':
            self.oiltype = adios.get_full_oil_from_id(oiltypeid)
            self.oil_name = self.oiltype.name
            if not self.oiltype.valid():
                logger.error(
                    f"{self.oiltype} is not a valid oil for Opendrift simulations"
                )
                raise ValueError()
        else:
            raise ValueError("unsupported oil weathering model")

    def set_oiltype_by_json(self, json):
        """
        Sets the oil type by specifing a JSON dict. The format should be the same as the ADIOS database. See the `ADIOS database <https://adios.orr.noaa.gov/oils>`_ for a list.
        """
        if self.oil_weathering_model == 'noaa':
            #o = { 'data': { 'attributes' : json } }
            #o['data']['_id'] = o['data']['attributes']['oil_id']
            #o['data']['attributes']['metadata']['location'] = 'Norway'

            self.oiltype = adios.oil.OpendriftOil(json)
            self.oil_name = self.oiltype.name
            if not self.oiltype.valid():
                logger.error(
                    f"{self.oiltype} is not a valid oil for Opendrift simulations"
                )
                raise ValueError()
        else:
            raise ValueError("unsupported oil weathering model")

    def set_oiltype_from_file(self, path):
        """
        Sets the oil type by specifing a JSON file. The format should be the same as the ADIOS database. See the `ADIOS database <https://adios.orr.noaa.gov/oils>`_ for a list.

        >>> o = OpenOil()
        >>> o.set_oiltype_from_file('opendrift/models/openoil/adios/extra_oils/AD04001.json')
        """
        if self.oil_weathering_model == 'noaa':
            import json
            with open(path, 'r') as fd:
                j = json.load(fd)

            self.set_oiltype_by_json(j)

        else:
            raise ValueError("unsupported oil weathering model")

    def store_oil_seed_metadata(self, **kwargs):
        for s in [
                'lon', 'lat', 'radius', 'time', 'number', 'z', 'm3_per_hour'
        ]:
            if not 'seed_' + s in self.metadata_dict:
                if s in kwargs:
                    val = kwargs[s]
                else:
                    if s == 'radius':
                        val = 0  # There is no default radius
                    elif s == 'z' and 'z' not in kwargs and self.get_config(
                            'seed:seafloor') is True:
                        val = 'seafloor'
                    else:
                        val = self.get_config('seed:' + s)
                if s == 'time':
                    if hasattr(kwargs[s], '__len__'):
                        self.add_metadata('seed_time', val[0])
                    else:
                        self.add_metadata('seed_time', val)
                elif isinstance(val, str):
                    self.add_metadata('seed_' + s, val)
                else:
                    self.add_metadata('seed_' + s, np.atleast_1d(val).mean())
        if not 'seed_oiltype' in self.metadata_dict:
            if 'oiltype' in kwargs:
                oiltype = kwargs['oiltype']
            elif 'oil_type' in kwargs:
                oiltype = kwargs['oil_type']
            else:
                oiltype = self.get_config('seed:oil_type')
            self.add_metadata('seed_oiltype', oiltype)

    def seed_elements(self, *args, **kwargs):

        if len(args) == 2:
            kwargs['lon'] = args[0]
            kwargs['lat'] = args[1]
            args = {}

        if 'number' not in kwargs:
            number = self.get_config('seed:number')
        else:
            number = kwargs['number']
        if 'diameter' in kwargs:
            logger.info('Droplet diameter is provided, and will '
                        'be kept constant during simulation')
            self.keep_droplet_diameter = True
        else:
            self.keep_droplet_diameter = False
        if 'z' not in kwargs or kwargs['z'] is None:
            if self.get_config('seed:seafloor') is True:
                kwargs['z'] = 'seafloor'
            else:
                kwargs['z'] = self.get_config('seed:z')
        if isinstance(kwargs['z'], str) and \
                kwargs['z'][0:8] == 'seafloor':
            z = -np.ones(number)
        else:
            z = np.atleast_1d(kwargs['z'])
        if len(z) == 1:
            z = z * np.ones(number)  # Convert scalar z to array
        subsea = z < 0
        if np.sum(subsea) > 0 and 'diameter' not in kwargs:
            dsd = self.get_config('seed:droplet_size_distribution')
            if dsd == 'uniform':
                # Droplet min and max for particles seeded below sea surface
                sub_dmin = self.get_config('seed:droplet_diameter_min_subsea')
                sub_dmax = self.get_config('seed:droplet_diameter_max_subsea')
                logger.info('Using uniform droplet size distribution between %s and %s m for '
                            'elements seeded below sea surface.' %
                            (sub_dmin, sub_dmax))
                kwargs['diameter'] = np.random.uniform(sub_dmin, sub_dmax, number)
            elif dsd == 'normal':
                # Droplet mu and sigma for particles seeded below sea surface
                sub_mu = self.get_config('seed:droplet_diameter_mu')
                sub_sigma = self.get_config('seed:droplet_diameter_sigma')
                logger.info('Using normal droplet size distribution with '
                            'mu = %s and sigma = %s m for elements seeded below sea surface.' %
                            (sub_mu, sub_sigma))
                kwargs['diameter'] = np.random.normal(sub_mu, sub_sigma, number)
            elif dsd == 'lognormal':
                # Droplet mu and sigma for particles seeded below sea surface
                sub_mu = self.get_config('seed:droplet_diameter_mu')
                sub_sigma = self.get_config('seed:droplet_diameter_sigma')
                logger.info('Using lognormal droplet size distribution with '
                            'mu = %s and sigma = %s m for elements seeded below sea surface.' %
                            (sub_mu, sub_sigma))
                # From numpy.random.lognormal:
                # "Note that the mean and standard deviation are not the values for the
                # distribution itself, but of the underlying normal distribution it is derived from."
                # So we need to compute the input to the function from the mean and
                # standard deviation of the data we want to generate (assumed as input)
                sub_sigma2 = sub_sigma**2
                sub_sigma2_lognormal = np.log(sub_sigma2/sub_mu**2 + 1)
                sub_mu_lognormal = np.log(sub_mu) - sub_sigma2_lognormal/2
                sub_sigma_lognormal = sub_sigma2_lognormal**0.5
                kwargs['diameter'] = np.random.lognormal(sub_mu_lognormal, sub_sigma_lognormal, number)
                # check it worked
                # print('median = '+str(np.median(kwargs['diameter'])))
                # print('mean = '+str(np.mean(kwargs['diameter'])))
                # print('sigma = '+str(np.std(kwargs['diameter'])))
            else:
                raise Exception("no valid initial subsea droplet size distribution specified")

        if 'oiltype' in kwargs:
            logger.warning(
                'Seed argument *oiltype* is deprecated, use *oil_type* instead'
            )
            kwargs['oil_type'] = kwargs['oiltype']
            del kwargs['oiltype']

        if 'oil_type' in kwargs:
            if self.get_config('seed:oil_type') != kwargs['oil_type']:
                self.__set_seed_config__('seed:oil_type', kwargs['oil_type'])
            del kwargs['oil_type']
        else:
            logger.info('Oil type not specified, using default: ' +
                        self.get_config('seed:oil_type'))
        self.set_oiltype(self.get_config('seed:oil_type'))

        if self.oil_weathering_model == 'noaa':
            self.Density = Density(self.oiltype.oil)
            self.KinematicViscosity = KinematicViscosity(self.oiltype.oil)
            oil_density = self.Density.at_temp(285)
            oil_viscosity = self.KinematicViscosity.at_temp(285)
            logger.info(
                'Using density %s and viscosity %s of oiltype %s' %
                (oil_density, oil_viscosity, self.get_config('seed:oil_type')))
            kwargs['density'] = oil_density
            kwargs['viscosity'] = oil_viscosity

        if 'm3_per_hour' in kwargs:
            m3_per_hour = kwargs['m3_per_hour']
            del kwargs['m3_per_hour']
        else:
            m3_per_hour = self.get_config('seed:m3_per_hour')

        if 'number' in kwargs:
            num_elements = kwargs['number']
        else:
            num_elements = self.get_config('seed:number')
        time = kwargs['time']
        if type(time) is list:
            duration_hours = ((time[1] - time[0]).total_seconds()) / 3600
            if duration_hours == 0:
                duration_hours = 1.
        else:
            duration_hours = 1.  # For instantaneous spill, we use 1h
        kwargs['mass_oil'] = (m3_per_hour * duration_hours / num_elements *
                              kwargs['density'])

        self.store_oil_seed_metadata(**kwargs)

        super(OpenOil, self).seed_elements(*args, **kwargs)

    def seed_cone(self, *args, **kwargs):
        if 'oiltype' in kwargs:
            logger.warning(
                'Seed argument *oiltype* is deprecated, use *oil_type* instead'
            )
            kwargs['oil_type'] = kwargs['oiltype']
            del kwargs['oiltype']

        self.store_oil_seed_metadata(**kwargs)

        super(OpenOil, self).seed_cone(*args, **kwargs)

    def seed_from_gml(self, gmlfile, num_elements=1000, *args, **kwargs):
        """Read oil slick contours from GML file, and seed particles within."""

        # Specific imports
        import datetime
        from matplotlib.path import Path
        from xml.etree import ElementTree
        from matplotlib.patches import Polygon
        import pyproj

        namespaces = {
            'od': 'http://cweb.ksat.no/cweb/schema/geoweb/oil',
            'gml': 'http://www.opengis.net/gml'
        }
        slicks = []

        with open(gmlfile, 'rt') as e:
            tree = ElementTree.parse(e)

        pos1 = 'od:oilDetectionMember/od:oilDetection/od:oilSpill/gml:Polygon'
        pos2 = 'gml:exterior/gml:LinearRing/gml:posList'

        # This retrieves some other types of patches, found in some files only
        # Should be combines with the above, to get all patches
        #pos1 = 'od:oilDetectionMember/od:oilDetection/od:oilSpill/
        # gml:Surface/gml:polygonPatches'
        #pos2 = 'gml:PolygonPatch/gml:exterior/gml:LinearRing/gml:posList'

        # Find detection time
        time_pos = 'od:oilDetectionMember/od:oilDetection/od:detectionTime'
        oil_time = datetime.datetime.strptime(
            tree.find(time_pos, namespaces).text, '%Y-%m-%dT%H:%M:%S.%fZ')

        for patch in tree.findall(pos1, namespaces):
            pos = patch.find(pos2, namespaces).text
            c = np.array(pos.split()).astype(float)
            lon = c[0::2]
            lat = c[1::2]
            slicks.append(Polygon(list(zip(lon, lat))))

        # Find boundary and area of all patches
        lons = np.array([])
        lats = lons.copy()
        for slick in slicks:
            ext = slick.get_extents()
            lons = np.append(lons, [ext.xmin, ext.xmax])
            lats = np.append(lats, [ext.ymin, ext.ymax])
            # Make a stereographic projection centred on the polygon
        lonmin = lons.min()
        lonmax = lons.max()
        latmin = lats.min()
        latmax = lats.max()

        # Place n points within the polygons
        proj = pyproj.Proj(
            '+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f' %
            (latmin, latmax, (latmin + latmax) / 2, (lonmin + lonmax) / 2))
        slickarea = np.array([])
        for slick in slicks:
            lonlat = slick.get_xy()
            lon = lonlat[:, 0]
            lat = lonlat[:, 1]
            x, y = proj(lon, lat)

            area_of_polygon = 0.0
            for i in range(-1, len(x) - 1):
                area_of_polygon += x[i] * (y[i + 1] - y[i - 1])
            area_of_polygon = abs(area_of_polygon) / 2.0
            slickarea = np.append(slickarea, area_of_polygon)  # in m2

        # Make points
        deltax = np.sqrt(np.sum(slickarea) / num_elements)

        lonpoints = np.array([])
        latpoints = np.array([])
        for i, slick in enumerate(slicks):
            lonlat = slick.get_xy()
            lon = lonlat[:, 0]
            lat = lonlat[:, 1]
            x, y = proj(lon, lat)
            xvec = np.arange(x.min(), x.max(), deltax)
            yvec = np.arange(y.min(), y.max(), deltax)
            x, y = np.meshgrid(xvec, yvec)
            lon, lat = proj(x, y, inverse=True)
            lon = lon.ravel()
            lat = lat.ravel()
            points = np.c_[lon, lat]
            ind = Path(slick.xy).contains_points(points)
            lonpoints = np.append(lonpoints, lon[ind])
            latpoints = np.append(latpoints, lat[ind])

        # Finally seed at found positions
        kwargs['lon'] = lonpoints
        kwargs['lat'] = latpoints
        self.seed_elements(time=oil_time, **kwargs)

    def seed_from_geotiff_thickness(self,
                                    filename,
                                    number=50000,
                                    *args,
                                    **kwargs):
        '''Seed from files as provided by Prof. Chuanmin Hu'''

        from osgeo import gdal, ogr, osr

        if not 'time' in kwargs:
            try:  # get time from filename
                time = datetime.strptime(filename[-28:-13], '%Y%m%d.%H%M%S')
                logger.info('Parsed time from filename: %s' % time)
            except:
                time = datetime.now()
                logger.warning('Could not pase time from filename, '
                               'using present time: %s' % time)

        ds = gdal.Open(filename)

        srcband = ds.GetRasterBand(1)
        data = srcband.ReadAsArray()

        thickness_microns = [0.04, 0.44, 4.4, 16]  # NB: approximate
        categories = [1, 2, 3, 4]  # categories

        # Make memory raster bands for each category
        memrastername = filename + '.mem'
        memdriver = gdal.GetDriverByName('MEM')
        mem_ds = memdriver.CreateCopy(memrastername, ds)
        for cat in categories:
            mem_ds.AddBand(gdal.GDT_Byte)
            mem_band = mem_ds.GetRasterBand(cat)
            mem_band.WriteArray(data == cat)

        # Make memory polygons for each category
        drv = ogr.GetDriverByName('MEMORY')
        mem_vector_ds = [0] * len(categories)
        mem_vector_layers = [0] * len(categories)
        for cat in categories:
            memshapename = filename + '%i.mem' % cat
            mem_vector_ds[cat - 1] = drv.CreateDataSource(memshapename)
            mem_vector_layers[cat-1] = \
                mem_vector_ds[cat-1].CreateLayer(
                    'thickness%i' % cat, srs=None)
            gdal.Polygonize(mem_ds.GetRasterBand(cat),
                            None,
                            mem_vector_layers[cat - 1],
                            -1, [],
                            callback=None)

        total_area = np.zeros(len(categories))
        layers = [0] * len(categories)

        src_srs = osr.SpatialReference()
        src_srs.ImportFromEPSG(4269)
        tgt_srs = osr.SpatialReference()
        tgt_srs.ImportFromEPSG(3857)
        transform = osr.CoordinateTransformation(src_srs, tgt_srs)
        for cat in categories:
            memshapename = filename + '%i.shp' % cat
            layers[cat - 1] = mem_vector_layers[cat - 1]
            areas = np.zeros(layers[cat - 1].GetFeatureCount())
            for i, feature in enumerate(layers[cat - 1]):
                geom = feature.GetGeometryRef()
                geom.Transform(transform)  # To get area in m2
                areas[i] = geom.GetArea()
            # Delete largest polygon, which is outer border
            outer = np.where(areas == max(areas))[0]
            areas[outer] = 0
            total_area[cat - 1] = np.sum(areas)
            layers[cat - 1].DeleteFeature(outer)
            layers[cat - 1].ResetReading()

        # Calculate how many elements to be seeded for each category
        areas_weighted = total_area * thickness_microns
        numbers = number * areas_weighted / np.sum(areas_weighted)
        numbers = np.round(numbers).astype(int)
        oil_density = 1000
        mass_oil = (total_area * thickness_microns / 1e6) * oil_density

        for i, num in enumerate(numbers):
            self.seed_from_shapefile([mem_vector_layers[i]],
                                     oil_film_thickness=thickness_microns[i] /
                                     1000000.,
                                     mass_oil=mass_oil[i] / num,
                                     number=num,
                                     time=time,
                                     *args,
                                     **kwargs)

    def _substance_name(self):
        return self.get_oil_name()
