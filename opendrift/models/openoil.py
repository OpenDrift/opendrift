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

from io import open
import os
import numpy as np
from datetime import datetime
import logging
import matplotlib.pyplot as plt

from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray
import opendrift.models.noaa_oil_weathering as noaa

try:
    from itertools import izip as zip
except ImportError:
    pass

# Defining the oil element properties
class Oil(LagrangianArray):
    """Extending LagrangianArray with variables relevant for oil particles."""

    variables = LagrangianArray.add_variables([
        ('mass_oil', {'dtype': np.float32,
                      'units': 'kg',
                      'default': 1}),
        ('viscosity', {'dtype': np.float32,
                       #'unit': 'mm2/s (centiStokes)',
                       'units': 'N s/m2 (Pa s)',
                       'default': 0.005}),
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 880}),
        ('wind_drift_factor', {'dtype': np.float32,
                               'units': '%',
                               'default': 0.03}),
        ('age_exposure_seconds', {'dtype': np.float32,
                                  'units': 's',
                                  'default': 0}),
        ('age_emulsion_seconds', {'dtype': np.float32,
                                  'units': 's',
                                  'default': 0}),
        ('bulltime', {'dtype': np.float32,
                      'units': 's',
                      'default': 0}),
        ('interfacial_area', {'dtype': np.float32,
                              'units': 'm2',
                              'default': 0}),
        ('mass_dispersed', {'dtype': np.float32,
                            'units': 'kg',
                            'default': 0}),
        ('mass_evaporated', {'dtype': np.float32,
                             'units': 'kg',
                             'default': 0}),
        ('fraction_evaporated', {'dtype': np.float32,
                                 'units': '%',
                                 'default': 0}),
        ('water_fraction', {'dtype': np.float32,
                            'units': '%',
                            'default': 0}),
        ('oil_film_thickness', {'dtype': np.float32,
                                'units': 'm',
                                'default': 0.001})])


class OpenOil(OpenDriftSimulation):
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_stokes_drift_x_velocity',
                          'sea_surface_wave_stokes_drift_y_velocity',
                          'sea_ice_area_fraction',
                          'sea_water_temperature',
                          'sea_floor_depth_below_sea_level',
                          'x_wind', 'y_wind', 'land_binary_mask']

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_stokes_drift_x_velocity': 0,
                       'sea_surface_wave_stokes_drift_y_velocity': 0,
                       'sea_floor_depth_below_sea_level': 0,
                       'sea_ice_area_fraction': 0,
                       'sea_water_temperature': 12,
                       'x_wind': 0, 'y_wind': 0}

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'missing_data': 'gray', 'stranded': 'red',
                     'evaporated': 'yellow', 'dispersed': 'magenta'}

    configspec = '''
        [processes]
            dispersion = boolean(default=True)
            evaporation = boolean(default=True)
            emulsification = boolean(default=True)
    '''

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

    # Workaround as ADIOS oil library uses
    # max water fraction of 0.9 for all crude oils
    max_water_fraction  = {
        'MARINE GAS OIL 500 ppm S 2017': 0.1}


    def __init__(self, weathering_model='default', *args, **kwargs):

        self._add_configstring(self.configspec)

        if weathering_model == 'noaa':
            try:
                from oil_library import _get_db_session
                from oil_library.models import Oil, ImportedRecord
            except:
                raise ImportError(
                    'NOAA oil library must be installed from: '
                    'https://github.com/NOAA-ORR-ERD/OilLibrary')
            # Get list of all oiltypes in NOAA database
            session = _get_db_session()
            if 'location' in kwargs:
                self.oiltypes = session.query(Oil.name).join(
                                ImportedRecord).filter(ImportedRecord.
                                location==kwargs['location']).all()
                del kwargs['location']
                all_oiltypes = session.query(Oil.name).all()
                generic_oiltypes = [o for o in all_oiltypes if o[0][0:2] == '*G']
                self.oiltypes.extend(generic_oiltypes)
            else:
                self.oiltypes = session.query(Oil.name).all()
            self.oiltypes = sorted([o[0] for o in self.oiltypes])
            self.oiltypes = [ot for ot in self.oiltypes if ot not in self.duplicate_oils]
        else:
            # Read oil properties from file
            self.oiltype_file = os.path.dirname(os.path.realpath(__file__)) + \
                '/oilprop.dat'
            oilprop = open(self.oiltype_file, 'r', encoding='utf-8')
            oiltypes = []
            linenumbers = []
            for i, line in enumerate(oilprop.readlines()):
                if line[0].isalpha():
                    oiltype = line.strip()[:-2].strip()
                    oiltypes.append(oiltype)
                    linenumbers.append(i)
            oiltypes, linenumbers = zip(*sorted(zip(oiltypes, linenumbers)))
            self.oiltypes = oiltypes
            self.oiltypes_linenumbers = linenumbers

        self.oil_weathering_model = weathering_model

        # Update config with oiltypes 
        oiltypes = [str(a) for a in self.oiltypes]
        self._add_config('seed:oil_type', oiltypes,
                         'Oil type', overwrite=True)

        # Calling general constructor of parent class
        super(OpenOil, self).__init__(*args, **kwargs)

    def evaporate(self):
        if self.get_config('processes:evaporation') is True:
            logging.debug('   Calculating: evaporation')
            windspeed = np.sqrt(self.environment.x_wind**2 +
                                self.environment.y_wind**2)

            # Store evaporation fraction at beginning of timestep
            fraction_evaporated_previous = self.elements.fraction_evaporated

            # Evaporate only elements at surface
            at_surface = (self.elements.z == 0)
            if np.isscalar(at_surface):
                at_surface = at_surface*np.ones(self.num_elements_active(),
                                                dtype=bool)
            Urel = windspeed/self.oil_data['reference_wind']  # Relative wind
            h = 2  # Film thickness in mm, harcoded for now

            # Calculate exposure time
            #   presently without compensation for sea temperature
            delta_exposure_seconds = \
                (self.oil_data['reference_thickness']/h)*Urel * \
                self.time_step.total_seconds()
            if np.isscalar(self.elements.age_exposure_seconds):
                self.elements.age_exposure_seconds += delta_exposure_seconds
            else:
                self.elements.age_exposure_seconds[at_surface] += \
                    delta_exposure_seconds[at_surface]

            self.elements.fraction_evaporated = np.interp(
                self.elements.age_exposure_seconds,
                self.oil_data['tref'], self.oil_data['fref'])
            self.mass_evaporated = \
                self.elements.mass_oil*self.elements.fraction_evaporated
            # Remove evaporated part from mass_oil
            mass_evaporated_timestep = self.elements.mass_oil*(
                self.elements.fraction_evaporated -
                fraction_evaporated_previous)
            self.elements.mass_oil -= mass_evaporated_timestep
            self.elements.mass_evaporated += mass_evaporated_timestep

            # Evaporation probability equals the difference in fraction
            # evaporated at this timestep compared to previous timestep,
            # divided by the remaining fraction of the particle at
            # previous timestep
            evaporation_probability = ((self.elements.fraction_evaporated -
                                        fraction_evaporated_previous) /
                                       (1 - fraction_evaporated_previous))
            evaporation_probability[~at_surface] = 0
            evaporated_indices = (evaporation_probability >
                                  np.random.rand(self.num_elements_active(),))
            self.deactivate_elements(evaporated_indices, reason='evaporated')

    def emulsification(self):
        if self.get_config('processes:emulsification') is True:
            logging.debug('   Calculating: emulsification')
            if 'reference_wind' not in self.oil_data.keys():
                logging.debug('Emulsification is currently only possible when'
                              'importing oil properties from file.')
            else:
                windspeed = np.sqrt(self.environment.x_wind**2 +
                                    self.environment.y_wind**2)
                # Apparent emulsion age of particles
                Urel = windspeed/self.oil_data['reference_wind']  # Relative wind
                self.elements.age_emulsion_seconds += \
                    Urel*self.time_step.total_seconds()

                self.elements.water_fraction = np.interp(
                    self.elements.age_emulsion_seconds,
                    self.oil_data['tref'], self.oil_data['wmax'])

    def disperse(self):
        if self.get_config('processes:dispersion') is True:

            logging.debug('   Calculating: dispersion')
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
            C_Roy = 2400.0 * np.exp(-73.682*np.sqrt(
                self.elements.viscosity/self.elements.density))

            q_disp = C_Roy * c_disp * v_entrain / self.elements.density

            oil_mass_loss = (q_disp * self.time_step.total_seconds() *
                             self.elements.density)*self.elements.mass_oil

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
            logging.debug('Skipping oil weathering for backwards run')
            return
        self.timer_start('main loop:updating elements:oil weathering')
        if self.oil_weathering_model == 'noaa':
            self.oil_weathering_noaa()
        else:
            self.oil_weathering_default()
        self.timer_end('main loop:updating elements:oil weathering')

    def oil_weathering_default(self):

        logging.debug('Default oil weathering')

        ## Evaporation
        self.evaporate()

        # Emulsification
        self.emulsification()

        # Dispersion
        self.disperse()

    def prepare_run(self):

        if self.oil_weathering_model == 'noaa':
            self.noaa_mass_balance = {}
            # Populate with seeded mass spread on oiltype.mass_fraction
            mass_oil = np.atleast_1d(self.elements_scheduled.mass_oil)
            if len(mass_oil) == 1:
                mass_oil = mass_oil*np.ones(self.num_elements_total())
            self.noaa_mass_balance['mass_components'] = \
                np.asarray(self.oiltype.mass_fraction)*(mass_oil.reshape(
                    (self.num_elements_total(), -1)))
            self.noaa_mass_balance['mass_evaporated'] = \
                self.noaa_mass_balance['mass_components']*0

            self.oil_water_interfacial_tension = \
                self.oiltype.oil_water_surface_tension()[0]
            logging.info('Oil-water surface tension is %f Nm' % 
                         self.oil_water_interfacial_tension)
        else:
            logging.info('Using default oil-water tension of 0.03Nm')
            self.oil_water_interfacial_tension = 0.03

    def oil_weathering_noaa(self):
        '''Oil weathering scheme adopted from NOAA PyGNOME model:
        https://github.com/NOAA-ORR-ERD/PyGnome
        '''
        logging.debug('NOAA oil weathering')
        # C to K
        self.environment.sea_water_temperature[
            self.environment.sea_water_temperature < 100] += 273.15

        #########################################################
        # Update density and viscosity according to temperature
        #########################################################
        try:  # New version of OilLibrary
            self.timer_start('main loop:updating elements:oil weathering:updating viscosities')
            oil_viscosity = self.oiltype.kvis_at_temp(
                self.environment.sea_water_temperature)
            self.timer_end('main loop:updating elements:oil weathering:updating viscosities')
            self.timer_start('main loop:updating elements:oil weathering:updating densities')
            oil_density = self.oiltype.density_at_temp(
                self.environment.sea_water_temperature)
            self.timer_end('main loop:updating elements:oil weathering:updating densities')
        except:  # Old version of OilLibrary
            oil_viscosity = np.array(
                [self.oiltype.get_viscosity(t) for t in
                 self.environment.sea_water_temperature])
            oil_density = np.array(
                [self.oiltype.get_density(t) for t in
                 self.environment.sea_water_temperature])

        # Calculate emulsion density
        self.elements.density = (
            self.elements.water_fraction*self.sea_water_density() +
           (1 - self.elements.water_fraction) * oil_density)

        # Calculate emulsion viscosity
        visc_f_ref = 0.84  # From PyGNOME
        visc_curvfit_param = 1.5e3 # units are sec^0.5 / m
        fw_d_fref = self.elements.water_fraction/visc_f_ref
        kv1 = np.sqrt(oil_viscosity)*visc_curvfit_param
        kv1[kv1<1] = 1
        kv1[kv1>10] = 10
        self.elements.fraction_evaporated = self.elements.mass_evaporated/(self.elements.mass_oil+self.elements.mass_evaporated)
        self.elements.viscosity = (
            oil_viscosity*np.exp(kv1*self.elements.fraction_evaporated)*
                (1 + (fw_d_fref / (1.187 - fw_d_fref))) ** 2.49)


        if self.get_config('processes:evaporation') is True:
            self.timer_start('main loop:updating elements:oil weathering:evaporation')
            self.evaporation_noaa()
            self.timer_end('main loop:updating elements:oil weathering:evaporation')

        if self.get_config('processes:emulsification') is True:
            self.timer_start('main loop:updating elements:oil weathering:emulsification')
            self.emulsification_noaa()
            self.timer_end('main loop:updating elements:oil weathering:emulsification')

        if self.get_config('processes:dispersion') is True:
            self.timer_start('main loop:updating elements:oil weathering:dispersion')
            self.disperse_noaa()
            self.timer_end('main loop:updating elements:oil weathering:dispersion')

    def disperse_noaa(self):
        logging.debug('    Calculating: dispersion - NOAA')
        # From NOAA PyGnome model:
        # https://github.com/NOAA-ORR-ERD/PyGnome/
        c_disp = np.power(self.wave_energy_dissipation(), 0.57) * \
            self.sea_surface_wave_breaking_fraction()
        # Roy's constant
        C_Roy = 2400.0 * np.exp(-73.682*np.sqrt(
            self.elements.viscosity/self.elements.density))
        v_entrain = 3.9E-8
        q_disp = C_Roy * c_disp * v_entrain / self.elements.density
        fraction_dispersed = (q_disp * self.time_step.total_seconds() *
                         self.elements.density)
        oil_mass_loss = fraction_dispersed*self.elements.mass_oil

        self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :] = \
            self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :]*(1-fraction_dispersed[:, np.newaxis])

        self.elements.mass_oil -= oil_mass_loss
        self.elements.mass_dispersed += oil_mass_loss

    def plot_droplet_spectrum(self):
        '''Plotting distribution of droplet radii, for debugging'''
        plt.hist(self.elements.diameter/2.0)
        plt.show()

    def evaporation_noaa(self):
        #############################################
        # Evaporation, for elements at surface only
        #############################################
        logging.debug('    Calculating evaporation - NOAA')
        surface = np.where(self.elements.z == 0)[0]  # of active elements
        if len(surface) == 0:
            logging.debug('All elements submerged, no evaporation')
            return
        surfaceID = self.elements.ID[surface] - 1    # of any elements
        # Area for each element, repeated for each component
        volume = (self.elements.mass_oil[surface] /
                  self.elements.density[surface])
        area = volume/self.elements.oil_film_thickness[surface]
        evap_decay_constant = noaa.evap_decay_constant(
            self.oiltype, self.wind_speed()[surface],
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
        # Emulsification (surface only?)
        #############################################
        logging.debug('    Calculating emulsification - NOAA')
        emul_time = self.oiltype.bulltime
        emul_constant = self.oiltype.bullwinkle
        # max water content fraction - get from database
        Y_max = self.oiltype.get('emulsion_water_fraction_max')
        if self.oil_name in self.max_water_fraction:
            max_water_fraction = self.max_water_fraction[self.oil_name]
            logging.debug('Overriding max water fraxtion with value %f instead of default %f'
                          % (max_water_fraction, Y_max))
            Y_max = max_water_fraction
        # emulsion
        if Y_max <= 0:
            logging.debug('Oil does not emulsify, returning.')
            return
        # Constants for droplets
        drop_min = 1.0e-6
        drop_max = 1.0e-5
        S_max = (6. / drop_min) * (Y_max / (1.0 - Y_max))
        S_min = (6. / drop_max) * (Y_max / (1.0 - Y_max))
        # Emulsify...
        fraction_evaporated = self.elements.mass_evaporated / (
            self.elements.mass_evaporated + self.elements.mass_oil)
        # f ((le_age >= emul_time && emul_time >= 0.) || frac_evap[i] >= emul_C && emul_C > 0.)

        start_emulsion = np.where(
            ((self.elements.age_seconds >= emul_time) & (emul_time >= 0)) |
            ((fraction_evaporated >= emul_constant) & (emul_constant > 0))
            )[0]
        if len(start_emulsion) == 0:
            logging.debug('        Emulsification not yet started')
            return

        if self.oiltype.bulltime > 0:  # User has set value
            start_time = self.oiltype.bulltime*np.ones(len(start_emulsion))
        else:
            start_time = self.elements.age_seconds[start_emulsion]
            start_time[self.elements.age_emulsion_seconds[start_emulsion]
                       >= 0] = self.elements.bulltime[start_emulsion]
        # Update droplet interfacial area
        k_emul = noaa.water_uptake_coefficient(
                    self.oiltype, self.wind_speed()[start_emulsion])
        self.elements.interfacial_area[start_emulsion] = \
            self.elements.interfacial_area[start_emulsion] + \
            (k_emul*self.time_step.total_seconds()*
             np.exp((-k_emul/S_max)*(
                self.elements.age_seconds[start_emulsion] - start_time)))
        self.elements.interfacial_area[self.elements.interfacial_area >
                                       S_max] = S_max
        # Update water fraction
        self.elements.water_fraction[start_emulsion] = (
            self.elements.interfacial_area[start_emulsion]*drop_max/
            (6.0 + (self.elements.interfacial_area[start_emulsion]
             *drop_max)))
        self.elements.water_fraction[self.elements.interfacial_area >=
            ((6.0 / drop_max)*(Y_max/(1.0 - Y_max)))] = Y_max

    def advect_oil(self):
        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Wind drag for elements at ocean surface
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Deactivate elements hitting sea ice
        if hasattr(self.environment, 'sea_ice_area_fraction'):
            self.deactivate_elements(
                self.environment.sea_ice_area_fraction > 0.6,
                reason='oil-in-ice')

    def update(self):
        """Update positions and properties of oil particles."""

        # Oil weathering
        self.oil_weathering()

        # Horizontal advection
        self.advect_oil()

    def get_oil_budget(self):

        if self.time_step.days < 0:  # Backwards simulation
            return None

        z, dummy = self.get_property('z')
        mass_oil, status = self.get_property('mass_oil')

        if 'stranded' not in self.status_categories:
            self.status_categories.append('stranded')
        mass_submerged = np.ma.masked_where(
            ((status == self.status_categories.index('stranded')) |
            (z == 0.0)), mass_oil)
        mass_submerged = np.ma.sum(mass_submerged, axis=1).filled(0)

        mass_surface = np.ma.masked_where(
            ((status == self.status_categories.index('stranded')) |
            (z < 0.0)), mass_oil)
        mass_surface = np.ma.sum(mass_surface, axis=1).filled(0)

        mass_stranded = np.ma.sum(np.ma.masked_where(
            status != self.status_categories.index('stranded'),
            mass_oil), axis=1).filled(0)
        mass_evaporated, status = self.get_property('mass_evaporated')
        mass_evaporated = np.sum(mass_evaporated, axis=1).filled(0)
        mass_dispersed, status = self.get_property('mass_dispersed')
        mass_dispersed = np.sum(mass_dispersed, axis=1).filled(0)

        oil_budget = {
            'mass_dispersed': mass_dispersed,
            'mass_submerged': mass_submerged,
            'mass_surface': mass_surface,
            'mass_stranded': mass_stranded,
            'mass_evaporated': mass_evaporated,
            'mass_total': (mass_dispersed + mass_submerged +
                           mass_surface + mass_stranded + mass_evaporated)
            }

        return oil_budget

    def plot_oil_budget(self, filename=None, ax=None):

        if ax==None:
            plt.close()

        if self.time_step.days < 0:  # Backwards simulation
            fig = plt.figure(figsize=(10, 6.))
            plt.text(0.1, 0.5, 'Oil weathering deactivated for '
                                'backwards simulations')
            plt.axis('off')
            if filename is not None:
                plt.savefig(filename)
                plt.close()
            plt.show()
            return

        b = self.get_oil_budget()

        oil_budget = np.row_stack(
            (b['mass_dispersed'], b['mass_submerged'],
             b['mass_surface'], b['mass_stranded'], b['mass_evaporated']))

        budget = np.cumsum(oil_budget, axis=0)

        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds()/3600. for t in time_relative])

        if ax==None:
            fig = plt.figure(figsize=(10, 6.))  # Suitable aspect ratio
            # Left axis showing oil mass
            ax1 = fig.add_subplot(111)
        else:
            ax1 = ax

        # Hack: make some emply plots since fill_between does not support label
        if np.sum(b['mass_dispersed']) > 0:
            ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                          color='darkslategrey', label='dispersed'))
            ax1.fill_between(time, 0, budget[0, :], facecolor='darkslategrey')
        if np.sum(b['mass_submerged']) > 0:
            ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                          color='darkblue', label='submerged'))
            ax1.fill_between(time, budget[0, :], budget[1, :],
                             facecolor='darkblue')
        if np.sum(b['mass_surface']) > 0:
            ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                          color='royalblue', label='surface'))
            ax1.fill_between(time, budget[1, :], budget[2, :],
                             facecolor='royalblue')
        if np.sum(b['mass_stranded']) > 0:
            ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                          color='black', label='stranded'))
            ax1.fill_between(time, budget[2, :], budget[3, :],
                             facecolor='black')
        if np.sum(b['mass_evaporated']) > 0:
            ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                          color='skyblue', label='evaporated'))
            ax1.fill_between(time, budget[3, :], budget[4, :],
                             facecolor='skyblue')

        ax1.set_ylim([0, budget.max()])
        ax1.set_xlim([0, time.max()])
        ax1.set_ylabel('Mass oil  [%s]' %
                       self.elements.variables['mass_oil']['units'])
        ax1.set_xlabel('Time  [hours]')
        # Right axis showing percent
        ax2 = ax1.twinx()
        ax2.set_ylim([0, 100])
        ax2.set_ylabel('Percent')
        if not hasattr(self, 'oil_name'):  # TODO
            self.oil_name = 'unknown oiltype'
            # TODO line below is dangerous when importing old files
            self.oil_name = self.get_config('seed:oil_type')
        plt.title('%s - %s to %s' %
                  (self.oil_name,
                   self.start_time.strftime('%Y-%m-%d %H:%M'),
                   self.time.strftime('%Y-%m-%d %H:%M')))
        # Shrink current axis's height by 10% on the bottom
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax1.legend(bbox_to_anchor=(0., -0.10, 1., -0.03), loc=1,
                   ncol=5, mode="expand", borderaxespad=0., fontsize=10)
        if filename is not None:
            plt.savefig(filename)
            plt.close()
        plt.show()

    def set_oiltype(self, oiltype):
        self.oil_name = oiltype
        self.set_config('seed:oil_type', oiltype)
        if self.oil_weathering_model == 'noaa':
            try:
                from oil_library import get_oil_props, _get_db_session
                from oil_library.models import Oil as ADIOS_Oil
                from oil_library.oil_props import OilProps
                session = _get_db_session()
                oils = session.query(ADIOS_Oil).filter(
                            ADIOS_Oil.name == oiltype)
                ADIOS_ids = [oil.adios_oil_id for oil in oils]
                if len(ADIOS_ids) == 0:
                    raise ValueError('Oil type "%s" not found in NOAA database' % oiltype)
                elif len(ADIOS_ids) == 1:
                    self.oiltype = get_oil_props(oiltype)
                else:
                    logging.warning('Two oils found with name %s (ADIOS IDs %s and %s). Using the first.' % (oiltype, ADIOS_ids[0], ADIOS_ids[1]))
                    self.oiltype = OilProps(oils[0])
            except Exception as e:
                logging.warning(e)
            return

        if oiltype not in self.oiltypes:
            raise ValueError('Oiltype %s is unknown. The following oiltypes are available: %s' % (oiltype, str(self.oiltypes)))
        indx = self.oiltypes.index(oiltype)
        linenumber = self.oiltypes_linenumbers[indx]
        oilfile = open(self.oiltype_file, 'r')
        for i in range(linenumber + 1):
            oilfile.readline()
        ref = oilfile.readline().split()
        self.oil_data = {}
        self.oil_data['reference_thickness'] = np.float(ref[0])
        self.oil_data['reference_wind'] = np.float(ref[1])
        tref = []
        fref = []
        wmax = []
        while True:
            line = oilfile.readline()
            if not line[0].isdigit():
                break
            line = line.split()
            tref.append(line[0])
            fref.append(line[1])
            wmax.append(line[3])
        self.oil_data['tref'] = np.array(tref, dtype='float')*3600.
        self.oil_data['fref'] = np.array(fref, dtype='float')*.01
        self.oil_data['wmax'] = np.array(wmax, dtype='float')
        self.oil_data['oiltype'] = oiltype  # Store name of oil type

    def seed_elements(self, *args, **kwargs):
        if 'oiltype' in kwargs:
            self.set_config('seed:oil_type', kwargs['oiltype'])
            del kwargs['oiltype']
        else:
            logging.info('Oil type not specified, using default: ' +
                         self.get_config('seed:oil_type'))
        self.set_oiltype(self.get_config('seed:oil_type'))

        if self.oil_weathering_model == 'noaa':
            try:  # Older version of OilLibrary
                oil_density = self.oiltype.get_density(285)  # 12 degrees
                oil_viscosity = self.oiltype.get_viscosity(285)
            except:  # Newer version of OilLibrary
                oil_density = self.oiltype.density_at_temp(285)
                oil_viscosity = self.oiltype.kvis_at_temp(285)
            logging.info('Using density %s and viscosity %s of oiltype %s' %
                         (oil_density, oil_viscosity, self.get_config('seed:oil_type')))
            kwargs['density'] = oil_density
            kwargs['viscosity'] = oil_viscosity

        if 'm3_per_hour' in kwargs:
            # From given volume rate, we calculate the mass per element
            num_elements = kwargs['number']
            time = kwargs['time']
            if type(time) is list:
                duration_hours = ((time[1] - time[0]).total_seconds())/3600
            else:
                duration_hours = 1  # For instantaneous spill, we use 1h
            kwargs['mass_oil'] = (kwargs['m3_per_hour']*duration_hours/
                                  num_elements*kwargs['density'])
            del kwargs['m3_per_hour']

        super(OpenOil, self).seed_elements(*args, **kwargs)

    def seed_from_gml(self, gmlfile, num_elements=1000, *args, **kwargs):
        """Read oil slick contours from GML file, and seed particles within."""

        # Specific imports
        import datetime
        # nxutils: Deprecated since version 1.2.0:
        # Use contains_points() instead.
        have_nx = True
        try:
            import matplotlib.nxutils as nx
        except:
            have_nx = False
            from matplotlib.path import Path
        from xml.etree import ElementTree
        from matplotlib.patches import Polygon
        from mpl_toolkits.basemap import pyproj

        namespaces = {'od': 'http://cweb.ksat.no/cweb/schema/geoweb/oil',
                      'gml': 'http://www.opengis.net/gml'}
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
            c = np.array(pos.split()).astype(np.float)
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
        proj = pyproj.Proj('+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f'
                           % (latmin, latmax,
                              (latmin+latmax)/2, (lonmin+lonmax)/2))
        slickarea = np.array([])
        for slick in slicks:
            lonlat = slick.get_xy()
            lon = lonlat[:, 0]
            lat = lonlat[:, 1]
            x, y = proj(lon, lat)

            area_of_polygon = 0.0
            for i in range(-1, len(x)-1):
                area_of_polygon += x[i] * (y[i+1] - y[i-1])
            area_of_polygon = abs(area_of_polygon) / 2.0
            slickarea = np.append(slickarea, area_of_polygon)  # in m2

        # Make points
        deltax = np.sqrt(np.sum(slickarea)/num_elements)

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
            if have_nx:
                ind = nx.points_inside_poly(points, slick.xy)
            else:
                ind = Path(slick.xy).contains_points(points)
            lonpoints = np.append(lonpoints, lon[ind])
            latpoints = np.append(latpoints, lat[ind])

        # Finally seed at found positions
        kwargs['lon'] = lonpoints
        kwargs['lat'] = latpoints
        self.seed_elements(time=oil_time, **kwargs)

    def seed_from_geotiff_thickness(self, filename, number=50000,
                                    *args, **kwargs):
        '''Seed from files as provided by Prof. Chuanmin Hu'''

        import gdal
        import ogr

        if not 'time' is kwargs:
            try:  # get time from filename
                timestr = filename[-28:-13]
                time = datetime.strptime(
                        filename[-28:-13], '%Y%m%d.%H%M%S')
                logging.info('Parsed time from filename: %s' % time)
            except:
                time = datetime.now()
                logging.warning('Could not pase time from filename, '
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
            mem_band.WriteArray(data==cat)

        # Make memory polygons for each category
        drv = ogr.GetDriverByName('MEMORY')
        mem_vector_ds = [0]*len(categories)
        mem_vector_layers = [0]*len(categories)
        for cat in categories:
            memshapename = filename + '%i.mem' % cat
            mem_vector_ds[cat-1] = drv.CreateDataSource(memshapename)
            mem_vector_layers[cat-1] = \
                mem_vector_ds[cat-1].CreateLayer(
                    'thickness%i' % cat, srs=None)
            gdal.Polygonize(mem_ds.GetRasterBand(cat), None,
                            mem_vector_layers[cat-1],
                            -1, [], callback=None)

        total_area = np.zeros(len(categories))
        layers = [0]*len(categories)
        for cat in categories:
            memshapename = filename + '%i.shp' % cat
            layers[cat-1] = mem_vector_layers[cat-1]
            areas = np.zeros(layers[cat-1].GetFeatureCount())
            for i, feature in enumerate(layers[cat-1]):
                areas[i] = feature.GetGeometryRef().GetArea()
            # Delete largest polygon, which is outer border
            outer = np.where(areas==max(areas))[0]
            areas[outer] = 0
            total_area[cat-1] = np.sum(areas)
            layers[cat-1].DeleteFeature(outer)
            layers[cat-1].ResetReading()
            
        # Calculate how many elements to be seeded for each category
        areas_weighted = total_area*thickness_microns
        numbers = number*areas_weighted/np.sum(areas_weighted)
        numbers = np.round(numbers).astype(int)

        for i, num in enumerate(numbers):
            self.seed_from_shapefile([mem_vector_layers[i]],
                oil_film_thickness=thickness_microns[i]/1000000.,
                number=num, time=time, *args, **kwargs)
