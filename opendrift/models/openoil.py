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

"""
OpenOil is an oil drift module bundled within the OpenDrift framework. There is both a 2D-version, and 3D-version which adds turbulent vertical mixing.

The oil weathering calculations is based on the NOAA-ERR-ERD OilLibrary package, which must be installed as a dependency. The code for evaporation and emulsification in OpenOil is borrowed from the NOAA PyGnome code, and adapted to the OpenDrift architecture.

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
The vertical (and thus indirectly also the horisontal) motion of oil (droplets) is affected by oil density and droplet diameters.

When using the NOAA oil weathering model (``o = OpenOil3D(weathering_model='noaa')``), the density is obtained from the NOAA database according to the oiltype selected when seeding. This value can not be overridden by the user, and it will also change during the simulation due to oil weathering processes (evaporation and emulsification). When using the default (primitive) weathering model (o = OpenOil3D()), the parameter 'density' may be set when seeding (default is 880 kg/m3).

The droplet diameter may be given explicitly when seeding, e.g.::

    o.seed_elements(4, 60, number=100, time=datetime.now(), diameter=1e-5)

In this case, the diameter will not change during the simulation, which is useful e.g. for sensitivity tests. The same diameter will be used for all elements for this example, but an array of the same length as the number of elements may also be provided.

If a constant droplet diameter is not given by the user, it will be chosen randomly within given config limits for a subsea spill ('blowout'), and modified after any later wave breaking event. Oil droplets seeded under sea surface (z<0) will be assigned initial diameters between the following limits, typical for a subsea blowout (Johansen, 2000)::

    o.config['input']['spill']['droplet_diameter_min_subsea'] = 0.0005  # 0.5 mm
    o.config['input']['spill']['droplet_diameter_max_subsea'] = 0.005   # 5 mm

After each wave breaking event, a new diameter will be chosen between configurable limits and exponent::

    o.config['turbulentmixing']['droplet_diameter_min_wavebreaking'] = 1e-5
    o.config['turbulentmixing']['droplet_diameter_max_wavebreaking'] = 1e-3
    o.config['turbulentmixing']['droplet_size_exponent] = 0

In this case a logarithmic distribution is used, with N ~ diameter^s, where s is the droplet_size_exponent. A droplet_size_exponent of 0 (default) gives uniform distribution. A value of -2.3 corresponds to the empirical number distribution as found by Delvigne and Sweeney, and 0.7 gives the corresponding volume distribution (number distribution multiplied by diameter3, since volume is proportinal to diameter3).

The droplet size limits could also have been calculated dynamically from oil viscosity, but this is not yet implemented. The limits above should eventually be adjusted before the time of seeding.
"""

from io import open
import os
import numpy as np
from datetime import datetime
import pyproj
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
        ('mass_biodegraded', {'dtype': np.float32,
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
            biodegradation = boolean(default=False)
            update_oilfilm_thickness = boolean(default=False)
        [drift]
            current_uncertainty = float(min=0, max=5, default=0.05)
            wind_uncertainty = float(min=0, max=5, default=.5)
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
            oiltypes = []
            linenumbers = []
            with open(self.oiltype_file) as f:
                for i, line in enumerate(f):
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

        # Overriding with specific configspec
        self._add_configstring(self.configspec)

    def update_surface_oilfilm_thickness(self):
        '''The mass of oil is summed within a grid of 20x20
        cells covering the oil at a given time. Each oil particle
        within each cell is given a film thickness as the amount of
        oil divided by the cell area.
        '''
        from scipy.stats import binned_statistic_2d
        surface = np.where(self.elements.z == 0)[0]
        if len(surface) == 0:
            print('No oil at surface, no film thickness to update')
            return
        print('Updating oil film thickness for %s of %s elements at surface' % (len(surface), self.num_elements_active()))
        meanlon = self.elements.lon[surface].mean()
        meanlat = self.elements.lat[surface].mean()
        # Using stereographic coordinates to get regular X and Y
        psproj = pyproj.Proj(
            '+proj=stere +lat_0=%s +lat_ts=%s +lon_0=%s' %
            (meanlat, meanlat, meanlon))
        X,Y = psproj(self.elements.lon[surface], self.elements.lat[surface])
        mass_bin, x_edge, y_edge, binnumber = binned_statistic_2d(
            X, Y, self.elements.mass_oil[surface],
            expand_binnumbers=True,
            statistic='sum', bins=100)
        bin_area = (x_edge[1]-x_edge[0])*(y_edge[1]-y_edge[0])
        oil_density = 1000  # ok approximation here
        film_thickness = (mass_bin/oil_density)/bin_area
        # Postulating min and max film thickness
        max_thickness = 0.01  # 1 cm
        min_thickness = 1e-9  # 1 nanometer
        if film_thickness.max() > max_thickness:
            print('Warning: decreasing thickness to %sm for %s of %s bins' % (max_thickness, np.sum(film_thickness>max_thickness), film_thickness.size))
            film_thickness[film_thickness>max_thickness] = max_thickness
        num_too_thin = np.sum((film_thickness<min_thickness) & (film_thickness>0))
        if num_too_thin > 0:
            print('Warning: increasing thickness to %sm for %s of %s bins' % (min_thickness, num_too_thin, film_thickness.size))
            film_thickness[film_thickness<min_thickness] = min_thickness

        # https://github.com/scipy/scipy/issues/7010
        binnumber = binnumber - 1

        bx = binnumber[0,:]
        by = binnumber[1,:]
        # Update thickness
        self.elements.oil_film_thickness[surface] = self.elements.oil_film_thickness[surface]*np.nan
        self.elements.oil_film_thickness[surface] = \
            film_thickness[bx, by]

    def evaporate(self):
        if self.get_config('processes:evaporation') is True:
            self.logger.debug('   Calculating: evaporation')
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
            self.logger.debug('   Calculating: emulsification')
            if 'reference_wind' not in self.oil_data.keys():
                self.logger.debug('Emulsification is currently only possible when'
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

    def biodegradation(self):
        if self.get_config('processes:biodegradation') is True:
            '''
            Oil biodegradation function based on the article:
            Adcroft et al. (2010), Simulations of underwater plumes of dissolved oil in the Gulf of Mexico.
            '''
            self.logger.debug('Calculating: biodegradation')

            swt = self.environment.sea_water_temperature.copy()
            swt[swt > 100] -= 273.15 # K to C
            age0 = self.time_step.total_seconds()/(3600*24)

            tau = (12)*(3**((20-swt)/10)) # Decay rate in days (temperature in Celsius)
            fraction_biodegraded = (1 - np.exp(-age0/tau))
            biodegraded_now = self.elements.mass_oil*fraction_biodegraded

            self.elements.mass_biodegraded = self.elements.mass_biodegraded + biodegraded_now
            self.elements.mass_oil = self.elements.mass_oil - biodegraded_now
            if self.oil_weathering_model == 'noaa':
                self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :] = \
                self.noaa_mass_balance['mass_components'][self.elements.ID - 1, :]*(1-fraction_biodegraded[:, np.newaxis])
            else:
                pass

    def disperse(self):
        if self.get_config('processes:dispersion') is True:

            self.logger.debug('   Calculating: dispersion')
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
            self.logger.debug('Skipping oil weathering for backwards run')
            return
        self.timer_start('main loop:updating elements:oil weathering')
        if self.oil_weathering_model == 'noaa':
            self.oil_weathering_noaa()
        else:
            self.oil_weathering_default()
        self.timer_end('main loop:updating elements:oil weathering')

    def oil_weathering_default(self):

        self.logger.debug('Default oil weathering')

        ## Evaporation
        self.evaporate()

        # Emulsification
        self.emulsification()

        # Dispersion
        self.disperse()

        # Biodegradation
        self.biodegradation()

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
            self.logger.info('Oil-water surface tension is %f Nm' %
                         self.oil_water_interfacial_tension)
        else:
            self.logger.info('Using default oil-water tension of 0.03Nm')
            self.oil_water_interfacial_tension = 0.03

    def oil_weathering_noaa(self):
        '''Oil weathering scheme adopted from NOAA PyGNOME model:
        https://github.com/NOAA-ORR-ERD/PyGnome
        '''
        self.logger.debug('NOAA oil weathering')
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

        if self.get_config('processes:biodegradation') is True:
            self.timer_start('main loop:updating elements:oil weathering:biodegradation')
            self.biodegradation()
            self.timer_end('main loop:updating elements:oil weathering:biodegradation')

    def disperse_noaa(self):
        self.logger.debug('    Calculating: dispersion - NOAA')
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
        self.logger.debug('    Calculating evaporation - NOAA')
        surface = np.where(self.elements.z == 0)[0]  # of active elements
        if len(surface) == 0:
            self.logger.debug('All elements submerged, no evaporation')
            return
        if self.elements.age_seconds[surface].min() > 3600*24:
            self.logger.debug('All surface oil elements older than 24 hours, ' +
                          'skipping further evaporation.')
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
        self.logger.debug('    Calculating emulsification - NOAA')
        emul_time = self.oiltype.bulltime
        emul_constant = self.oiltype.bullwinkle
        # max water content fraction - get from database
        Y_max = self.oiltype.get('emulsion_water_fraction_max')
        if self.oil_name in self.max_water_fraction:
            max_water_fraction = self.max_water_fraction[self.oil_name]
            self.logger.debug('Overriding max water fraxtion with value %f instead of default %f'
                          % (max_water_fraction, Y_max))
            Y_max = max_water_fraction
        # emulsion
        if Y_max <= 0:
            self.logger.debug('Oil does not emulsify, returning.')
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
            self.logger.debug('        Emulsification not yet started')
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
            mass_emulsion = mass_oil / (1 - water_fraction)
        I.e. water_fraction = 0 means pure oil, water_fraction = 0.5 means mixture of
        50% oil and 50% water, and water_fraction = 0.9 (which is maximum)
        means 10% oil and 90% water.
        """

        if self.time_step.days < 0:  # Backwards simulation
            return None

        z, dummy = self.get_property('z')
        mass_oil, status = self.get_property('mass_oil')
        density = self.get_property('density')[0][0,0]

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
        mass_biodegraded, status = self.get_property('mass_biodegraded')
        mass_biodegraded = np.sum(mass_biodegraded, axis=1).filled(0)

        oil_budget = {
            'oil_density': density,
            'mass_dispersed': mass_dispersed,
            'mass_submerged': mass_submerged,
            'mass_surface': mass_surface,
            'mass_stranded': mass_stranded,
            'mass_evaporated': mass_evaporated,
            'mass_biodegraded': mass_biodegraded,
            'mass_total': (mass_dispersed + mass_submerged +
                           mass_surface + mass_stranded + mass_evaporated + mass_biodegraded)
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
             b['mass_surface'], b['mass_stranded'], b['mass_evaporated'], b['mass_biodegraded']))
        oil_density = b['oil_density']

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
        if np.sum(b['mass_biodegraded']) > 0:
            ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                          color='indigo', label='biodegraded'))
            ax1.fill_between(time, budget[4, :], budget[5, :],
                             facecolor='indigo')

        ax1.set_ylim([0, budget.max()])
        ax1.set_xlim([0, time.max()])
        ax1.set_ylabel('Mass oil  [%s]' %
                       self.elements.variables['mass_oil']['units'])
        ax1.set_xlabel('Time  [hours]')
        # Right axis showing volume
        ax2 = ax1.twinx()
        mass_total = b['mass_total'][-1]
        ax2.set_ylim([0, mass_total/oil_density])
        ax2.set_ylabel('Volume oil [m3]')
        if not hasattr(self, 'oil_name'):  # TODO
            self.oil_name = 'unknown oiltype'
            # TODO line below is dangerous when importing old files
            self.oil_name = self.get_config('seed:oil_type')
        plt.title('%s (%.1f kg/m3) - %s to %s' %
                  (self.oil_name,
                   oil_density,
                   self.start_time.strftime('%Y-%m-%d %H:%M'),
                   self.time.strftime('%Y-%m-%d %H:%M')))
        # Shrink current axis's height by 10% on the bottom
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax1.legend(bbox_to_anchor=(0., -0.10, 1., -0.03), loc=1,
                   ncol=6, mode="expand", borderaxespad=0., fontsize=10)
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
                    self.logger.warning('Two oils found with name %s (ADIOS IDs %s and %s). Using the first.' % (oiltype, ADIOS_ids[0], ADIOS_ids[1]))
                    self.oiltype = OilProps(oils[0])
            except Exception as e:
                self.logger.warning(e)
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
        oilfile.close()
        self.oil_data['tref'] = np.array(tref, dtype='float')*3600.
        self.oil_data['fref'] = np.array(fref, dtype='float')*.01
        self.oil_data['wmax'] = np.array(wmax, dtype='float')
        self.oil_data['oiltype'] = oiltype  # Store name of oil type

    def seed_elements(self, *args, **kwargs):
        if 'oiltype' in kwargs:
            self.set_config('seed:oil_type', kwargs['oiltype'])
            del kwargs['oiltype']
        else:
            self.logger.info('Oil type not specified, using default: ' +
                         self.get_config('seed:oil_type'))
        self.set_oiltype(self.get_config('seed:oil_type'))

        if self.oil_weathering_model == 'noaa':
            try:  # Older version of OilLibrary
                oil_density = self.oiltype.get_density(285)  # 12 degrees
                oil_viscosity = self.oiltype.get_viscosity(285)
            except:  # Newer version of OilLibrary
                oil_density = self.oiltype.density_at_temp(285)
                oil_viscosity = self.oiltype.kvis_at_temp(285)
            self.logger.info('Using density %s and viscosity %s of oiltype %s' %
                         (oil_density, oil_viscosity, self.get_config('seed:oil_type')))
            kwargs['density'] = oil_density
            kwargs['viscosity'] = oil_viscosity

        if 'm3_per_hour' in kwargs:
            # From given volume rate, we calculate the mass per element
            num_elements = kwargs['number']
            time = kwargs['time']
            if type(time) is list:
                duration_hours = ((time[1] - time[0]).total_seconds())/3600
                if duration_hours == 0:
                    duration_hours = 1.
            else:
                duration_hours = 1.  # For instantaneous spill, we use 1h
            kwargs['mass_oil'] = (kwargs['m3_per_hour']*duration_hours/
                                  num_elements*kwargs['density'])
            del kwargs['m3_per_hour']

        super(OpenOil, self).seed_elements(*args, **kwargs)

    def seed_from_gml(self, gmlfile, num_elements=1000, *args, **kwargs):
        """Read oil slick contours from GML file, and seed particles within."""

        # Specific imports
        import datetime
        from matplotlib.path import Path
        from xml.etree import ElementTree
        from matplotlib.patches import Polygon
        import pyproj

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
        import osr

        if not 'time' is kwargs:
            try:  # get time from filename
                timestr = filename[-28:-13]
                time = datetime.strptime(
                        filename[-28:-13], '%Y%m%d.%H%M%S')
                self.logger.info('Parsed time from filename: %s' % time)
            except:
                time = datetime.now()
                self.logger.warning('Could not pase time from filename, '
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

        src_srs = osr.SpatialReference()
        src_srs.ImportFromEPSG(4269)
        tgt_srs = osr.SpatialReference()
        tgt_srs.ImportFromEPSG(3857)
        transform = osr.CoordinateTransformation(src_srs, tgt_srs)
        for cat in categories:
            memshapename = filename + '%i.shp' % cat
            layers[cat-1] = mem_vector_layers[cat-1]
            areas = np.zeros(layers[cat-1].GetFeatureCount())
            for i, feature in enumerate(layers[cat-1]):
                geom = feature.GetGeometryRef()
                geom.Transform(transform)  # To get area in m2
                areas[i] = geom.GetArea()
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
        oil_density = 1000
        mass_oil = (total_area*thickness_microns/1e6)*oil_density

        for i, num in enumerate(numbers):
            self.seed_from_shapefile([mem_vector_layers[i]],
                oil_film_thickness=thickness_microns[i]/1000000.,
                mass_oil=mass_oil[i]/num,
                number=num, time=time, *args, **kwargs)

