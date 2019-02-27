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

import os
import numpy as np
from datetime import datetime
import logging
import json

from opendrift.models.openoil import OpenOil, Oil
from opendrift.models.opendrift3D import OpenDrift3DSimulation

try:
    basestring
except NameError:
    basestring = str


# Defining the oil element properties
class Oil3D(Oil):
    """Extending Oil class with variables relevant for the vertical."""

    variables = Oil.add_variables([
        # Entrainment length scale, see Tkalich and Chan (2002)
        ('entrainment_length_scale', {'dtype': np.float32,
                                      'units': 'm',
                                      'default': 0.03}),
        ('diameter', {'dtype': np.float32,  # Particle diameter
                      'units': 'm',
                      'default': 0.})
        ])


class OpenOil3D(OpenDrift3DSimulation, OpenOil):  # Multiple inheritance
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil3D

    required_variables = [
        'x_sea_water_velocity', 'y_sea_water_velocity',
        'sea_surface_wave_significant_height',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'sea_ice_area_fraction',
        'x_wind', 'y_wind', 'land_binary_mask',
        'sea_floor_depth_below_sea_level',
        'ocean_vertical_diffusivity',
        'sea_water_temperature',
        'sea_water_salinity',
        'upward_sea_water_velocity'
        ]

    # Desired variables do not require initialisation of Lazy readers
    desired_variables = [
        'sea_surface_wave_significant_height',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'sea_ice_area_fraction',
        'ocean_vertical_diffusivity',
        'upward_sea_water_velocity'
        ]

    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    fallback_values = {
        #'x_sea_water_velocity': 0,
        #'y_sea_water_velocity': 0,
        'sea_surface_wave_significant_height': 0,
        'sea_surface_wave_stokes_drift_x_velocity': 0,
        'sea_surface_wave_stokes_drift_y_velocity': 0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
        'sea_ice_area_fraction': 0,
        #'x_wind': 0, 'y_wind': 0,
        'sea_floor_depth_below_sea_level': 10000,
        'ocean_vertical_diffusivity': 0.02,  # m2s-1
        'sea_water_temperature': 10.,
        'sea_water_salinity': 34.,
        'upward_sea_water_velocity': 0
        }

    max_speed = 1.3  # m/s

    # Read oil types from file (presently only for illustrative effect)
    oil_types = str([str(l.strip()) for l in open(
                    os.path.dirname(os.path.realpath(__file__)) +
                    '/oil_types.txt').readlines()])[1:-1]
    default_oil = oil_types.split(',')[0].strip()

    # Configuration
    configspecOO3D = '''
        [input]
            [[spill]]
                oil_type = option(%s, default=%s)
                droplet_diameter_min_subsea = float(min=1e-8, max=1, default=0.0005)
                droplet_diameter_max_subsea = float(min=1e-8, max=1, default=0.005)
        [drift]
            wind_drift_depth = float(min=0, max=10, default=0.1)
            verticaladvection = boolean(default=False)
        [wave_entrainment]
            droplet_size_distribution = option('Exponential', 'Johansen et al. (2015)', 'Li et al. (2017)', default='Johansen et al. (2015)')
            entrainment_rate = option('Tkalich & Chan (2002)', 'Li et al. (2017)', default='Li et al. (2017)')
        [turbulentmixing]
            droplet_diameter_min_wavebreaking = float(default=1e-5, min=1e-8, max=1)
            droplet_diameter_max_wavebreaking = float(default=2e-3, min=1e-8, max=1)
            droplet_size_exponent = float(default=0, min=-10, max=10)
    ''' % (oil_types, default_oil)

    def __init__(self, *args, **kwargs):

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

        self._add_configstring(self.configspecOO3D)

        # Calling general constructor of parent class
        super(OpenOil3D, self).__init__(*args, **kwargs)

    def seed_elements(self, *args, **kwargs):

        if len(args) == 2:
            kwargs['lon'] = args[0]
            kwargs['lat'] = args[1]
            args = {}

        seed_json = {'start':{}, 'end':{}}
        for kw in kwargs:
            data = kwargs[kw]
            if not isinstance(data, basestring):
                data = np.atleast_1d(data)
            if isinstance(data[0], datetime):
                data[0] = str(data[0])
                if len(data) == 2:
                    data[1] = str(data[1])
            if not isinstance(kwargs[kw], basestring):
                if kw in ['lon', 'lat', 'z', 'radius', 'time']:
                    pointer = seed_json['start']
                    pointer2 = seed_json['end']
                else:
                    pointer = seed_json
                if len(data) == 1:
                    self.add_metadata('seed_' + kw, data[0])
                    pointer[kw] = data[0]
                elif len(kwargs[kw]) == 2:
                    self.add_metadata('seed_' + kw + '_start', str(kwargs[kw][0]))
                    self.add_metadata('seed_' + kw + '_end', str(kwargs[kw][1]))
                    pointer[kw] = data[0]
                    pointer2[kw] = data[1]
                else:
                    pass
                    #logging.info('Not adding array %s to metadata' % kw)
            else:
                self.add_metadata('seed_' + kw, str(data))
                seed_json[kw] = data

        if 'number' not in kwargs:
            number = 1
        else:
            number = kwargs['number']
        if 'diameter' in kwargs:
            logging.info('Droplet diameter is provided, and will '
                         'be kept constant during simulation')
            self.keep_droplet_diameter = True
        else:
            self.keep_droplet_diameter = False
        if 'z' not in kwargs:
            kwargs['z'] = 0
        if isinstance(kwargs['z'], basestring) and \
                kwargs['z'][0:8] == 'seafloor':
            z = -np.ones(number)
        else:
            z = np.atleast_1d(kwargs['z'])
        if len(z) == 1:
            z = z*np.ones(number)  # Convert scalar z to array
        subsea = z < 0
        if np.sum(subsea) > 0 and 'diameter' not in kwargs:
            # Droplet min and max for particles seeded below sea surface
            sub_dmin = self.get_config('input:spill:droplet_diameter_min_subsea')
            sub_dmax = self.get_config('input:spill:droplet_diameter_max_subsea')
            logging.info('Using particle diameters between %s and %s m for '
                         'elements seeded below sea surface.' %
                         (sub_dmin, sub_dmax))
            kwargs['diameter'] = np.random.uniform(sub_dmin, sub_dmax, number)

        super(OpenOil3D, self).seed_elements(*args, **kwargs)

        # Add oil metadata
        try:
            self.add_metadata('seed_oil_density', self.oiltype.get_density())
            seed_json['oil_density'] = self.oiltype.get_density()
        except:
            try:
                self.add_metadata('seed_oil_density',
                                  self.oiltype.density_at_temp(283))
                seed_json['oil_density'] = self.oiltype.density_at_temp(283)
            except:
                pass
        try:
            self.add_metadata('seed_oil_viscosity',
                              self.oiltype.get_viscosity(283))
            seed_json['oil_viscosity'] = self.oiltype.get_viscosity(283)
        except:
            try:
                self.add_metadata('seed_oil_viscosity',
                                  self.oiltype.kvis_at_temp(283))
                seed_json['oil_viscosity'] = self.oiltype.kvis_at_temp(283)
            except:
                pass

        if not hasattr(self, 'seed_json'):
            self.seed_json = []
        self.seed_json.append(seed_json)

        class MyEncoder(json.JSONEncoder):  # Serialisation of json
            def default(self, obj):
                if isinstance(obj, np.bool_) :
                    return str(obj)
                elif isinstance(obj, bool) :
                    return str(obj)
                elif isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                else:
                    return super(MyEncoder, self).default(obj) 

        self.add_metadata('seed_json', json.dumps(self.seed_json,
                                            cls=MyEncoder))

    def prepare_run(self):
        super(OpenOil3D, self).prepare_run()

    def particle_radius(self):
        """Calculate radius of entained particles.

        Per now a fixed radius, should later use a distribution.
        """

        # Delvigne and Sweeney (1988)
        # rmax = 1818*np.power(self.wave_energy_dissipation(), -0.5)* \
        #             np.power(self.elements.viscosity, 0.34) / 1000000

        # r = np.random.uniform(0, rmax, self.num_elements_active())
        return self.elements.diameter/2.0  # Hardcoded diameter

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for oil droplets

        according to
        Tkalich et al. (2002): Vertical mixing of oil droplets
                               by breaking waves
        Marine Pollution Bulletin 44, 1219-1229

        If profiles of temperature and salt are passed into this function,
        they will be interpolated from the profiles.
        if not, T,S will be fetched from reader.
        """
        g = 9.81  # ms-2

        r = self.particle_radius()*2.0

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

        rho_oil = self.elements.density
        rho_water = self.sea_water_density(T=T0, S=S0)

        # dynamic water viscosity
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0)
        # ~0.0014 kg m-1 s-1
        # kinemativ water viscosity
        ny_w = my_w / rho_water
        rhopr = rho_oil/rho_water

        # terminal velocity for low Reynolds numbers
        kw = 2*g*(1-rhopr)/(9*ny_w)
        W = kw * r**2

        # check if we are in a high Reynolds number regime
        Re = 2*r*W/ny_w
        highRe = np.where(Re > 50)

        # Terminal velocity in high Reynolds numbers
        kw = (16*g*(1-rhopr)/3)**0.5
        W2 = kw*r**0.5

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    def oil_wave_entrainment_rate(self):
        er = self.get_config('wave_entrainment:entrainment_rate')
        if er == 'Tkalich & Chan (2002)':
            entrainment_rate = self.oil_wave_entrainment_rate_tkalich2002()
        elif er == 'Li et al. (2017)':
            entrainment_rate = self.oil_wave_entrainment_rate_li2017()
        return entrainment_rate

    def oil_wave_entrainment_rate_li2017(self):
        # Z. Li, M.L. Spaulding, D. French McCay, J. Mar. Pollut. Bull. (2016):
        # An algorithm for modeling entrainment and naturally and chemically dispersed
        # oil droplet size distribution under surface breaking wave conditions
        g = 9.81
        interfacial_tension = self.oil_water_interfacial_tension
        delta_rho = self.sea_water_density() - self.elements.density
        d_o = 4 * (interfacial_tension / (delta_rho*g))**0.5
        we = ( self.sea_water_density() * g * self.significant_wave_height() * d_o ) / interfacial_tension
        oh = self.elements.viscosity * self.elements.density * (self.elements.density * interfacial_tension * d_o )**-0.5 # From kin. to dyn. viscosity by * density
        entrainment_rate = 4.604e-10 * we**1.805 *oh**-1.023 * self.sea_surface_wave_breaking_fraction()
        return entrainment_rate

    def oil_wave_entrainment_rate_tkalich2002(self):
        # Tkalich P. and Chan E.S.  
        # Vertical mixing of oil droplets by breaking waves
        # Marine Pollution Bulletin. 2002, V.44 (11), pp. 1219-1229
        kb = 0.4
        omega = (2.*np.pi)/self.wave_period()
        gamma = self.wave_damping_coefficient()
        alpha = 1.5
        Low = self.elements.entrainment_length_scale
        entrainment_rate = \
            kb*omega*gamma*self.significant_wave_height() / \
            (16*alpha*Low)
        return entrainment_rate

    def prepare_vertical_mixing(self):
        '''Calculate entrainment probability before main loop'''
        self.oil_entrainment_probability = \
            self.oil_wave_entrainment_rate()*\
                self.get_config('turbulentmixing:timestep')
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
        entrained = np.logical_and(surface,
                        random_number<self.oil_entrainment_probability)

        # Intrusion depth for wave entrainment from Delvigne and Sweeney (1988), Li et al. (2017):
        if entrained.sum() > 0:
            logging.debug('Entraining %i of %i surface elements' %
						  (entrained.sum(), surface.sum()))
            zb = 1.5 * self.significant_wave_height() # between 0 and zb
            intrusion_depth = np.random.uniform(0, np.mean(zb), entrained.sum())
            self.elements.z[entrained] = - intrusion_depth
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
            d = self.get_wave_breaking_droplet_diameter_johansen2015()
        elif dm == 'Li et al. (2017)':
            d = self.get_wave_breaking_droplet_diameter_liz2017()
        elif dm == 'Exponential':
            d = self.get_wave_breaking_droplet_diameter_exponential()
        return d

    def get_wave_breaking_droplet_diameter_exponential(self):
        if not hasattr(self, 'droplet_spectrum_pdf'):
            # Generate droplet spectrum, if not already done
            logging.debug('Generating wave breaking droplet size spectrum')
            s = self.get_config('turbulentmixing:droplet_size_exponent')
            dmax = self.get_config('turbulentmixing:droplet_diameter_max_wavebreaking')
            dmin = self.get_config('turbulentmixing:droplet_diameter_min_wavebreaking')
            # Note: a long array of diameters is required for 
            # sufficient resolution at both ends of logarithmic scale.
            # Could perhaps use logspace instead of linspace(?)
            self.droplet_spectrum_diameter = np.linspace(dmin, dmax, 1000000)
            spectrum = self.droplet_spectrum_diameter**s
            self.droplet_spectrum_pdf = spectrum/np.sum(spectrum)

        return np.random.choice(self.droplet_spectrum_diameter,
                                size=self.num_elements_active(),
                                p=self.droplet_spectrum_pdf)

    def get_wave_breaking_droplet_diameter_liz2017(self):
        # Li,Zhengkai, M. Spaulding, D. French-McCay, D. Crowley, J.R. Payne: "Development of a unified oil droplet size distribution model 
        # with application to surface breaking waves and subsea blowout releases considering dispersant effects" Mar. Pol. Bul.
        # DOI: 10.1016/j.marpolbul.2016.09.008
        # Should be prefered when the oil film thickness is unknown.
        if not hasattr(self, 'droplet_spectrum_pdf'):
            # Generate droplet spectrum as in Li (Zhengkai) et al. (2017)
            logging.debug('Generating wave breaking droplet size spectrum')
            dmax = self.get_config('turbulentmixing:droplet_diameter_max_wavebreaking')
            dmin = self.get_config('turbulentmixing:droplet_diameter_min_wavebreaking')
            self.droplet_spectrum_diameter = np.linspace(dmin, dmax, 1000000)
            g = 9.81
            interfacial_tension = self.oil_water_interfacial_tension
            delta_rho = self.sea_water_density() - self.elements.density
            d_o = 4 * (interfacial_tension / (delta_rho*g))**0.5
            we = ( self.sea_water_density() * g * self.significant_wave_height() * d_o ) / interfacial_tension
            oh = self.elements.viscosity * self.elements.density * (self.elements.density * interfacial_tension * d_o )**-0.5 # From kin. to dyn. viscosity by * density
            r = 1.791
            p = 0.460
            q = -0.518
            dV_50 = d_o * r * (1+10*oh)**p * we**q # median droplet diameter in volume distribution
            sd = 0.4 # log standard deviation in log10 units
            Sd = np.log(10) *sd # log standard deviation in natural log units
            # TODO: calculation below with scalars, but we have arrays, with varying oil properties
            # treat all particle in one go:
            dV_50 = np.mean(dV_50) # mean log diameter
            dN_50 = np.exp( np.log(dV_50) - 3*Sd**2 ) # convert number distribution to volume distribution
            logging.debug('Droplet distribution median diameter dV_50: %f, dN_50: %f ' %( dV_50, np.mean(dN_50)))
            spectrum = (np.exp(-(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 / (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd * np.sqrt(2 * np.pi))
            self.droplet_spectrum_pdf = spectrum/np.sum(spectrum)
        if ~np.isfinite(np.sum(self.droplet_spectrum_pdf)) or \
                np.abs(np.sum(self.droplet_spectrum_pdf) - 1) > 1e-6:
            logging.warning('Could not update droplet diameters.')
            return self.elements.diameter
        else:
            return np.random.choice(self.droplet_spectrum_diameter,
                                    size=self.num_elements_active(),
                                    p=self.droplet_spectrum_pdf)


    def get_wave_breaking_droplet_diameter_johansen2015(self):
        # Johansen O, Reed M, Bodsberg NR, Natural dispersion revisited
        # DOI: 10.1016/j.marpolbul.2015.02.026
        # requires oil film thickness
        if not hasattr(self, 'droplet_spectrum_pdf') or self.get_config('processes:update_oilfilm_thickness') is True:
            # Generate droplet spectrum as in Johansen et al. (2015)
            logging.debug('Generating wave breaking droplet size spectrum')
            dmax = self.get_config('turbulentmixing:droplet_diameter_max_wavebreaking')
            dmin = self.get_config('turbulentmixing:droplet_diameter_min_wavebreaking')
            self.droplet_spectrum_diameter = np.linspace(dmin, dmax, 1000000)
            g = 9.81
            interfacial_tension = self.oil_water_interfacial_tension
            #
            #A = self.significant_wave_height()/2. # wave amplitude
            #re = (self.elements.density*self.elements.oil_film_thickness*(2*g*A)**0.5) / (self.elements.viscosity*self.elements.density) # Reyolds number
            #we = (self.elements.density*self.elements.oil_film_thickness*2*g*A) / interfacial_tension # Weber number
            #
            H = self.significant_wave_height() # fall height = 2 * wave amplitude
            # Reyolds number (Eq. 7a from Johansen et al. 2015)
            re = (self.elements.density*self.elements.oil_film_thickness*(g*H)**0.5) / (self.elements.viscosity*self.elements.density) 
            # Weber number (Eq. 7b from Johansen et al.2015)
            we = (self.elements.density*self.elements.oil_film_thickness*g*H) / interfacial_tension # Weber number
            A = 2.251 # parameters from Johansen et al. 2015
            Bp = 0.027
            B = A*Bp  
            dN_50 = (A*self.elements.oil_film_thickness*we**-0.6) + (B*self.elements.oil_film_thickness* re**-0.6) # median droplet diameter in number distribution
            sd = 0.4 # log standard deviation in log10 units
            Sd = np.log(10) *sd # log standard deviation in natural log units
            dV_50 = np.exp( np.log(dN_50) + 3*Sd**2 ) # convert number distribution to volume distribution
            # TODO: calculation below with scalars, but we have
            # arrays, with varying oil properties
            # treat all particle in one go:
            dV_50 = np.mean(dV_50) # mean log diameter
            logging.debug('Droplet distribution median diameter dV_50: %f, dN_50: %f ' %( dV_50, np.mean(dN_50)))
            spectrum = (np.exp(-(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 / (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd * np.sqrt(2 * np.pi))
            self.droplet_spectrum_pdf = spectrum/np.sum(spectrum)
        if ~np.isfinite(np.sum(self.droplet_spectrum_pdf)) or \
                np.abs(np.sum(self.droplet_spectrum_pdf) - 1) > 1e-6:
            logging.warning('Could not update droplet diameters.')
            return self.elements.diameter
        else:
            return np.random.choice(self.droplet_spectrum_diameter,
                                    size=self.num_elements_active(),
                                    p=self.droplet_spectrum_pdf)

    def resurface_elements(self, minimum_depth=None):
        """Oil elements reaching surface (or above) form slick, not droplet"""
        surface = np.where(self.elements.z >= 0)[0]
        self.elements.z[surface] = 0

    def update(self):
        """Update positions and properties of oil particles."""

        if self.get_config('processes:update_oilfilm_thickness') is True:
            self.update_surface_oilfilm_thickness()

        # Oil weathering (inherited from OpenOil)
        self.oil_weathering()

        # Turbulent Mixing
        if self.get_config('processes:turbulentmixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
            del self.droplet_spectrum_pdf

        # Vertical advection
        if self.get_config('processes:verticaladvection') is True:
            self.vertical_advection()

        # Horizontal advection (inherited from OpenOil)
        self.advect_oil()

