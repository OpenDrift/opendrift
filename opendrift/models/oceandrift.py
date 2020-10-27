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

import sys
from datetime import timedelta
import numpy as np
from scipy.interpolate import interp1d
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray
from opendrift.models.physics_methods import verticaldiffusivity_Large1994, verticaldiffusivity_Sundby1983, gls_tke

# Defining the oil element properties
class Lagrangian3DArray(LagrangianArray):
    """Extending LagrangianArray for elements moving in 3 dimensions
    The Particle may be buoyant and/or subject to vertical mixing
    buoyant bahaviour is described by terminal velocity
    """

    variables = LagrangianArray.add_variables([
        ('wind_drift_factor', {'dtype': np.float32,
                               'units': '1',
                               'default': 0.02}),
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': 0.})])


class OceanDrift(OpenDriftSimulation):
    """Open source buoyant particle trajectory model based on OpenDrift.

        Developed at MET Norway

        Generic module for particles that move in 3 dimensions
        and may be to vertical turbulent mixing
        with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, nutrients or sediments,
        Model may be subclassed for more specific behaviour.

    """

    ElementType = Lagrangian3DArray

    max_speed = 1  # m/s

    required_variables = [
        'x_sea_water_velocity',
        'y_sea_water_velocity',
        'x_wind', 'y_wind',
        'upward_sea_water_velocity',
        'ocean_vertical_diffusivity',
        'sea_surface_wave_significant_height',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'surface_downward_x_stress',
        'surface_downward_y_stress',
        'turbulent_kinetic_energy',
        'turbulent_generic_length_scale',
        'sea_floor_depth_below_sea_level',
        'land_binary_mask'
        ]

    required_profiles = ['ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-20, 0]

    fallback_values = {
        'x_sea_water_velocity': 0,
        'y_sea_water_velocity': 0,
        'upward_sea_water_velocity': 0,
        'sea_surface_wave_significant_height': 0,
        'sea_surface_wave_stokes_drift_x_velocity': 0,
        'sea_surface_wave_stokes_drift_y_velocity': 0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
        'x_wind': 0,
        'y_wind': 0,
        'ocean_vertical_diffusivity': 0,
        'surface_downward_x_stress': 0,
        'surface_downward_y_stress': 0,
        'turbulent_kinetic_energy': 0,
        'turbulent_generic_length_scale': 0,
        'sea_floor_depth_below_sea_level': 10000
        }


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(OceanDrift, self).__init__(*args, **kwargs) 

        self._add_config({
            'drift:vertical_advection': {'type': 'bool', 'default': True, 'description': 
                'Advect elements with vertical component of ocean current.',
                'level': self.CONFIG_LEVEL_BASIC},
            'drift:vertical_mixing': {'type': 'bool', 'default': False, 'level': self.CONFIG_LEVEL_BASIC,
                'description': 'Activate vertical mixing scheme with inner loop'},
            'vertical_mixing:timestep': {'type': 'float', 'min': 0.1, 'max': 3600, 'default': 60,
                'level': self.CONFIG_LEVEL_ADVANCED, 'units': 'seconds', 'description':
                'Time step used for inner loop of vertical mixing.'},
            'vertical_mixing:diffusivitymodel': {'type': 'enum', 'default': 'environment',
                'enum': ['environment', 'stepfunction', 'windspeed_Sundby1983',
                 'windspeed_Large1994', 'gls_tke','constant'], 'level': self.CONFIG_LEVEL_ADVANCED,
                 'units': 'seconds', 'description': 'Time step used for inner loop of vertical mixing.'},
            'vertical_mixing:TSprofiles': {'type': 'bool', 'default': False, 'level':
                self.CONFIG_LEVEL_ADVANCED,
                'description': 'Update T and S profiles within inner loop of vertical mixing.'},
            'drift:wind_drift_depth': {'type': 'float', 'default': 0.1,
                'min': 0, 'max': 10, 'units': 'meters',
                'description': 'The direct wind drift (windage) is linearly decreasing from the surface value (wind_drift_factor) until 0 at this depth.',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:stokes_drift': {'type': 'bool', 'default': True,
                'description': 'Advection elements with Stokes drift (wave orbital motion).',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:use_tabularised_stokes_drift': {'type': 'bool', 'default': False,
                'description': 'If True, Stokes drift is estimated from wind based on look-up-tables for given fetch (drift:tabularised_stokes_drift_fetch).',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:tabularised_stokes_drift_fetch': {'type': 'enum', 'enum': ['5000', '25000', '50000'], 'default': '25000',
                'level': self.CONFIG_LEVEL_ADVANCED, 'description':
                'The fetch length when using tabularised Stokes drift.'},
            'drift:lift_to_seafloor': {'type': 'bool', 'default': True,
                'description': 'If True, elements hitting/penetrating seafloor, are lifted to seafloor height. The alternative (False) is to deactivate elements).',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:truncate_ocean_model_below_m': {'type': 'float', 'default': None,
                'min': 0, 'max': 10000, 'units': 'm',
                'description': 'Ocean model data are only read down to at most this depth, and extrapolated below. May be specified to read less data to improve performance.',
                'level': self.CONFIG_LEVEL_ADVANCED},

            })

    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

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
            self.vertical_buoyancy()

        # Vertical advection
        self.vertical_advection()

    def disable_vertical_motion(self):
        """Deactivate any vertical processes/advection"""
        conf = {
                'drift:vertical_advection': False,
                'drift:vertical_mixing': False}
        for co, va in conf.items():
            self.logger.info('Setting config: %s -> %s' % (co, va))
            self.set_config(co, va)

    def update_terminal_velocity(self, Tprofiles=None, Sprofiles=None,
                                 z_index=None):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables. Sub-modules should overload
        this method for particle-specific behaviour
        """
        pass 

    def prepare_vertical_mixing(self):
        pass  # To be implemented by subclasses as needed

    def vertical_advection(self):
        """Move particles vertically according to vertical ocean current

            Vertical advection by ocean currents is normally small
            compared to termical velocity
        """
        if self.get_config('drift:vertical_advection') is False:
            self.logger.debug('Vertical advection deactivated')
            return

        in_ocean = np.where(self.elements.z<0)[0]
        if len(in_ocean) > 0:
            w = self.environment.upward_sea_water_velocity[in_ocean]
            self.elements.z[in_ocean] = np.minimum(0,
                self.elements.z[in_ocean] + w * self.time_step.total_seconds())
        else:
            self.logger.debug('No vertical advection for elements at surface')

    def vertical_buoyancy(self):
        """Move particles vertically according to their buoyancy"""
        in_ocean = np.where(self.elements.z<0)[0]
        if len(in_ocean) > 0:
            self.elements.z[in_ocean] = np.minimum(0,
                self.elements.z[in_ocean] + self.elements.terminal_velocity[in_ocean] * self.time_step.total_seconds())
	
        # check for minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # Let particles stick to bottom 
        bottom = np.where(self.elements.z < Zmin)
        if len(bottom[0]) > 0:
            self.logger.debug('%s elements reached seafloor, set to bottom' % len(bottom[0]))
            self.elements.z[bottom] = Zmin[bottom]
            self.bottom_interaction(Zmin)
	
    def surface_stick(self):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        # keep particle just below the surface
        surface = np.where(self.elements.z >= 0)
        if len(surface[0]) > 0:
            self.elements.z[surface] = -0.01
            
    def bottom_interaction(self, Zmin=None):
        '''To be overloaded by subclasses, e.g. radionuclides in sediments'''
        pass

    def surface_wave_mixing(self, time_step_seconds):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        pass

    def get_diffusivity_profile(self, model):
        depths = self.environment_profiles['z']
        wind, depth = np.meshgrid(self.wind_speed(), depths)

        if model == 'windspeed_Large1994':
            return verticaldiffusivity_Large1994(wind, depth)
        elif model == 'windspeed_Sundby1983':
            return verticaldiffusivity_Sundby1983(wind)
        elif model == 'gls_tke':
            if not hasattr(self, 'gls_parameters'):
                self.logger.info('Searching readers for GLS parameters...')
                for reader_name, reader in self.readers.items():
                    if hasattr(reader, 'gls_parameters'):
                        self.gls_parameters = reader.gls_parameters
                        self.logger.info('Found gls-parameters in ' + reader_name)
                        break  # Success
                if not hasattr(self, 'gls_parameters'):
                    self.logger.info('Did not find gls-parameters in any readers.')
                    self.gls_parameters = None
            windstress = np.sqrt(self.environment.surface_downward_x_stress**2 +
                                 self.environment.surface_downward_y_stress**2)
            return gls_tke(windstress, depth, self.sea_water_density(),
                           self.environment.turbulent_kinetic_energy,
                           self.environment.turbulent_generic_length_scale,
                           gls_parameters)

        else:
            raise ValueError('Unknown diffusivity model: ' + model)

    def vertical_mixing(self, store_depths=False):
        """Mix particles vertically according to eddy diffusivity and buoyancy

            Buoyancy is expressed as terminal velocity, which is the
            steady-state vertical velocity due to positive or negative
            buoyant behaviour. It is usually a function of particle density,
            diameter, and shape.

            Vertical particle displacemend du to turbulent mixing is
            calculated using a random walk scheme" (Visser et al. 1996)
        """

        if self.get_config('drift:vertical_mixing') is False:
            self.logger.debug('Turbulent mixing deactivated')
            return

        self.timer_start('main loop:updating elements:vertical mixing')

        dt_mix = self.get_config('vertical_mixing:timestep')

        # minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # Eventual model specific preparions
        self.prepare_vertical_mixing()

        # get profile of eddy diffusivity
        # get vertical eddy diffusivity from environment or specific model
        diffusivity_model = self.get_config('vertical_mixing:diffusivitymodel')
        if diffusivity_model == 'environment':
            if 'ocean_vertical_diffusivity' in self.environment_profiles and not (self.environment_profiles['ocean_vertical_diffusivity'].min() == self.fallback_values['ocean_vertical_diffusivity'] and self.environment_profiles['ocean_vertical_diffusivity'].max() == self.fallback_values['ocean_vertical_diffusivity']):
                Kprofiles = self.environment_profiles[
                    'ocean_vertical_diffusivity']
                self.logger.debug('Using diffusivity from ocean model')
            else:
                self.logger.debug('Using diffusivity from Large1994 since model diffusivities not available')
                # Using higher vertical resolution when analytical
                self.environment_profiles['z'] = -np.arange(0, 50)
                Kprofiles = self.get_diffusivity_profile('windspeed_Large1994')
        elif diffusivity_model == 'constant':
            self.logger.debug('Using constant diffusivity specified by fallback_values[''ocean_vertical_diffusivity''] = %s m2.s-1' % (self.fallback_values['ocean_vertical_diffusivity']))
            Kprofiles = self.fallback_values['ocean_vertical_diffusivity']*np.ones(
                    self.environment_profiles['ocean_vertical_diffusivity'].shape) # keep constant value for ocean_vertical_diffusivity
        else:
            self.logger.debug('Using functional expression for diffusivity')
            # Using higher vertical resolution when analytical
            self.environment_profiles['z'] = -np.arange(0, 50)
            # Note: although analytical functions, z is discretised
            Kprofiles = self.get_diffusivity_profile(diffusivity_model)

        self.logger.debug('Diffusivities are in range %s to %s' %
                      (Kprofiles.min(), Kprofiles.max()))

        # get profiles of salinity and temperature
        # (to save interpolation time in the inner loop)
        if (self.get_config('vertical_mixing:TSprofiles') is True 
            and 'sea_water_salinity' in self.required_variables):
            Sprofiles = self.environment_profiles['sea_water_salinity']
            Tprofiles = \
                self.environment_profiles['sea_water_temperature']
            if ('sea_water_salinity' in self.fallback_values and
                Sprofiles.min() == Sprofiles.max()):
                self.logger.debug('Salinity and temperature are fallback'
                              'values, skipping TSprofile')
                Sprofiles = None
                Tprofiles = None
            else:
                self.logger.debug('Using TSprofiles for vertical mixing')
        else:
            self.logger.debug('TSprofiles deactivated for vertical mixing')
            Sprofiles = None
            Tprofiles = None

        # prepare vertical interpolation coordinates
        z_i = range(self.environment_profiles['z'].shape[0])
        if len(z_i) == 1:
            z_index = 0
        else:
            z_index = interp1d(-self.environment_profiles['z'],
                               z_i, bounds_error=False,
                               fill_value=(0,len(z_i)-1))  # Extrapolation

        # Internal loop for fast time step of vertical mixing model.
        # Random walk needs faster time step than horizontal advection.
        self.logger.debug('Vertical mixing module:' +
            self.get_config('vertical_mixing:diffusivitymodel'))
        ntimes_mix = np.abs(int(self.time_step.total_seconds()/dt_mix))
        self.logger.debug('Turbulent diffusion with random walk '
                      'scheme using ' + str(ntimes_mix) +
                      ' fast time steps of dt=' + str(dt_mix) + 's')

        if store_depths is not False:
            depths = np.zeros((ntimes_mix, self.num_elements_active()))
            depths[0, :] = self.elements.z

        # Calculating dK/dz for all profiles before the loop
        gradK = -np.gradient(Kprofiles, self.environment_profiles['z'], axis=0)
        gradK[np.abs(gradK)<1e-10] = 0

        for i in range(0, ntimes_mix):
            #remember which particles belong to the exact surface
            surface = self.elements.z == 0

            # Update the terminal velocity of particles
            self.update_terminal_velocity(Tprofiles=Tprofiles, Sprofiles=Sprofiles, z_index=z_index)
            w = self.elements.terminal_velocity

            # Diffusivity and its gradient at z
            zi = np.round(z_index(-self.elements.z)).astype(np.int)
            Kz = Kprofiles[zi, range(Kprofiles.shape[1])]
            dKdz = gradK[zi, range(Kprofiles.shape[1])]

            # Visser et al. 1996 random walk mixing
            # requires an inner loop time step dt such that
            # dt << (d2K/dz2)^-1, e.g. typically dt << 15min
            R = 2*np.random.random(self.num_elements_active()) - 1            
            r = 1.0/3
            # New position  =  old position   - up_K_flux   + random walk
            self.elements.z = self.elements.z - self.elements.moving*(
                dKdz*dt_mix - R*np.sqrt((Kz*dt_mix*2/r)))
 
            # Reflect from surface 
            reflect = np.where(self.elements.z >= 0)
            if len(reflect[0]) > 0:
                self.elements.z[reflect] = -self.elements.z[reflect]

            # Reflect elements going below seafloor
            bottom = np.where(self.elements.z < Zmin)
            if len(bottom[0]) > 0:
                self.logger.debug('%s elements penetrated seafloor, lifting up' % len(bottom[0]))
                self.elements.z[bottom] = 2*Zmin[bottom] - self.elements.z[bottom]
           
            # Advect due to buoyancy
            self.elements.z = self.elements.z + w*dt_mix*self.elements.moving

            # Put the particles that belonged to the surface slick
            # (if present) back to the surface
            self.elements.z[surface] = 0.

            # Formation of slick and wave mixing for surfaced particles
            # if implemented for this class
            self.surface_stick()
            self.surface_wave_mixing(dt_mix)

            # Let particles stick to bottom 
            bottom = np.where(self.elements.z < Zmin)
            if len(bottom[0]) > 0:
                self.logger.debug('%s elements reached seafloor, set to bottom' % len(bottom[0]))
                self.elements.z[bottom] = Zmin[bottom]
                self.bottom_interaction(Zmin)

            if store_depths is not False:
                depths[i, :] = self.elements.z
 
        self.timer_end('main loop:updating elements:vertical mixing')

        if store_depths is not False:
            return depths
        else:
            return None

    def animate_vertical_distribution(self, depths=None, maxdepth=None, bins=50, filename=None):
        """Function to animate vertical distribution of particles"""
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation

        fig, (axk, axn) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]})
        if depths is not None:  # Debug mode, output from one cycle has been provided
            z = depths
            time_step = self.get_config('vertical_mixing:timestep')
            times = [self.time + i*timedelta(seconds=time_step) for i in range(z.shape[0])]
        else:
            z = self.get_property('z')[0]
            z = np.ma.filled(z, np.nan)
            K = self.get_property('ocean_vertical_diffusivity')[0]
            time_step = self.time_step.total_seconds()
            times = self.get_time_array()[0]
        if maxdepth is None:
            maxdepth = np.nanmin(z)
        if maxdepth > 0:
            maxdepth = -maxdepth  # negative z

        if depths is not None:
            axk.plot(np.nanmean(self.environment_profiles['ocean_vertical_diffusivity'], 1), self.environment_profiles['z'])
            xmax = self.environment_profiles['ocean_vertical_diffusivity'].max()
        else:
            axk.plot(K, z, 'k.')
            xmax = np.nanmax(K)
        axk.set_ylim([maxdepth, 0])
        axk.set_xlim([0, xmax*1.1])
        axk.set_ylabel('Depth [m]')
        axk.set_xlabel('Vertical diffusivity [$m^2/s$]')

        hist_series = np.zeros((bins, len(times)))
        bin_series = np.zeros((bins+1, len(times)))
        for i in range(len(times)):
            hist_series[:,i], bin_series[:,i] = np.histogram(z[i,:][np.isfinite(z[i,:])], bins=bins)
        maxnum = hist_series.max()

        def update_histogram(i):
            axn.clear()
            axn.barh(bin_series[0:-1,i], hist_series[:,i], height=-maxdepth/bins, align='edge')
            axn.set_ylim([maxdepth, 0])
            axn.set_xlim([0, maxnum])
            axn.set_title('%s UTC' % times[i])
            axn.set_xlabel('Number of particles')
            #axn.set_ylabel('Depth [m]')

        animation = animation.FuncAnimation(fig, update_histogram, len(times))
        if filename is not None or 'sphinx_gallery' in sys.modules:
            self._save_animation(animation, filename, fps=10)
        else:
            plt.show()

    def plot_vertical_distribution(self):
        """Function to plot vertical distribution of particles"""
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider, Button, RadioButtons
        from pylab import axes, draw
        from matplotlib import dates

        fig = plt.figure()
        mainplot = fig.add_axes([.15, .3, .8, .5])
        sliderax = fig.add_axes([.15, .08, .75, .05])
        data = self.history['z'].T[1, :]
        tslider = Slider(sliderax, 'Timestep', 0, self.steps_output-1,
                         valinit=self.steps_output-1, valfmt='%0.0f')
        try:
            dz = self.get_config('vertical_mixing:verticalresolution')
        except:
            dz = 1.
        maxrange = -100

        def update(val):
            tindex = np.int(tslider.val)
            mainplot.cla()
            mainplot.grid()
            mainplot.hist(self.history['z'].T[tindex, :], bins=int(-maxrange/dz),
                          range=[maxrange, 0], orientation='horizontal')
            mainplot.set_ylim([maxrange, 0])
            mainplot.set_xlabel('number of particles')
            mainplot.set_ylabel('depth [m]')
            x_wind = self.history['x_wind'].T[tindex, :]
            y_wind = self.history['y_wind'].T[tindex, :]
            windspeed = np.mean(np.sqrt(x_wind**2 + y_wind**2))
            mainplot.set_title(str(self.get_time_array()[0][tindex]) +
                               #'   Percent at surface: %.1f %' % percent_at_surface)
                               '   Mean windspeed: %.1f m/s' % windspeed)
            draw()

        update(0)  # Plot initial distribution
        tslider.on_changed(update)
        plt.show()

    def plotter_vertical_distribution_time(self, ax=None, mask=None,
            dz=1., maxrange=-100, bins=None, step=1):
        """Function to plot vertical distribution of particles.
	
	Use mask to plot any selection of particles.
	"""
        from pylab import axes, draw
        from matplotlib import dates, pyplot

        if ax is None:
            fig = pyplot.figure()
            ax = fig.gca()
            show = True
        else:
            show = False

        if mask is None: # create a mask that is True for all particles
            mask = self.history['z'].T[0] == self.history['z'].T[0] 

        if bins is None:
            bins=int(-maxrange/dz)

        ax.hist(self.history['z'].T[step,mask], bins=bins,
                range=[maxrange, 0], orientation='horizontal')
        ax.set_ylim([maxrange, 0])
        ax.grid()
        #ax.set_xlim([0, mask.sum()*.15])
        ax.set_xlabel('Number of particles')
        ax.set_ylabel('Depth [m]')
        x_wind = self.history['x_wind'].T[step, :]
        y_wind = self.history['x_wind'].T[step, :]
        windspeed = np.mean(np.sqrt(x_wind**2 + y_wind**2))
        ax.set_title(str(self.get_time_array()[0][step]) +
                     '   Mean windspeed: %.1f m/s' % windspeed)
        if show is True:
            pyplot.show()
