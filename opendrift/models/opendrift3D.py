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

import numpy as np
import logging
from scipy.interpolate import interp1d
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray


# Defining the oil element properties
class Lagrangian3DArray(LagrangianArray):
    """Extending LagrangianArray for elements moving in 3 dimensions
    The Particle may be buoyant and/or subject to vertical mixing
    buoyant bahaviour is described by terminal velocity
    """

    variables = LagrangianArray.add_variables([
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': 0.})])


class OpenDrift3DSimulation(OpenDriftSimulation):
    """Open source buoyant particle trajectory model based on OpenDrift.

        Developed at MET Norway

        Generic module for particles that move in 3 dimensions
        and may be to vertical turbulent mixing
        with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, nutrients or sediments

        Under construction.
    """

    max_speed = 1  # m/s

    def __init__(self, *args, **kwargs):

        configspec_oceandrift3D = '''
            [processes]
                turbulentmixing = boolean(default=True)
                verticaladvection = boolean(default=True)
            [turbulentmixing]
                timestep = float(min=0.1, max=3600, default=60.)
                verticalresolution = float(min=0.01, max=10, default = 1.)
                diffusivitymodel = option('environment', 'stepfunction', 'windspeed_Sundby1983', 'windspeed_Large1994', 'gls_tke', default='environment')
                TSprofiles = boolean(default=False)
                '''
        self._add_configstring(configspec_oceandrift3D)

        # Calling general constructor of parent class
        super(OpenDrift3DSimulation, self).__init__(*args, **kwargs) 

    def update_terminal_velocity(self, Tprofiles=None, Sprofiles=None,
                                 z_index=None):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables
        overload this function to create particle-specific behaviour
        """
        #self.elements.terminal_velocity = 0.00

    def vertical_advection(self):
        """Move particles vertically according to vertical ocean current

            Vertical advection by ocean currents is small compared to
            termical velocity
        """
        if self.get_config('processes:verticaladvection') is False:
            logging.debug('Vertical advection deactivated')
            return

        in_ocean = np.where(self.elements.z<0)[0]
        if len(in_ocean) > 0:
            w = self.environment.upward_sea_water_velocity[in_ocean]
            self.elements.z[in_ocean] = np.minimum(0,
                self.elements.z[in_ocean] + w * self.time_step.total_seconds())
        else:
            logging.debug('No vertical advection for elements at surface')

    def prepare_vertical_mixing(self):
        pass  # To be implemented by subclasses as needed

    def surface_stick(self):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        # keep particle just below the surface
        surface = np.where(self.elements.z >= 0)
        if len(surface[0]) > 0:
            self.elements.z[surface] = -0.01

    def surface_wave_mixing(self, time_step_seconds):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        # do nothing 
        pass

    def vertical_mixing(self):
        """Mix particles vertically according to eddy diffusivity and buoyancy

            Buoyancy is expressed as terminal velocity, which is the
            steady-state vertical velocity due to positive or negative
            buoyant behaviour. It is usually a function of particle density,
            diameter, and shape.

            Vertical particle displacemend du to turbulent mixing is
            calculated using a random walk scheme" (Visser et al. 1996)
        """

        if self.get_config('processes:turbulentmixing') is False:
            logging.debug('Turbulent mixing deactivated')
            return

        self.timer_start('main loop:updating elements:vertical mixing')
        from opendrift.models import eddydiffusivity

        dt_mix = self.get_config('turbulentmixing:timestep')

        # minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # Eventual model specific preparions
        self.prepare_vertical_mixing()

        # get profile of eddy diffusivity
        # get vertical eddy diffusivity from environment or specific model
        if (self.get_config('turbulentmixing:diffusivitymodel') ==
                'environment'):
            if 'ocean_vertical_diffusivity' in self.environment_profiles:
                Kprofiles = self.environment_profiles[
                    'ocean_vertical_diffusivity']
                logging.debug('Using diffusivity from ocean model')
            else:
                # NB: using constant diffusivity, and value from first
                # element only - this should be checked/improved!
                Kprofiles = \
                    self.environment.ocean_vertical_diffusivity[0] * \
                    np.ones((len(self.environment_profiles['z']),
                             self.num_elements_active()))
                logging.debug('Using constant diffusivity')
        else:
            logging.debug('Using functional expression for diffusivity')
            Kprofiles = getattr(
                eddydiffusivity,
                self.get_config('turbulentmixing:diffusivitymodel'))(self)

        logging.debug('Diffiusivities are in range %s to %s' %
                      (Kprofiles.min(), Kprofiles.max()))

        # get profiles of salinity and temperature
        # (to save interpolation time in the inner loop)
        if (self.get_config('turbulentmixing:TSprofiles') is True 
            and 'sea_water_salinity' in self.required_variables):
            Sprofiles = self.environment_profiles['sea_water_salinity']
            Tprofiles = \
                self.environment_profiles['sea_water_temperature']
            if ('sea_water_salinity' in self.fallback_values and
                Sprofiles.min() == Sprofiles.max()):
                logging.debug('Salinity and temperature are fallback'
                              'values, skipping TSprofile')
                Sprofiles = None
                Tprofiles = None
            else:
                logging.debug('Using TSprofiles for vertical mixing')
        else:
            logging.debug('TSprofiles deactivated for vertical mixing')
            Sprofiles = None
            Tprofiles = None

        # prepare vertical interpolation coordinates
        #z_i = range(Kprofiles.shape[0])
        z_i = range(self.environment_profiles['z'].shape[0])
        #print len(self.environment_profiles['z']), len(z_i)
        z_index = interp1d(-self.environment_profiles['z'],
                           z_i, bounds_error=False,
                           fill_value=(0,len(z_i)-1))  # Extrapolation

        # internal loop for fast time step of vertical mixing model
        # random walk needs faster time step compared
        # to horizontal advection
        logging.debug('Vertical mixing module:')
        ntimes_mix = np.abs(int(self.time_step.total_seconds()/dt_mix))
        logging.debug('Turbulent diffusion with random walk '
                      'scheme using ' + str(ntimes_mix) +
                      ' fast time steps of dt=' + str(dt_mix) + 's')

        for i in range(0, ntimes_mix):
            #remember which particles belong to the exact surface
            surface = self.elements.z == 0

            # update terminal velocity according to environmental variables
            if self.get_config('turbulentmixing:TSprofiles') is True:
                self.update_terminal_velocity(Tprofiles=Tprofiles,
                                              Sprofiles=Sprofiles,
                                              z_index=z_index)
            else:
                # this is faster, but ignores density gradients in
                # water column for the inner loop
                self.update_terminal_velocity()

            w = self.elements.terminal_velocity

            # diffusivity K at depth z+dz


            dz = 1e-3
            zi = z_index(-self.elements.z+0.5*dz)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Kprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
            weight_upper[np.isnan(weight_upper)] = 1
            K1 = Kprofiles[upper, range(Kprofiles.shape[1])] * \
                weight_upper + \
                Kprofiles[lower, range(Kprofiles.shape[1])] * \
                (1-weight_upper)

            # diffusivity K at depth z-dz
            zi = z_index(-self.elements.z-0.5*dz)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Kprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
            weight_upper[np.isnan(weight_upper)] = 1
            K2 = Kprofiles[upper, range(Kprofiles.shape[1])] * \
                weight_upper + \
                Kprofiles[lower, range(Kprofiles.shape[1])] * \
                (1-weight_upper)

            # diffusivity gradient
            dKdz = (K1 - K2) / dz

            # K at depth z+dKdz*dt/2 
            zi = z_index(-(self.elements.z+dKdz*dt_mix/2))
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Kprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
            weight_upper[np.isnan(weight_upper)] = 1
            K3 = Kprofiles[upper, range(Kprofiles.shape[1])] * \
                weight_upper + \
                Kprofiles[lower, range(Kprofiles.shape[1])] * \
                (1-weight_upper)


            # Visser et al. 1996 random walk mixing
            # requires an inner loop time step dt such that
            # dt << (d2K/dz2)^-1, e.g. typically dt << 15min
            R = 2*np.random.random(self.num_elements_active()) - 1            
            r = 1.0/3
            # new position  =  old position   - up_K_flux   + random walk
            self.elements.z = self.elements.z - dKdz*dt_mix + R*np.sqrt(( K3*dt_mix*2/r))
 
            # Reflect from surface 
            reflect = np.where(self.elements.z >= 0)
            if len(reflect[0]) > 0:
                self.elements.z[reflect] = -self.elements.z[reflect]

            # reflect elements going below seafloor
            bottom = np.where(self.elements.z < Zmin)
            if len(bottom[0]) > 0:
                logging.debug('%s elements penetrated seafloor, lifting up' % len(bottom[0]))
                self.elements.z[bottom] = 2*Zmin[bottom] - self.elements.z[bottom]
           
            # advect due to buoyancy
            self.elements.z = self.elements.z + w*dt_mix

            # put the particles that belonged to the surface slick (if present) back to the surface
            self.elements.z[surface] = 0.

            # formation of slick and wave mixing for surfaced particles if implemented for this class
            self.surface_stick()
            self.surface_wave_mixing(dt_mix)

            # let particles stick to bottom 
            bottom = np.where(self.elements.z < Zmin)
            if len(bottom[0]) > 0:
                logging.debug('%s elements reached seafloor, set to bottom' % len(bottom[0]))
                self.elements.z[bottom] = Zmin[bottom]
 
        self.timer_end('main loop:updating elements:vertical mixing')


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
            dz = self.get_config('turbulentmixing:verticalresolution')
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

    def plotter_vertical_distribution_time(self, ax=None, mask=None, dz=1., maxrange=-100, bins=None, step=1):
        """Function to plot vertical distribution of particles
	
	use mask to plot any selection of particles
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


