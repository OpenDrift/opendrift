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
            [drift]
                max_age_seconds = float(min=0, default=None)
            [processes]
                turbulentmixing = boolean(default=True)
                verticaladvection = boolean(default=True)
            [turbulentmixing]
                timestep = float(min=0.1, max=3600, default=4.)
                verticalresolution = float(min=0.01, max=10, default = 2.)
                diffusivitymodel = option('environment', 'stepfunction', 'windspeed_Sundby1983', 'gls_tke', default='environment')
                TSprofiles = boolean(default=True)
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
            logging.debug('Vertical advection deactivated.')
            return

        w = self.environment.upward_sea_water_velocity
        self.elements.z = np.minimum(0,
            self.elements.z + w * self.time_step.total_seconds())

    def prepare_vertical_mixing(self):
        pass  # To be implemented by subclasses as needed

    def surface_interaction(self, time_step_seconds):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        # Place particles above surface into the uppermost layer
        surface = self.elements.z >= 0
        self.elements.z[surface] = \
            -self.get_config('turbulentmixing:verticalresolution')/2.

    def vertical_mixing(self):
        """Mix particles vertically according to eddy diffusivity and buoyancy

            Buoyancy is expressed as terminal velocity, which is the
            steady-state vertical velocity due to positive or negative
            buoyant behaviour. It is usually a function of particle density,
            diameter, and shape.

            Vertical particle displacemend du to turbulent mixing is
            calculated using the "binned random walk scheme" (Thygessen and
            Aadlandsvik, 2007).
            The formulation of this scheme is copied from LADIM (IMR).
        """

        if self.get_config('processes:turbulentmixing') is False:
            logging.debug('Turbulent mixing deactivated.')
            return

        self.timer_start('main loop:updating elements:vertical mixing')
        from opendrift.models import eddydiffusivity

        dz = self.get_config('turbulentmixing:verticalresolution')
        dz = np.float32(dz)  # Convert to avoid error for older numpy
        dt_mix = self.get_config('turbulentmixing:timestep')

        # minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # place particle in center of bin
        surface = self.elements.z == 0
        self.elements.z[~surface] = np.round(self.elements.z[~surface]/dz)*dz

        #avoid that elements are below bottom
        bottom = np.where(self.elements.z < Zmin)
        if len(bottom[0]) > 0:
            self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

        # Eventual model specific preparions
        self.prepare_vertical_mixing()

        # get profile of eddy diffusivity
        # get vertical eddy diffusivity from environment or specific model
        if (self.get_config('turbulentmixing:diffusivitymodel') ==
                'environment'):
            if 'ocean_vertical_diffusivity' in self.environment_profiles:
                Kprofiles = self.environment_profiles[
                    'ocean_vertical_diffusivity']
                logging.debug('use diffusivity from ocean model')
            else:
                # NB: using constant diffusivity, and value from first
                # element only - this should be checked/improved!
                Kprofiles = \
                    self.environment.ocean_vertical_diffusivity[0] * \
                    np.ones((len(self.environment_profiles['z']),
                             self.num_elements_active()))
                logging.debug('use constant diffusivity')
        else:
            logging.debug('use functional expression for diffusivity')
            Kprofiles = getattr(
                eddydiffusivity,
                self.get_config('turbulentmixing:diffusivitymodel'))(self)

        logging.debug('Diffiusivities are in range %s to %s.' %
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
                logging.debug('Salinity and temperature are fallback values, '
                              'skipping TSprofile')
                Sprofiles = None
                Tprofiles = None
        else:
            Sprofiles = None
            Tprofiles = None

        # prepare vertical interpolation coordinates
        z_i = range(Kprofiles.shape[0])
        z_index = interp1d(-self.environment_profiles['z'],
                           z_i, bounds_error=False,
                           fill_value=(0,len(z_i)-1))  # Extrapolation
        # internal loop for fast time step of vertical mixing model
        # binned random walk needs faster time step compared
        # to horizontal advection
        ntimes_mix = np.abs(int(self.time_step.total_seconds()/dt_mix))
        logging.debug('Vertical mixing module:')
        logging.debug('turbulent diffusion with binned random walk scheme')
        logging.debug('using ' + str(ntimes_mix) + ' fast time steps of dt=' +
                      str(dt_mix) + 's')
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

            # diffusivity K at depth z
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Kprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
            weight_upper[np.isnan(weight_upper)] = 1
            K1 = Kprofiles[upper, range(Kprofiles.shape[1])] * \
                weight_upper + \
                Kprofiles[lower, range(Kprofiles.shape[1])] * \
                (1-weight_upper)

            # K at depth z-dz ; gradient of K is required for correct
            # solution with random walk scheme
            zi = z_index(-(self.elements.z-dz))
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Kprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
            weight_upper[np.isnan(weight_upper)] = 1
            K2 = Kprofiles[upper, range(Kprofiles.shape[1])] * \
                weight_upper + \
                Kprofiles[lower, range(Kprofiles.shape[1])] * \
                (1-weight_upper)

            # calculate rise/sink probability dependent on K and w
            p = dt_mix * (2.0*K1 + dz*w)/(2.0*dz*dz)  # probability to rise
            q = dt_mix * (2.0*K2 - dz*w)/(2.0*dz*dz)  # probability to sink

            # check if probabilities are reasonable or wrong; which can happen if K is very high (K>0.1)
            wrong = p+q > 1.00002
            if wrong.sum() > 0:
                logging.info('WARNING! '+str(wrong.sum())+' elements have p+q>1; you might need a smaller mixing time step')
                # fixing p and q by scaling them to assure p+q<1:
                norm = p+q
                p[wrong] = p[wrong]/norm[wrong] 
                q[wrong] = q[wrong]/norm[wrong]

            # use probabilities to mix some particles up or down
            RandKick = np.random.random(self.num_elements_active())           
            up = np.where(RandKick < p)
            down = np.where(RandKick > 1.0 - q)           
            self.elements.z[up] = self.elements.z[up] + dz # move to layer above
            self.elements.z[down] = self.elements.z[down] - dz # move to layer underneath

            # put the particles that belong to the surface slick (if present) back to the surface
            self.elements.z[surface] = 0.

            #avoid that elements are below bottom
            bottom = np.where(self.elements.z < Zmin)
            if len(bottom[0]) > 0:
                self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

            # Call surface interaction:
            # reflection at surface or formation of slick and wave mixing if implemented for this class
            self.surface_interaction(dt_mix)

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
            bins=-maxrange/dz

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
