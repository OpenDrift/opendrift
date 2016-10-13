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
        w = self.environment.upward_sea_water_velocity
        self.elements.z = np.minimum(0,
            self.elements.z + w * self.time_step.total_seconds())

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

        if self.config['processes']['turbulentmixing'] is False:
            logging.debug('Turbulent mixing deactivated.')
            return

        # if terminal_velocity is None:
        #     w = self.config['drift']['terminal_velocity']
        #else:
        #    w = terminal_velocity
        from opendrift.models import eddydiffusivity

        dz = self.config['turbulentmixing']['verticalresolution']
        dz = np.float32(dz)  # Convert to avoid error for older numpy
        dt_mix = self.config['turbulentmixing']['timestep']

        # minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # place particle in center of bin
        self.elements.z = np.round(self.elements.z/dz)*dz

        #avoid that elements are below bottom
        bottom = np.where(self.elements.z < Zmin)
        self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

        # get profile of eddy diffusivity
        # get vertical eddy diffusivity from environment or specific model
        if (self.config['turbulentmixing']['diffusivitymodel'] ==
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
                self.config['turbulentmixing']['diffusivitymodel'])(self)

        # get profiles of salinity and temperature
        # (to save interpolation time in the inner loop)
        if 'TSprofiles' in self.config['turbulentmixing']:
            if self.config['turbulentmixing']['TSprofiles'] is True:
                Sprofiles = self.environment_profiles['sea_water_salinity']
                Tprofiles = \
                    self.environment_profiles['sea_water_temperature']

        # prepare vertical interpolation coordinates
        z_i = range(Kprofiles.shape[0])
        z_index = interp1d(-self.environment_profiles['z'],
                           z_i, bounds_error=False)

        # internal loop for fast time step of vertical mixing model
        # binned random walk needs faster time step compared
        # to horizontal advection
        ntimes_mix = int(self.time_step.total_seconds()/dt_mix)
        logging.debug('Vertical mixing module:')
        logging.debug('turbulent diffusion with binned random walk scheme')
        logging.debug('using ' + str(ntimes_mix) + ' fast time steps of dt=' +
                      str(dt_mix) + 's')
        for i in range(0, ntimes_mix):

            # update terminal velocity according to environmental variables
            if 'TSprofiles' in self.config['turbulentmixing']:
                if self.config['turbulentmixing']['TSprofiles'] is True:
                    self.update_terminal_velocity(Tprofiles=Tprofiles,
                                                  Sprofiles=Sprofiles,
                                                  z_index=z_index)
                else:
                    self.update_terminal_velocity()
            else:
                # this is faster, but ignores density gradients in
                # water column for the inner loop
                self.update_terminal_velocity()

            w = self.elements.terminal_velocity

            # K at depth z
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Kprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
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
            K2 = Kprofiles[upper, range(Kprofiles.shape[1])] * \
                weight_upper + \
                Kprofiles[lower, range(Kprofiles.shape[1])] * \
                (1-weight_upper)

            p = dt_mix * (2.0*K1 + dz*w)/(2.0*dz*dz)  # probability to rise
            q = dt_mix * (2.0*K2 - dz*w)/(2.0*dz*dz)  # probability to sink

            # check if probabilities are reasonable or wrong; which can happen if K is very high (K>0.1)
            wrong = p+q > 1.00002
            if wrong.sum() > 0:
                if i==0:
                    logging.info('WARNING! '+str(wrong.sum())+' elements have p+q>1; you might need a smaller mixing time step')
                else: 
                    logging.debug('WARNING! '+str(wrong.sum())+' elements have p+q>1;step')
		# fixing p and q by scaling them to match p+q=1:
		norm = p+q
		p[wrong] = p[wrong]/norm[wrong]
		q[wrong] = q[wrong]/norm[wrong]

            RandKick = np.random.random(self.num_elements_active())
            
            #up = np.where(RandKick < p)
            #down = np.where(RandKick > 1.0 - q)
            # Modified lines: do not mix particles which have resurfaced
            up = ((RandKick < p) & (self.elements.z < 0))
            down = ((RandKick > (1.0 - q)) & (self.elements.z < 0))

            self.elements.z[up] = self.elements.z[up] + dz
            self.elements.z[down] = self.elements.z[down] - dz

            #avoid that elements are above surface / below bottom
            self.resurface_elements(minimum_depth=-dz/10)
            bottom = np.where(self.elements.z < Zmin)
            self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

            # Call wave mixing, if implemented for this class
            self.wave_mixing(dt_mix)

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
        dz = 1.
        maxrange = -100

        def update(val):
            tindex = np.int(tslider.val)
            mainplot.cla()
            mainplot.grid()
            mainplot.hist(self.history['z'].T[tindex, :], bins=-maxrange/dz,
                          range=[maxrange, 0], orientation='horizontal')
            mainplot.set_ylim([maxrange, 0])
            mainplot.set_xlim([0, self.num_elements_total()*.1])
            mainplot.set_xlabel('number of particles')
            mainplot.set_ylabel('depth [m]')
            x_wind = self.history['x_wind'].T[tindex, :]
            y_wind = self.history['y_wind'].T[tindex, :]
            windspeed = np.mean(np.sqrt(x_wind**2 + y_wind**2))
            mainplot.set_title(str(self.get_time_array()[0][tindex]) +
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
        print x_wind
        y_wind = self.history['x_wind'].T[step, :]
        windspeed = np.mean(np.sqrt(x_wind**2 + y_wind**2))
        ax.set_title(str(self.get_time_array()[0][step]) +
                     '   Mean windspeed: %.1f m/s' % windspeed)
        if show is True:
            pyplot.show()
