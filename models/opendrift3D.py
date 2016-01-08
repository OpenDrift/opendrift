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
import logging
from datetime import datetime

from opendrift import OpenDriftSimulation
from elements import LagrangianArray


# Defining the oil element properties
class Lagrangian3DArray(LagrangianArray):
    """Extending LagrangianArray to a passive particle that moves in 3dimensions
    The Particle may be buoyant and/or subject to vertical mixing
    buoyant bahaviour is described by terminal velocity
    """

    variables = LagrangianArray.add_variables([
        ('terminal_velocity', {'dtype': np.float32,
                     'units': 'm/s',
                     'default': 0.})])

class OpenDrift3DSimulation(OpenDriftSimulation):
    """Open source buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway
        
        Generic module for particles that move in 3 dimensions 
        and may be to vertical turbulent mixing
        with the possibility for positive or negative buoyancy

        Particles could be pollution (e.g. oil droplets), plankton, nutrients or sediments

        Under construction.
    """

    def update_terminal_velocity(self):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables
        overload this function to create particle-specific behaviour
        """
        #self.elements.terminal_velocity = 0.00
        

    def vertical_advection(self, w_vel):
        """Move particles vertically according to vertical ocean current

            Vertical advection by ocean currents is small compared to termical velocity
            to do buoyancy and turbulent mixing
        """

    def vertical_mixing(self):
        """Mix particles vertically according to eddy diffusivity and buoyancy

            buoyancy is expressed as terminal velocity, which is the steady-state
            vertical velocity due to positive or negative buoyant behaviour.
            It is usually a function of particle density, diameter, and shape.
            
            Vertical particle displacemend du to turbulent mixing is calculated using 
            the "binned random walk scheme" (Thygessen and Aadlandsvik, 2007).
            The formulation of this scheme is copied from LADIM (IMR).
        """
        # if terminal_velocity is None:
        #     w = self.config['drift']['terminal_velocity']
        #else:
        #    w = terminal_velocity
        from models import eddydiffusivity

        dz = self.config['turbulentmixing']['verticalresolution']
        dt_mix = self.config['turbulentmixing']['timestep']

        # minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # place particle in center of bin
        self.elements.z = np.round(self.elements.z/dz)*dz

        #avoid that elements are above surface / below bottom
        surface = np.where(self.elements.z > dz/2.)
        self.elements.z[surface] = dz/2.
        bottom = np.where(self.elements.z < Zmin)
        self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

        # internal loop for fast time step of vertical mixing model
        # binned random walk needs faster time step compared to horizontal advection
        ntimes_mix = int(self.time_step.total_seconds()/dt_mix)
        logging.debug('Vertical mixing module:')
        logging.debug('turbulent diffusion with binned random walk scheme')
        logging.debug('using '+str(ntimes_mix)+' fast time steps of dt='+str(dt_mix)+'s')
        for i in range(0, ntimes_mix):
            # update terminal velocity according to environmental variables 
            self.update_terminal_velocity()
            w = self.elements.terminal_velocity

            # get vertical eddy diffusivity from environment or specific model
            if self.config['turbulentmixing']['diffusivitymodel'] == 'environment':
                K = self.environment.ocean_vertical_diffusivity
            else:
                K = getattr(eddydiffusivity, self.config['turbulentmixing']['diffusivitymodel'])(self)

            #K1 = K # K at current depth
            #K2 = K # K at depth z-dz

            p = dt_mix * (2.0*K + dz*w)/(2.0*dz*dz)  #probability to rise
            q = dt_mix * (2.0*K - dz*w)/(2.0*dz*dz)  #probability to sink

            # check if probabilities are reasonable or wrong
            wrong = p+q > 1.
            if wrong.sum() > 0:
                logging.debug('WARNING! some elements have p+q>1.')

            RandKick = np.random.random(self.num_elements_active())
            up   = np.where(RandKick < p)
            down = np.where(RandKick > 1.0 - q)

            self.elements.z[up]   = self.elements.z[up]   + dz
            self.elements.z[down] = self.elements.z[down] - dz

            #avoid that elements are above surface / below bottom
            surface = np.where(self.elements.z > dz/2.)
            self.elements.z[surface] = dz/2.
            bottom = np.where(self.elements.z < Zmin)
            self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

    def plot_vertical_distribution(self):
        """Function to plot vertical distribution of particles"""
        import matplotlib.pyplot as plt
        from matplotlib import dates
        dz=1.
        maxrange=-100
        hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        fig = plt.figure()
        ax = fig.gca()
        #ax.xaxis.set_major_formatter(hfmt)
        #plt.xticks(rotation='vertical')
        #times = [self.start_time + n*self.time_step
        #         for n in range(self.steps + 1)]
        #data = self.history[prop].T[0:len(times), :]
        plt.hist(self.elements.z, bins=-maxrange/dz, range=[maxrange,0], orientation='horizontal')
        plt.title('Vertical particle distribution')
        plt.xlabel('number of particles')
        plt.ylabel('depth [m]')
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        plt.show()


