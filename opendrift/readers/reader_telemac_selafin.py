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
# Copyright 2021, Julien Moreau, Plastic@Bay CIC (UK)

import numpy as np
from datetime import datetime, timedelta
import sys, os, logging

logger = logging.getLogger(__name__)
from opendrift.readers.basereader import BaseReader, UnstructuredReader

try:
    from data_manip.formats.selafin import Selafin
except ImportError:
    # Activating PYTEL
    try:
        pytel = os.path.join(os.environ['HOMETEL'], 'scripts', 'python3')
        if pytel not in sys.path:
            logger.warning('adding telemac python scripts to PYTHONPATH: %s' %
                           pytel)
            sys.path.append(pytel)

        from data_manip.formats.selafin import Selafin
    except (KeyError, ImportError):
        logger.error(
            'Telemac python scripts cannot be found. These are distributed together with the Telemac source code. This reader will not work.'
        )


class Reader(BaseReader, UnstructuredReader):
    """
    A reader for unstructured (irregularily gridded) `Telemac3D` files.
    Args:
        :param filename: A single Selafin file
        :type filename: string, required.
        :param name: Name of reader
        :type name: string, optional
        :param proj4: PROJ.4 string describing projection of data.
        :type proj4: string, optional
    .. seealso::
        py:mod:`opendrift.readers.basereader.unstructured`.
    """
    def __init__(self, filename=None, name=None, proj4=None, start_time=None):
        def vardic(vars_slf):
            """
            Match the selafin variables from Telemac 3D to the variables used in
            OpenDrift.
            """
            # Define all the variables used in OpenDrift as a dictionary
            # This is done to the best of our knowledge
            Vars_OD = {
                'VELOCITY U      ': 'x_sea_water_velocity',
                'VELOCITY V      ': 'y_sea_water_velocity',
                'VELOCITY W      ': 'upward_sea_water_velocity',
                'TURBULENT ENERGY': 'turbulent_kinetic_energy',
                'TEMPERATURE     ': 'sea_water_temperature',
                'SALINITY        ': 'sea_water_salinity',
                'NUZ FOR VELOCITY': 'ocean_vertical_diffusivity',
                'NUX FOR VELOCITY': 'horizontal_diffusivity',
                'ELEVATION Z     ': 'sea_floor_depth_below_sea_level',
            }

            No_OD_equiv = {
                'x_wind', 'y_wind', 'wind_speed',
                'sea_floor_depth_below_sea_level', 'wind_from_direction',
                'sea_ice_x_velocity', 'sea_ice_y_velocity',
                'sea_surface_wave_significant_height',
                'sea_surface_wave_stokes_drift_x_velocity',
                'sea_surface_wave_stokes_drift_y_velocity',
                'sea_surface_wave_period_at_variance_spectral_density_maximum',
                'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
                'sea_ice_area_fraction', 'surface_downward_x_stress',
                'surface_downward_y_stress', 'turbulent_generic_length_scale'
            }
            No_Telemac_equiv = {
                'NUY FOR VELOCITY',
                'DISSIPATION     ',
            }
            variables = []
            var_idx = []
            for i, var in enumerate(vars_slf):
                try:
                    variables.append(Vars_OD[var])
                    var_idx.append(i)
                except:
                    logger.info(
                        "Selafin variable {} has no equivalent in OpenDrift".
                        format(var))
            return np.array(variables), np.array(var_idx)

        self.name = name if name is not None else filename

        self.timer_start("open dataset")
        logger.info('Opening dataset: %s' % filename)
        self.slf = Selafin(filename)

        logger.info("File:\n{}\nTitle:\n{}".format(self.slf.file, \
                                self.slf.title))
        logger.info('Using projection: %s.' % proj4)
        self.proj4 = proj4

        logger.info('Reading 2D grid')

        # Run constructor of parent Reader class
        super().__init__()
        self.boundary = self._build_boundary_polygon_( \
                        self.slf.meshx,self.slf.meshy)

        self.timer_start("build index")
        logger.debug("building index of nodes...")
        # using selafin method (scipy)
        self.x, self.y = self.slf.meshx, self.slf.meshy

        logger.debug('nodes: %d' % len(self.x))

        # using scipy directly
        self.tree = self._build_ckdtree_(self.slf.meshx, self.slf.meshy)
        #bounds
        self.xmin,self.ymin,self.xmax, self.ymax= self.slf.meshx.min(),\
                self.slf.meshy.min(), self.slf.meshx.max(),self.slf.meshy.max()

        ### time management
        # Slf files have a wrong start time due to a formating error in Telemac
        if start_time is not None:
            if type(start_time) is datetime:
                self.start_time = start_time
            else:
                logger.warning(
                    "The specified start time is not a datetime object")
        else:
            logger.info("loading the datetime from the selafin file")
            self.start_time = datetime(self.slf.datetime[0],
                                       self.slf.datetime[1],
                                       self.slf.datetime[2],
                                       self.slf.datetime[3],
                                       self.slf.datetime[4])
        self.times = []
        for i in range(len(self.slf.tags['times'])):
            self.times.append(self.start_time +
                              timedelta(seconds=self.slf.tags['times'][i]))
        self.end_time = self.times[-1]
        self.altitude_ID=np.array(self.slf.varindex)[np.array( \
                            self.slf.varnames)=='ELEVATION Z     '].tolist()
        self.meshID=(np.arange(self.slf.nplan)[:,None] \
                     *self.slf.npoin2).astype(int)
        self.variables, self.var_idx = vardic(self.slf.varnames)

        self.timer_end("build index")
        self.timer_end("open dataset")

    def plot_mesh(self):
        """
        Plot the grid mesh. Does not automatically show the figure.
        """
        title = 'Unstructured grid: %s\n%s' % (self.name, self.proj)
        from importlib.util import find_spec
        if find_spec("pyvista") is not None:
            import pyvista as pv
            cells = np.hstack(
                ((np.ones(len(self.slf.ikle2), dtype=np.int32) * 3)[:, None],
                 self.slf.ikle2))
            points = np.vstack((self.slf.meshx, self.slf.meshy,
                                np.zeros(len(self.slf.meshx)))).T
            u = pv.PolyData(points, cells)
            plotter = pv.Plotter()
            plotter.add_mesh(u, show_edges=True)
            plotter.show_bounds(mesh=u)
            plotter.view_xy()
            plotter.show(title=title, window_size=[800, 640])
        else:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.scatter(self.slf.meshx,
                        self.slf.meshy,
                        marker='x',
                        color='blue',
                        label='nodes')
            x, y = getattr(self.boundary, 'context').exterior.xy
            plt.plot(x, y, color='green', label='boundary')

            plt.legend()
            plt.title(title)
            plt.xlabel('x [m]')
            plt.ylabel('y [m]')

    def get_variables(self,
                      requested_variables,
                      time=None,
                      x=None,
                      y=None,
                      z=None):
        """
        - Query variables based on the particle coordinates x, y, z
        - find the nearest node in the KD tree
        - extract the z array corresponding.
        - extract the index of the node within the 3D mesh
        - extract the variables at the point
        Args:
            x,y,z: np.arrays(float)
                3D coordinates of the particles
            time: np.datetime64
                age of the particle set
            variables: np.array(int)
                indexes of variables
        Returns:
            variables: dictionary of numpy arrays
        """
        ### nearest time tupple
        frames, duration=self.__nearest_idx__(np.array(self.times).astype( \
                         'datetime64[s]'),np.datetime64(time))
        ### nearest node in 2D
        iii = self.__nearest_ckdtree__(self.tree, x, y)
        # build depth ndarrays of each fibre
        niii = len(iii)
        idx_3D = self.meshID + iii
        # locate the profile dimension
        pm = self.__extractslf__(frames, duration, self.altitude_ID, idx_3D)[0]
        # calculate distance from particles to nearest point altitude
        idx_layer = np.abs(pm - z).argmin(axis=0)
        indices = self.variables == np.array(requested_variables)[:, None]
        Idxs = []
        for line in indices:
            Idxs.append(self.var_idx[line])
        idx_nodes = idx_3D[idx_layer, np.arange(len(idx_layer))]
        vectors = self.__extractslf__(frames, duration, Idxs,
                                      idx_nodes.ravel())
        vars = {}
        for i in range(len(requested_variables)):
            if requested_variables[i] == 'sea_floor_depth_below_sea_level':
                bottom = self.__extractslf__(frames, duration,
                                             self.altitude_ID, idx_3D[0])
                surface = self.__extractslf__(frames, duration,
                                              self.altitude_ID, idx_3D[-1])
                vars['sea_floor_depth_below_sea_level'] = surface - bottom
            else:
                vars[requested_variables[i]] = vectors[i]
        return vars

    @staticmethod
    def __nearest_idx__(array, value):
        """
        we are looking for a tuple describing where the sample is and at which
        distance. So we can calculate the FE solution of the variable value.
        input:
            array: a 1D numpy array
            monotonic array
        output:
            bounds: a tupple
            tupple describing the bounding indexes
            dist = tupple
            the distance between the sample and the index values
            so that dist[0]+dist[1]=0
        """
        distance = (array - value).astype(float)
        nearest = np.argsort(abs(distance))[:2]
        # test exact match and out of bounds
        if (distance == 0).any() | (distance > 0).all() | (distance < 0).all():
            bounds = (nearest[0], None)
            dist = (1, 0)
        else:
            bounds = nearest
            if distance[nearest[0]] == distance[nearest[1]] * -1:
                dist = (.5, .5)
            else:
                prop= abs(distance[nearest[0]])/ \
                   (abs(distance[nearest[0]])+abs(distance[nearest[1]]))
                dist = (1 - prop, prop)
        return bounds, dist

    def __extractslf__(self, frames, duration, index_var, index_nodes):
        """
        extract variables from slf files
        index_var must be a list of integer
        """
        vector1 = self.slf.get_variables_at(frames[0], index_var)[:,
                                                                  index_nodes]
        if frames[1] is not None:
            vector2 = self.slf.get_variables_at(frames[1],
                                                index_var)[:, index_nodes]
        else:
            vector2 = 0

        vector = duration[0] * vector1 + duration[1] * vector2
        vector1, vector2 = None, None
        return vector
