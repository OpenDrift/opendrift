import numpy as np
from datetime import datetime, timedelta
import sys
from os import path, environ, sep
from scipy.spatial import cKDTree
#Activating PYTEL
dir = environ['HOMETEL']+'{}scripts{}python3'.format(sep,sep)
sys.path.append(path.join( path.dirname(sys.argv[0]),dir))#pytel path use environment variables
from data_manip.formats.selafin import Selafin
import logging
logger = logging.getLogger(__name__)

from opendrift.readers.basereader import BaseReader, UnstructuredReader


class Reader(BaseReader, UnstructuredReader):
    """
    A reader for unstructured (irregularily gridded) `Telemac3D files.

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
    # node_variables = ['z','salinity','temperature',
    #     'eastward_sea_water_velocity',
    #     'northward_sea_water_velocity',
    #     'upward_sea_water_velocity'
    def __init__(self, filename=None, name=None, proj4=None):
        try:
            filestr = str(filename)
        except Exception as e:
                logger.info(e)
        ### that's not a needed test?
        if name is None:
            self.name = filestr
        else:
            self.name = name
        self.workers = -1
        self.timer_start("open dataset")
        logger.info('Opening dataset: ' + filestr)

        logger.info('Opening file with Dataset')
        self.slf = Selafin(filename)

        logger.info("File:\n{}\nTitle:\n{}".format(self.slf.file, \
                                self.slf.title))
        try:
            str(proj4)
        except Exception as e:
                logger.info(e)
        logger.info('Using custom projection: %s.' % proj4)
        self.proj4 = proj4

        logger.info('Reading 2D grid')
        # Run constructor of parent Reader class
        super().__init__()
        self.boundary = self._build_boundary_polygon_( \
                        self.slf.meshx,self.slf.meshy)
        self.timer_start("build index")
        logger.debug("building index of nodes..")
        # using selafin method (scipy)
        #self.slf.tree= set_kd_tree()
        # using scipy directly
        self.tree = self._build_ckdtree_(self.slf.meshx, self.slf.meshy)
        #bounds
        self.xmin,self.ymin,self.xmax, self.ymax= self.slf.meshx.min(),\
                self.slf.meshy.min(), self.slf.meshx.max(),self.slf.meshy.max()

        # time management
        self.start_time=datetime(self.slf.datetime[0],self.slf.datetime[1],
                self.slf.datetime[2],self.slf.datetime[3],self.slf.datetime[4])
        # self.start_time=np.datetime64(datetime(self.slf.datetime[0],self.slf.datetime[1],
        #         self.slf.datetime[2],self.slf.datetime[3],self.slf.datetime[4]))
        # self.time = self.start_time + self.slf.tags['times'].astype('timedelta64[s]')
        self.times= []
        for i in range(len(self.slf.tags['times'])):
            self.times.append(self.start_time+timedelta(seconds= self.slf.tags['times'][i]))
        self.end_time = self.times[-1]
        # sorry hardcoded
        self.variables=['altitude','x_sea_water_velocity',
        'y_sea_water_velocity', 'upward_sea_water_velocity',
        'specific_turbulent_kinetic_energy',
        'specific_turbulent_kinetic_energy_dissipation_in_sea_water' ]


        self.timer_end("build index")
        self.timer_end("open dataset")

    def plot_mesh(self):
        """
        Plot the grid mesh. Does not automatically show the figure.
        """
        #### needs to read the architecture of the nodes
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(self.slf.meshx, self.slf.meshy, marker='x', color='blue', label='nodes')
        x, y = getattr(self.boundary, 'context').exterior.xy
        plt.plot(x, y, color='green', label='boundary')

        plt.legend()
        plt.title('Unstructured grid: %s\n%s' % (self.name, self.proj))
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')

    def get_variables(self,
                      requested_variables,
                      time=None,
                      x=None,
                      y=None,
                      z=None):
        """
        Querry variables based on the particle coordinates x, y, z
        find the nearest node in the KD tree
        extract the z array corresponding.
        extract the index of the node within the 3D mesh
        extract the variables at the point
        input:
            x,y,z: np.arrays(float)
                3D coordinates of the particles
            time: np.datetime64
                age of the particle set
            variables: np.array(int)
                indexes of variables
        returns:
            variables: dictionary of numpy arrays
        """
        def cfvar(self, requested_variables):
            """
            setup a dictionary of equivalent names between Telemac and CF standard.
            return index of variable within selafin file to allow a query.
            """
            defaultvars =['ELEVATION Z     ', 'VELOCITY U      ',
            'VELOCITY V      ', 'VELOCITY W      ', 'NUX FOR VELOCITY',
            'NUY FOR VELOCITY', 'NUZ FOR VELOCITY', 'TURBULENT ENERGY',
            'DISSIPATION     ']
            #eastward northward removed for x y
            standard_names=['altitude','x_sea_water_velocity',
            'y_sea_water_velocity', 'upward_sea_water_velocity',None,None,
            None,'specific_turbulent_kinetic_energy',
            'specific_turbulent_kinetic_energy_dissipation_in_sea_water' ]
            equiv={}
            vars = np.array(self.slf.varnames)
            for i in range(len(defaultvars)):
                if standard_names is not None:
                    equiv[standard_names[i]]=defaultvars[i]
                    eq_r = [equiv.get(v) for v in requested_variables]
            var_i =[]
            for i in range(len(eq_r)):
                    var_i.append(np.where(vars==eq_r[i])[0][0])
            return var_i

        def nearest_idx(array, value):
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
                distance = (array - value).astype(np.float)
                nearest =  np.argsort(abs(distance))[:2]
                # test exact match and out of bounds
                if (distance == 0).any() | (distance>0).all()|(distance<0).all():
                    bounds = (nearest[0],None)
                    dist = (1,0)
                else :
                    prop= abs(distance[nearest[0]]/(distance[nearest[0]]+distance[nearest[1]]))
                    dist = (1-prop, prop)
                    bounds=nearest
                return bounds, dist

        ### nearest time tupple
        frames, duration=nearest_idx(np.array(self.times).astype('datetime64[s]'),np.datetime64(time))
        ### nearest node in 2D
        # will have to be tested for scipy>1.6.0 upgrades n_jobs=workers
        # scipy 1.3.3 fails as there is no argument n_jobs
        # will have to be improved for a real finite element solution (k=3 pts).
        _,iii = self.tree.query(np.vstack((x,y)).T,k=1, workers=self.workers)
        # build depth ndarrays of each fibre
        niii = len(iii)
        idx_3D = np.arange(self.slf.nplan).reshape(self.slf.nplan,1)*self.slf.npoin2+iii
        depths1=self.slf.get_variables_at(frames[0],[0])[0,idx_3D]
        if frames[1] is not None:
            depths2=self.slf.get_variables_at(frames[1],[0])[0,idx_3D]
        else:
            depths2=0
        # locate the profile dimension
        pm=duration[0]*depths1+duration[1]*depths2
        # calculate distance from particles to nearest point depth
        # check function is fully vectorial
        # depth, _= nearest_idx(pm, z)
        idx_layer = np.abs(pm.flatten()-z[:,None]).argmin(axis=1)
        # need optimisation:
        #idx = idx_3D[np.arange(niii), idx_layer]
        variables={}
        var_i=cfvar(self, requested_variables)
        for i in range(len(var_i)):
            vectors1=self.slf.get_variables_at(frames[0],[var_i[i]])[0,idx_layer].ravel()
            if frames[1] is not None:
                vectors2=self.slf.get_variables_at(frames[1],[var_i[i]])[0,idx_layer].ravel()
            else:
                vectors2=0

            variables[requested_variables[i]]=duration[0]*vectors1+duration[1]*vectors2
        return variables


    def filter_points(self, indexes):
        """
        Filter points that are not within the grid.
        use finite element method to evaluate properties
        ~~~ To be continued ~~~
        """
        # test they correspond to faces:
        ifaces= np.where((np.sort(iii, axis=1)[:,None] ==np.sort(self.slf.ikle2, axis=1)).all(-1).any(-1))[0]
        # in the future
        # extract profile from the 2 frames bounding t.
        # z is always the variable idx 0
        p1= self.slf.get_variables_at(frames[0], [0])[0,self.slf.ikle3[ifaces]]
        p2= self.slf.get_variables_at(frames[1], [0])[0,self.slf.ikle3[ifaces]]
        x1=slef.slf.meshx[ifaces]
        y1=slef.slf.meshy[ifaces]
