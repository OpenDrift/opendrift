###
# Exporting a Telemac 3D selafin file to an unstructured Netcdf file
# Following CF 1-8 and ugrid conventions
import numpy as np
import netCDF4 as nC
from datetime import datetime, timedelta
import sys
from os import path, environ, sep
from data_manip.formats.selafin import Selafin
from scipy.interpolate import LinearNDInterpolator
#import scipy, scipy.spatial
from scipy.spatial import cKDTree
#import pyproj
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
        self.boundary = self._build_boundary_polygon_(
                        self.slf.meshx.compressed(),self.slf.meshx.compressed())
        self.timer_start("build index")
        logger.debug("building index of nodes..")
        # using selafin method (scipy)
        self.slf.set_kd_tree()

        # time management
        self.start_time=np.datetime64(datetime(self.slf.datetime[0],self.slf.datetime[1],
                self.slf.datetime[2],self.slf.datetime[3],self.slf.datetime[4]))
        self.time = self.start_time) + self.slf.tags['times'].astype('timedelta64[s]')
        self.end_time = self.time[-1]

        #self.variable_mapping = self.slf.varnames
        # for information


        self.timer_end("build index")
        self.timer_end("open dataset")

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
            t: np.datetime64
                age of the particle set
            variables: np.array(int)
                indexes of variables
        returns:
            vectors: np.array
            structured array of vectors u,v,w for each particles
        """
        ### nearest time tupple
        frames, duration = nearest_idx(self.time, t)

        ### 3 nearest nodes in 2D (prism face)
        _,iii = self.slf.tree.query(np.vstack((x,y)).T,k=1)
        # build depth ndarrays of each fibre
        niii = len(iii)
        idx_3D= np.tile(np.arange(self.slf.nplan),niii).reshape( \
                niii,self.slf.nplan)*iii[:,None]+iii[:,None]
        depths1=self.slf.get_variables_at(frames[0],[0])[0,idx_3D]
        depths2=self.slf.get_variables_at(frames[1],[0])[0,idx_3D]
        # locate the profile dimension
        pm=duration[0]*p1+duration[1]*p2
        # calculate distance from particles to nearest point depth
        # check function is fully vectorial
        # depth, _= nearest_idx(pm, z)
        idx_layer = np.abs(pm-z[:,None]).argmin(axis=1)
        # need optimisation:
        idx = idx3D[np.arange(niii), idx_layer]
        variables={}
        var_i=fvar(self, requested_variables):
        for i in var_i:
            vectors1=self.slf.get_variables_at(frames[0],[i])
            vectors2=self.slf.get_variables_at(frames[1],[i])
            variables[requested_variables[i]]=duration[0]*vectors1+duration[1]*vectors2
        return variables

    def cfvar(self, requested_variables):
        """
        setup a dictionary of equivalent names between Telemac and CF standard.
        return index of variable within selafin file to allow a query.
        """
        defaultvars =['ELEVATION Z     ', 'VELOCITY U      ',
                'VELOCITY V      ', 'VELOCITY W      ', 'NUX FOR VELOCITY',
                'NUY FOR VELOCITY', 'NUZ FOR VELOCITY', 'TURBULENT ENERGY',
                'DISSIPATION     ']
        standard_names=['altitude','eastward_sea_water_velocity',
        'northward_sea_water_velocity', 'upward_sea_water_velocity',None,None,
        None,'specific_turbulent_kinetic_energy_of_sea_water',
        'specific_turbulent_kinetic_energy_dissipation_in_sea_water' ]
        equiv={}
        for i in range(len(defaultvars)):
            if standard_names is not None:
                equiv[standard_names[i]]=defaultvars[i]
        vars = np.array(self.slf.varnames)
        eq_r = [equiv.get(v) for v in req]
        var_i =[]
        for i in range(len(eq_r)):
            var_i.append(np.where(vars==eq_r[i])[0][0])
        return var_i

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
        distance = array - value
        nearest = np.abs(dist).argmin()
        # test exact match and out of bounds
        if distance == 0 | np.all(distance>0)|np.all(distance>0):
            bounds = (array[nearest],None)
            dist = (1,None)
        else :
            if distance(nearest)>0:
                bounds = (nearest, nearest+1)
                prop= distance(nearest)/(distance(nearest)+distance(nearest+1))
                dist = (prop, 1-prop)
            else:
                bounds = (nearest-1, nearest)
                dist = (distance(nearest)+1, -1*distance(nearest))
        # test out of borders
        return bounds, dist
