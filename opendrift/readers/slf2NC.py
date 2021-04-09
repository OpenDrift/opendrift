###
# Exporting a Telemac 3D selafin file to an unstructured Netcdf file
# Following CF 1-8 and ugrid conventions
import numpy as np
import netCDF4 as nC
from datetime import datetime
import sys
from os import path, environ, sep
from data_manip.formats.selafin import Selafin

# Selafin variables
local='/home/julien/DATA/Marine plastics/PlasticsATBay/Simulations/oceano/NW_Scot/'#Directory
f=local+'r3d_tide-ES_real_gen_new_7.slf' #Input file
slf = Selafin(f)

steps= len(slf.tags['times'])
#timestep= slf.tags['times'][1]
starttime=datetime(slf.datetime[0],slf.datetime[1],slf.datetime[2],slf.datetime[3],slf.datetime[4])

d = nC.Dataset('somefile.nc', 'w')

d.setncattr('Conventions','CF-1.8') #https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.pdf
d.setncattr('title',slf.title)
d.setncattr('history',str(datetime.now())[:16])
d.setncattr('institution',"Plastic@bay - MTS-CFD") #Specifies where the original data was produced
d.setncattr('source',"Telemac 3D v8p2") # The  method  of  production  of  the  original  data.
d.setncattr('CoordinateSystem', 'Cartesian')
d.setncattr('CoordinateProjection', 'UTM30')

###make topology
node= d.createDimension('nnode',slf.npoin2)
node3= d.createDimension('nnode3',slf.npoin3)
nele= d.createDimension('nface',slf.nelem2)
maxnode=d.createDimension('nMaxface_nodes',slf.ndp2)
siglay= d.createDimension('layers',slf.nplan)
time= d.createDimension('time',None)#len(slf.tags['times'])
Two = d.createDimension('Two',2)

# time constrain
time= d.createVariable('time',np.float,('time'))
time.long_name='time'
time.standard_name='time'
time.axis="T"
time.units='seconds since '+str(starttime)
time.base_date=str(starttime)
time.calendar='standard'
time[:] = slf.tags['times']
### TO DO check conventions in netcdf for time .time_zone .format


## variables.
# mesh topology
mesh2= d.createVariable('mesh2', np.int)
mesh2.cf_role = "mesh_topology"
mesh2.long_name = "Topology data of 2D unstructured mesh"
mesh2.topology_dimension = 2
mesh2.node_coordinates = "x y"
mesh2.face_node_connectivity = "face_nodes"
mesh2.face_dimension = "nele"

face_nodes =  d.createVariable('face_nodes',np.int,('nface','nMaxface_nodes'), fill_value=999999)
face_nodes.cf_role = "face_node_connectivity"
face_nodes.long_name = "Maps every face to its corner nodes."
face_nodes.start_index = 0
face_nodes[:]=slf.ikle2

#creating a var of projection
proj = d.createVariable('crs',np.int)
proj.grid_mapping_name = 'transverse_mercator'
proj.scale_factor_at_central_meridian =0.9996
proj.longitude_of_central_meridian = 3
proj.latitude_of_projection_origin = 0
proj.false_easting =500000.
proj.false_northing=0

# mesh node coordinates (needs changed from crs)
x = d.createVariable('x', np.float, ('nnode'))
x.standard_name = "projection_x_coordinate"
x.axis="X"
x.long_name = "x-coordinate in Cartesian system"
x.units = "m"
x[:]=slf.meshx

y = d.createVariable('y', np.float, ('nnode'))
y.standard_name = "projection_y_coordinate"
y.axis="Y"
y.long_name = "y-coordinate in Cartesian system"
y.units = "m"
y[:]=slf.meshy

z= d.createVariable('z', np.float, ('nnode3'))
z.axis="Z"
z.long_name = "z-coordinate in Cartesian system"
z.units = "m"
for i in range(steps):
    z[i]=slf.get_values(i)[0]

# Vertical coordinate
layers= d.createVariable('layers', np.float,('layers'))
layers.standard_name = "ocean_sigma_coordinate"
layers.long_name = "sigma at layer"
layers.positive = "up"
layers.formula_terms = "sigma. layers eta. surface depth. depth"
layers[:]= np.linspace(0,1, slf.nplan) ###<= needs to be checked for FE

depth= d.createVariable('depth', np.float,('nnode'))
depth.standard_name = "sea_floor_depth_below_geoid"
depth.units = "m"
depth.positive = "up"
depth.mesh = "mesh2"
depth.location = "node"
depth.coordinates = "time x y"
for i in range(steps):
    depth[i]=slf.get_values(i)[0][:slf.npoin2]

surface= d.createVariable('surface', np.float,('nnode'))
surface.standard_name = "sea_surface_height_above_geoid"
surface.units = "m"
surface.mesh = "mesh2"
surface.location = "node"
surface.coordinates = "time x y"
for i in range(steps):
    surface[i]=slf.get_values(i)[0][-slf.npoin2:]

#creating variable currents X
u = d.createVariable('u','f4',('time','layers','nnode'))
u.standard_name='eastward_sea_water_velocity'
u.units='m s-1'
u.coordinates = 'time siglay x y'
for i in range(steps):
    u[i]=slf.get_values(i)[1]

#creating variable currents Y
v = d.createVariable('v','f4',('time','layers','nnode'))
v.standard_name='northward_sea_water_velocity'
v.units='m s-1'
v.coordinates = 'time siglay x y'
for i in range(steps):
    v[i]=slf.get_values(i)[2]

#creating variable currents Z
w = d.createVariable('w','f4',('time','layers','nnode'))
w.standard_name='upward_sea_water_velocity'
w.units='m s-1'
w.coordinates = 'time siglay x y'
for i in range(steps):
    w[i]=slf.get_values(i)[3]

d.close()
### Optional stuff not used
#mesh2.edge_node_connectivity = "edge_nodes"  // attribute required if variables will be defined on edges
#mesh2.edge_dimension = "nedge"
#mesh2.edge_coordinates = "edge_x edge_y"  // optional attribute (requires edge_node_connectivity)
#mesh2.face_coordinates = "face_x face_y"  // optional attribute
#mesh2.face_edge_connectivity = "face_edges"  // optional attribute (requires edge_node_connectivity)
#mesh2.face_face_connectivity = "face_links"  // optional attribute
#mesh2.edge_face_connectivity = "edge_face_links"  // optional attribute (requires edge_node_connectivity)
# integer edge_nodes(nedge, Two)
# edge_nodes.cf_role = "edge_node_connectivity"
# edge_nodes.long_name = "Maps every edge to the two nodes that it connects."
# edge_nodes.start_index = 1

# // Optional mesh topology variables
# integer face_edges(nface, nMaxface_nodes)
# face_edges.cf_role = "face_edge_connectivity"
# face_edges.long_name = "Maps every face to its edges."
# face_edges._FillValue = 999999
# face_edges.start_index = 1
# integer face_links(nface, nMaxface_nodes)
# face_links.cf_role = "face_face_connectivity"
# face_links.long_name = "neighbor faces for faces"
# face_links.start_index = 1
# face_links._FillValue = -999
# face_links.comment = "missing edges as well as missing neighbor faces are indicated using _FillValue"
# integer edge_face_links(nedge, Two)
# edge_face_links.cf_role = "edge_face_connectivity"
# edge_face_links.long_name = "neighbor faces for edges"
# edge_face_links.start_index = 1
# edge_face_links._FillValue = -999
# edge_face_links.comment = "missing neighbor faces are indicated using _FillValue"
# // Optional mesh face and edge coordinate variables
# double face_x(nface)
# face_x.standard_name = "longitude"
# face_x.long_name = "Characteristics longitude of 2D mesh face."
# face_x.units = "degrees_east"
# face_x.bounds = "face_xbnds"
# double face_y(nface)
# face_y.standard_name = "latitude"
# face_y.long_name = "Characteristics latitude of 2D mesh face."
# face_y.units = "degrees_north"
# face_y.bounds = "face_ybnds"
# double face_xbnds(nface,nMaxface_nodes)
# face_xbnds.standard_name = "longitude"
# face_xbnds.long_name = "Longitude bounds of 2D mesh face (i.e. corner coordinates)."
# face_xbnds.units = "degrees_east"
# face_xbnds._FillValue = 9.9692099683868690E36
# double face_ybnds(nface,nMaxface_nodes)
# face_ybnds.standard_name = "latitude"
# face_ybnds.long_name = "Latitude bounds of 2D mesh face (i.e. corner coordinates)."
# face_ybnds.units = "degrees_north"
# face_ybnds._FillValue = 9.9692099683868690E36
# double edge_x(nedge)
# edge_x.standard_name = "longitude"
# edge_x.long_name = "Characteristic longitude of 2D mesh edge (e.g. midpoint of the edge)."
# edge_x.units = "degrees_east"
# double edge_y(nedge)
# edge_y.standard_name = "latitude"
# edge_y.long_name = "Characteristic latitude of 2D mesh edge (e.g. midpoint of the edge)."
# edge_y.units = "degrees_north"
# // bounds variables for edges skipped
