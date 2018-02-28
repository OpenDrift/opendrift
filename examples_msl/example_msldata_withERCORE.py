#!/usr/bin/env python
import numpy,datetime,sys,os,time
sys.path.append(os.path.join('..'))
from ercore import ERcore,dt2ncep
from ercore.fields import *
from ercore.materials import BuoyantTracer
from ercore.materials import PassiveTracer
from ercore.shoreline import Boundary
from ercore.shoreline import Shoreline

#
# The script aims at running ercore for the exact same conditions as in example_grid_msldata.py
#
# To start ERcore docker :
# docker run -it -v /home/simon:/home/simon metocean/ercore_opensrc
# run folder is here : /home/simon/COMPARISON_ERCORE_OPENDRIFT/

site_names = ['test']
x = [174.5133]
y=  [-41.2348] 
# path of folder with all necessary run files 

# STICKERS--------------------------------------------------------------------------------------------------

bnd=Boundary('bnd',[168.0,177.46,-42.84,-37.9])
shoreline = Shoreline('shoreline', '/home/simon/0387_OMV_drillcut_oilspill/drill_cuttings/' + 'CNZ_coast.bna')
# MOVERS----------------------------------------------------------
fname_flow = 'cnz_surf_res200401.nc' # just for testing - need to update
current = GriddedMover('current',['uso','vso'],file = fname_flow)

# OUTPUTS----------------------------------------------------------------
output_folder = './outputs/'
if not(os.path.exists(output_folder)):
    os.mkdir(output_folder)
#-----------------------------------------------------------------------
# RUN TIMING -------------------------------------------------------------------------------
sedclass_name = ['passive']
pdt = 1800
tstart = datetime.datetime(2004,1,1,0,0,0)
tend = datetime.datetime(2004,1,4,0,0,0)
dt_out  =  1800
npart_per_release = 1000 # used to be one but these sites are shallow..

ercore = ERcore(geod = True,tout = dt_out)


material = PassiveTracer(id = 'test',
                                  nbuff = npart_per_release,
                                  movers = [current], #current3d
                                  stickers = [bnd,shoreline],
                                  unstick = [0,0],
                                  tstart = tstart, # single release
                                  tend = tstart,
                                  reln = int(npart_per_release),
                                  # tstep_release = dt_release_hours,
                                  P0 = [x[0],y[0],-1],
                                  save_initial_positions = True)

ercore.materials = [material] 
current.reset()             # Reset Time cursor in flow files
ercore.run(tstart,tend,pdt) 
