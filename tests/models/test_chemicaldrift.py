#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
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

from datetime import datetime, timedelta
import numpy as np
import xarray as xr
import opendrift
from opendrift.models.chemicaldrift import ChemicalDrift
from opendrift.readers.reader_constant import Reader as ConstantReader

"""Tests for ChemicalDrift module."""

def test_chemicaldrift_partitioning_organics():
    """ Test dynamic partitioning for organic compunds with different values of
        fOC (organic carbon fraction) in suspended particles and KOW.
        High SPM concentration, 80 g/m3."""

    days=200
    ntraj=200

    # If one parameter is given as a list, different instances of
    # ChemicalDrift will be created
    configs = [
    [ "seed:LMM_fraction"                       , 1                     ],
    [ "seed:particle_fraction"                  , 0                     ],
    [ "chemical:species:Sediment_reversible"    , True                  ],
    [ "chemical:dynamic_partitioning"           , True                  ],
    [ "chemical:transfer_setup"                 , 'organics'            ],
    [ "chemical:transformations:fOC_SPM"        , [0.01, 0.05]          ],
    [ "chemical:transformations:LogKOW"         , [4, 6.618, 7]         ],
    [ "chemical:transformations:DeltaH_KOC_Sed" , 0                     ],
    [ "chemical:transformations:Setchenow"      , 0                     ],
    [ "chemical:particle_diameter"              , 0                     ],
    [ "chemical:particle_diameter_uncertainty"  , 0                     ],
    [ "drift:vertical_mixing"                   , False                 ],
    [ "environment:constant:land_binary_mask"   , 0                     ],
    [ "general:use_auto_landmask"               , False                 ],
    ]

    # possible to specify lists of different constant readers values
    # that will be used in different longitude intervals
    readers = [
    [ "spm"                                     , 80                    ],
    [ "x_sea_water_velocity"                    , 0                     ],
    [ "y_sea_water_velocity"                    , 0                     ]]

    # Build list of ChemicalDrift objecs with different parameters
    o=list()
    for i in range(len(configs)):
        if type(configs[i][1])==list:
            for j in range(len(configs[i][1])):
                o.append(ChemicalDrift(loglevel=20))
                #o.append(ChemicalDrift(loglevel=20,logfile="logfile_"+str(j)+".txt"))
                o[-1].set_config(configs[i][0],configs[i][1][j])

    if o==[]:
        j=0
        o.append(ChemicalDrift(loglevel=0))

    for i in range(len(configs)):
        if type(configs[i][1])!=list:
            for j in range(len(o)):
                o[j].set_config(configs[i][0],configs[i][1])

    # seed in the middle of the Altantic Ocean
    lonmin=-50
    lonmax=-50+1
    lat=25

    r=list()
    for i in range(len(readers)):
        if type(readers[i][1])==list:
            for j in range(len(readers[i][1])):
                r.append(ConstantReader({readers[i][0]: readers[i][1][j]}))
                r[-1].xmin = lonmin + j
                r[-1].xmax = lonmax = lonmin + j + 1
        else:
            r.append(ConstantReader({readers[i][0]: readers[i][1]}))

    for j in range(len(o)):
        o[j].add_reader(r)

    for j in range(len(o)):
        iniz=np.random.rand(ntraj) * -10. # seeding the chemicals in the upper 10m
        for k in range(lonmax-lonmin):
            o[j].seed_elements(lon=lonmin+0.5+k, lat=lat, time=datetime(2015, 1, 1, 0), z=iniz, number=ntraj, radius=20000)

        o[j].run(duration=timedelta(hours=24*days),time_step=60*60*24*100,time_step_output=60*60*24*100)

    # Show actual computed values (use --show-capture=stdout option in pytest)
    for j in range(len(o)):
        print("run ", j)
        print("fOC : ", o[j].get_config("chemical:transformations:fOC_SPM"))
        print("LogKOW : ", o[j].get_config("chemical:transformations:LogKOW"))
        print("#SPM: ", sum(o[j].elements.specie==2 ))

    assert sum(o[0].elements.specie==2) == 1
    assert sum(o[1].elements.specie==2) == 1
    assert sum(o[2].elements.specie==2) == 0
    assert sum(o[3].elements.specie==2) == 18
    assert sum(o[4].elements.specie==2) == 64

def test_chemicaldrift_degradation():
    """ Test degradation for organic compounds with different values of
        t12 (half life) """

    days=200
    ntraj=200

    # If one parameter is given as a list, different instances of
    # ChemicalDrift will be created
    configs = [
    [ "seed:LMM_fraction"                       , 1                     ],
    [ "seed:particle_fraction"                  , 0                     ],
    [ "chemical:species:Sediment_reversible"    , True                  ],
    [ "chemical:dynamic_partitioning"           , True                  ],
    [ "chemical:transfer_setup"                 , 'organics'            ],
    [ "chemical:transformations:LogKOW"         , 6                     ],
    [ "chemical:transformations:fOC_SPM"        , 0.1                   ],
    [ "chemical:transformations:LogKOW"         , 2                     ],
    [ "chemical:transformations:t12_W_tot"      , [500, 1000]           ],
    [ "chemical:transformations:DeltaH_KOC_Sed" , 0                     ],
    [ "chemical:transformations:Setchenow"      , 0                     ],
    [ "chemical:particle_diameter"              , 0                     ],
    [ "chemical:particle_diameter_uncertainty"  , 0                     ],
    [ "drift:vertical_mixing"                   , False                 ],
    [ "environment:constant:land_binary_mask"   , 0                     ],
    [ "general:use_auto_landmask"               , False                 ],
    [ "chemical:transformations:degradation"    , True                  ],
    [ "chemical:transformations:volatilization" , False                 ],
    ]

    # possible to specify lists of different constant readers values
    # that will be used in different longitude intervals
    readers = [
    [ "spm"                                     , 0.1                   ],
    [ "x_sea_water_velocity"                    , 0                     ],
    [ "y_sea_water_velocity"                    , 0                     ]]

    # Build list of ChemicalDrift objecs with different parameters
    o=list()
    for i in range(len(configs)):
        if type(configs[i][1])==list:
            for j in range(len(configs[i][1])):
                o.append(ChemicalDrift(loglevel=0))
                #o.append(ChemicalDrift(loglevel=0,logfile="logfile_"+str(j)+".txt"))
                o[-1].set_config(configs[i][0],configs[i][1][j])

    if o==[]:
        j=0
        o.append(ChemicalDrift(loglevel=0))

    for i in range(len(configs)):
        if type(configs[i][1])!=list:
            for j in range(len(o)):
                o[j].set_config(configs[i][0],configs[i][1])

    # seed in the middle of the Atlantic Ocean
    lonmin=-50
    lonmax=-50+1
    lat=25

    r=list()
    for i in range(len(readers)):
        if type(readers[i][1])==list:
            for j in range(len(readers[i][1])):
                r.append(ConstantReader({readers[i][0]: readers[i][1][j]}))
                r[-1].xmin = lonmin + j
                r[-1].xmax = lonmax = lonmin + j + 1
        else:
            r.append(ConstantReader({readers[i][0]: readers[i][1]}))

    for j in range(len(o)):
        o[j].add_reader(r)

    for j in range(len(o)):
        iniz=np.random.rand(ntraj) * -10. # seeding the chemicals in the upper 10m
        for k in range(lonmax-lonmin):
            o[j].seed_elements(lon=lonmin+0.5+k, lat=lat, time=datetime(2015, 1, 1, 0), z=iniz, number=ntraj, radius=20000)

        o[j].run(duration=timedelta(hours=24*days),time_step=60*60*24*100,time_step_output=60*60*24*100)

    # Show actual computed values (use --show-capture=stdout option in pytest)
    for j in range(len(o)):
        print("run ", j)
        print("t12 : ", o[j].get_config("chemical:transformations:t12_W_tot"))
        print("mass: ", sum(o[j].result.mass.values ))
        print("degraded: ", sum(o[j].result.mass_degraded.values ))
        print("total: ", sum(o[j].result.mass.values) + sum(o[j].result.mass_degraded.values ))
        print("dissolved: ", sum(o[j].result.specie.values==0 ))        

    #assert 0

    # Check t12 500
    assert abs(sum(o[0].result.mass.values)[-1] - 20771.232) < 0.01
    assert abs(sum(o[0].result.mass_degraded.values)[-1] - 179228.38) < 0.01

    # Check mass conservation residual + degraded = total
    assert abs(sum(o[0].result.mass.values)[-1] + sum(o[0].result.mass_degraded.values)[-1] - 200000) < 0.5

    # Check mass conservation residual + degraded = initial = ( 1000 ug default value )
    assert np.all(o[0].result.mass.values + o[0].result.mass_degraded.values == 1000)

    # Check t12 1000
    # Residual mass larger, degraded mass smaller
    assert abs(sum(o[1].result.mass.values)[-1] - 64513.277) < 0.01
    assert abs(sum(o[1].result.mass_degraded.values)[-1] - 135486.4) < 0.01

    # Check mass conservation residual + degraded = initial = ( 1000 ug default value )
    assert np.all(o[1].result.mass.values + o[1].result.mass_degraded.values == 1000)

def test_chemicaldrift_volatilization():
    """ Test volatilization for different values of wind speed. """

    days=2
    ntraj=2

    # If one parameter is given as a list, different instances of
    # ChemicalDrift will be created
    configs = [
    [ "seed:LMM_fraction"                       , 1                     ],
    [ "seed:particle_fraction"                  , 0                     ],
    [ "chemical:species:Sediment_reversible"    , True                  ],
    [ "chemical:dynamic_partitioning"           , True                  ],
    [ "chemical:transfer_setup"                 , 'organics'            ],
    [ "chemical:transformations:fOC_SPM"        , 0.1                   ],
    [ "chemical:transformations:LogKOW"         , 2                     ],
    [ "chemical:transformations:DeltaH_KOC_Sed" , 0                     ],
    [ "chemical:transformations:Setchenow"      , 0                     ],
    [ "chemical:particle_diameter"              , 0                     ],
    [ "chemical:particle_diameter_uncertainty"  , 0                     ],
    [ "drift:vertical_mixing"                   , False                 ],
    [ "environment:constant:land_binary_mask"   , 0                     ],
    [ "general:use_auto_landmask"               , False                 ],
    [ "chemical:transformations:degradation"    , False                 ],
    [ "chemical:transformations:volatilization" , True                  ],
    ]

    # possible to specify lists of different constant readers values
    # that will be used in different longitude intervals
    readers = [
    [ "spm"                                     , 0.1                   ],
    [ "x_sea_water_velocity"                    , 0                     ],
    [ "y_sea_water_velocity"                    , 0                     ],
    [ "x_wind"                                  , [5,10]                ],
    [ "y_wind"                                  , 0                     ],
    ]

    # Build list of ChemicalDrift objecs with different parameters
    o=list()
    for i in range(len(configs)):
        if type(configs[i][1])==list:
            for j in range(len(configs[i][1])):
                o.append(ChemicalDrift(loglevel=0))
                #o.append(ChemicalDrift(loglevel=0,logfile="logfile_"+str(j)+".txt"))
                o[-1].set_config(configs[i][0],configs[i][1][j])

    if o==[]:
        j=0
        o.append(ChemicalDrift(loglevel=0))

    for i in range(len(configs)):
        if type(configs[i][1])!=list:
            for j in range(len(o)):
                o[j].set_config(configs[i][0],configs[i][1])

    # seed in the middle of the Atlantic Ocean
    lonmin=0
    lonmax=0+1
    lat=0

    r=list()
    for i in range(len(readers)):
        if type(readers[i][1])==list:
            for j in range(len(readers[i][1])):
                r.append(ConstantReader({readers[i][0]: readers[i][1][j]}))
                r[-1].xmin = lonmin + j
                r[-1].xmax = lonmax = lonmin + j + 1
        else:
            r.append(ConstantReader({readers[i][0]: readers[i][1]}))

    for j in range(len(o)):
        o[j].add_reader(r)

    for j in range(len(o)):
        iniz=np.random.rand(ntraj) * -10. # seeding the chemicals in the upper 10m
        for k in range(lonmax-lonmin):
            o[j].seed_elements(lon=lonmin+0.5+k, lat=lat, time=datetime(2015, 1, 1, 0), z=iniz, number=ntraj, radius=1000., mass=1.)

        o[j].run(duration=timedelta(hours=24*days),time_step=60*60*24*1,time_step_output=60*60*24*1)

    # Show actual computed values (use --show-capture=stdout option in pytest)
    for j in range(len(o)):
        print("run ", j)
        print("mass: ", o[j].result.mass.values )
        print("volatilized : ", o[j].result.mass_volatilized.values )
        print("total: ", o[j].result.mass.values + o[j].result.mass_volatilized.values)
        print("wind: ", o[j].result.x_wind.values )

    #assert 0
    
    # Check trajectory 0, wind 5 m/s
    assert abs( o[0].result.mass_volatilized.sel(trajectory=0,time='2015-01-03') - 0.03400866) < 1e-5

    # Check trajectory 2, wind 10 m/2, stronger volatilization
    assert abs( o[0].result.mass_volatilized.sel(trajectory=2,time='2015-01-03') - 0.11619269) < 1e-5
    
    # Check mass conservation mass residual + mass volatilized = initial = 1
    # for all trajectories
    assert np.all(o[0].result.mass.values + o[0].result.mass_volatilized.values == 1 )


def test_chemicaldrift_write_concentration_map():
    import os

    def delete_file(file_path):
        if os.path.exists(file_path):
            os.remove(file_path)

    o = ChemicalDrift(loglevel=0, seed=0)

    mixed_layer = ConstantReader({'ocean_mixed_layer_thickness': 40})
    u = ConstantReader({'x_sea_water_velocity' : 0})
    v = ConstantReader({'y_sea_water_velocity' : 0})
    o.add_reader([u,v,mixed_layer])

    o.set_config('environment:constant:land_binary_mask',0)
    o.set_config('general:use_auto_landmask', False)
    o.set_config("drift:horizontal_diffusivity", 10)
    o.set_config("seed:LMM_fraction",1)
    o.set_config("seed:particle_fraction",0)
    o.set_config("chemical:dynamic_partitioning", False)
    o.init_chemical_compound("Copper")

    # Seeding 500 lagrangian elements each representign 2mg og target chemical

    time = datetime(1976, 4, 14, 0, 0)

    ntraj=500
    mass=2e3
    import numpy as np
    iniz=np.random.rand(ntraj) * -10. # seeding the chemicals in the upper 10m

    o.seed_elements(0., 0., z=iniz, radius=2000,number=ntraj,time=time, mass=mass) 

    # Running model
    #
    # if we want to test that it works also after saving the output file and reloading
    #
    # o.run(steps=10, time_step=1800, time_step_output=1800, outfile="tmp_test_chemicaldrift_output.nc",export_buffer_length=5)
    # del(o)
    # import opendrift
    # o=opendrift.open("tmp_test_chemicaldrift_output.nc")

    o.run(steps=10, time_step=60*60*24, time_step_output=60*60*24)

    # Define the grid
    lat = np.linspace(-2., 2., 100)
    lon = np.linspace(-2., 2., 200)

    # Create constant depth data (e.g., 1500 meters everywhere)
    depth = np.full((len(lat), len(lon)), 1500.0)

    # Create the dataset
    ds_depth = xr.Dataset(
        {
            "sea_floor_depth_below_sea_level": (["lat", "lon"], depth)
        },
        coords={
            "lat": lat,
            "lon": lon,
        },
    )

    # Add attributes
    ds_depth["sea_floor_depth_below_sea_level"].attrs["units"] = "m"
    ds_depth["sea_floor_depth_below_sea_level"].attrs["standard_name"] = "sea_floor_depth_below_sea_level"

    bathymetry_file="tmp_test_chemicaldrift_bathymetry.nc"
    delete_file(bathymetry_file)
    ds_depth.to_netcdf(bathymetry_file)
    
    concentration_file="tmp_test_chemicaldrift_concentgration.nc"
    delete_file(concentration_file)
    o.write_netcdf_chemical_density_map(concentration_file, pixelsize_m=1000.,
                                             zlevels=None,
                                             mass_unit='ug',
                                             horizontal_smoothing=False,
                                             smoothing_cells=1,
                                             time_avg_conc=True,
                                             deltat=24, # hours
                                             llcrnrlon=-1., llcrnrlat=-1.,
                                             urcrnrlon=1., urcrnrlat=1.,
                                             reader_sea_depth='./'+bathymetry_file)
    
    # Check that the concentration file is created

    assert os.path.exists(concentration_file)
    ds=xr.open_dataset(concentration_file)

    # Check for each timestep that the mass in conserved
    # total mass = sum_for_all_cells of concentration x volume )

    for time in ds.concentration_avg.avg_time:
        assert np.sum(ds.concentration_avg.sel(avg_time=time).sum(dim='specie')*ds.volume) == (ntraj*mass)
    
    # Check that, due to diffusion, the number of cells with concentration > 0 increase at each time step

    cells=[]
    for time in ds.avg_time:
        cells.append(np.sum(ds.sel(avg_time=time,depth=0,specie=0).concentration_avg.data>0))
    assert all(x < y for x, y in zip(cells, cells[1:]))

    delete_file(concentration_file)
    delete_file(bathymetry_file)
