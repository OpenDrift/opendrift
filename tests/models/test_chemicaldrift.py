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

    # Check t12 1000
    # Residual mass larger, degraded mass smaller
    assert abs(sum(o[1].result.mass.values)[-1] - 64513.277) < 0.01
    assert abs(sum(o[1].result.mass_degraded.values)[-1] - 135486.4) < 0.01

    # Check mass conservation residual + degraded = total
    assert abs(sum(o[1].result.mass.values)[-1] + sum(o[1].result.mass_degraded.values)[-1] - 200000) < 0.5

    # TODO: Check why conserved mass is not exactly 200000 - some rounding issues?
    # Or is it because of deactivated elements that are not counted?
    # This should be investigated

def test_chemicaldrift_volatilization():
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
    [ "chemical:transformations:LogKOW"         , 6                     ],
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
    [ "chemical:transformations:volatilization" , True                 ],
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
        print("mass: ", sum(o[j].result.mass.values ))
        print("volatilized : ", sum(o[j].result.mass_volatilized.values ))
        print("total: ", sum(o[j].result.mass.values) + sum(o[j].result.mass_volatilized.values ))
        print("dissolved: ", sum(o[j].result.specie.values==0 ))        

    #assert 0
    
    # Check t12 500
    assert abs(sum(o[0].result.mass.values)[-1] - 59328.223) < 0.01
    assert abs(sum(o[0].result.mass_volatilized.values)[-1] - 140672.05) < 0.01
    
    # Check mass conservation residual + degraded = total
    assert abs(sum(o[0].result.mass.values)[-1] + sum(o[0].result.mass_volatilized.values)[-1] - 200000) < 0.5
