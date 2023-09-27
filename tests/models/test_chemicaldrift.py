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
        fOC, organic carbon fraction, in suspended particles.
        High SPM concentration, 80 g/m3."""

    days=300
    ntraj=200

    # If one parameter is given as a list, different instances of
    # ChemicalDrift will be created
    configs = [
    [ "seed:LMM_fraction"                       , 1                     ],
    [ "seed:particle_fraction"                  , 0                     ],
    [ "chemical:species:Sediment_reversible"    , True                  ],
    [ "chemical:dynamic_partitioning"           , True                  ],
    [ "chemical:transfer_setup"                 , 'organics'            ],
    [ "chemical:transformations:fOC_SPM"        , [0.01, 0.05 ,0.1]     ],
    [ "chemical:transformations:LogKOW"         , 6.618                 ],
    [ "chemical:transformations:DeltaH_KOC_Sed" , 0                     ],
    [ "chemical:transformations:Setchenow"      , 0                     ],
    [ "chemical:particle_diameter"              , 0                     ],
    [ "chemical:particle_diameter_uncertainty"  , 0                     ],
    [ "drift:vertical_mixing"                   , False                 ],
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
                o.append(ChemicalDrift(loglevel=100))
                o[-1].set_config(configs[i][0],configs[i][1][j])

    if o==[]:
        j=0
        o.append(ChemicalDrift(loglevel=100))

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

        o[j].run(duration=timedelta(hours=24*days),time_step=60*60*168,time_step_output=60*60*168)

    assert sum(o[0].elements.specie==2)/o[0].num_elements_total()*100 == 35
    assert sum(o[1].elements.specie==2)/o[1].num_elements_total()*100 == 72.5
    assert sum(o[2].elements.specie==2)/o[2].num_elements_total()*100 == 84
