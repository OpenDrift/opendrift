#!/usr/bin/env python
"""
Seeding from shapefile
==================================
"""

from datetime import datetime
import geopandas as gpd
from opendrift.models.oceandrift import OceanDrift


o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Seed particles within all polygons from shapefile
o.seed_from_shapefile(o.test_data_folder() +
                      'shapefile_spawning_areas/Torsk.shp',
                      number=1000, time=datetime.utcnow())


#%%
# Alternative by using GeoPandas and selecting a subset of two polygons
gdf = gpd.read_file(o.test_data_folder() + 'shapefile_spawning_areas/Torsk.shp')
gdf = gdf.iloc[[1, 5]]
o.seed_from_geopandas(gdf, number=2000, time=datetime.utcnow())

o.plot(fast=True)
