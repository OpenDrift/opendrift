#!/usr/bin/env python
"""
Plotting map
===============
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import shapely

print(cartopy.__version__)
print(shapely.__version__)

fig = plt.figure()
sp = ccrs.Stereographic(central_longitude=0, central_latitude=60)
ax = fig.add_subplot(1, 1, 1, projection=sp)

corners = [-10, 10, 55, 65]
f = cfeature.GSHHSFeature(scale='h', levels=[1])
ax.add_geometries(
    f.intersecting_geometries(corners),
    ccrs.PlateCarree(),
    facecolor=cfeature.COLORS['land'],
    edgecolor='black')

ax.set_extent(corners, crs=ccrs.PlateCarree())
gl = ax.gridlines(ccrs.PlateCarree())

plt.show()
