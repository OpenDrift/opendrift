#!/usr/bin/env python

### Sea Surface Slope Computation
import numpy as np
from netCDF4 import Dataset


### Get ssh from TP5
TP5path = '/cluster/work/users/achoth/MonthlyData/TP5/Hyc2proj/'
TP5file = TP5path + 'archm_20190112_00.nc'
TP5data = Dataset(TP5file, mode='r')
ssh = TP5data.variables['zos'][:].filled(np.nan)
plon = TP5data.variables['longitude'][:].filled(np.nan)
plat = TP5data.variables['latitude'][:].filled(np.nan)
ssh=np.squeeze(ssh)
I,J=plon.shape


def find_nearest_grid_point(plon, plat, lon_ib, lat_ib):
    dist = np.sqrt((plon - lon_ib) ** 2 + (plat - lat_ib) ** 2)
    II, JJ = np.unravel_index(dist.argmin(), dist.shape)
    return II, JJ


### Find nearest grid point to the Iceberg's position
lon_ib= -57
lat_ib= 69
### Verification
II,JJ= find_nearest_grid_point(plon, plat, lon_ib, lat_ib)
lon=plon[II,JJ]
lat=plat[II,JJ]
print(lon, lat)


def sea_surface_slope_force(mass, ssh, lat):
    # Constants
    g = 9.81  # Acceleration due to gravity in m/s²

    dx_lon_deg = 0.01  # Longitude spacing in degrees (fixed for testing)
    dy_lat_deg = 0.01  # Latitude spacing in degrees (fixed for testing)

    # Convert to meters
    dy_lat_m = dy_lat_deg * 111320  # Approx. meters for latitude
    dx_lon_m = dx_lon_deg * (111320 * np.cos(np.radians(lat)))

    # Get neighboring SSH values
    ssh_north = ssh[II, JJ + 1]
    ssh_south = ssh[II, JJ - 1] 
    ssh_east = ssh[II + 1, JJ]
    ssh_west = ssh[II - 1, JJ]
    print('ssh_north =', ssh_north, '\nssh_south =', ssh_south,
          '\nssh_east =', ssh_east, '\nssh_west =', ssh_west)

    # Calculate the slope in the x and y directions
    slope_x = (ssh_east - ssh_west) / (2 * dx_lon_m)
    slope_y = (ssh_north - ssh_south) / (2 * dy_lat_m)
    print('slope_x =', slope_x, '\nslope_y =', slope_y)

    # Calculate the forces in x and y directions
    F_slope_x = mass * g * slope_x
    F_slope_y = mass * g * slope_y

    print('Slope in the x direction =', F_slope_x)
    print('Slope in the y direction =', F_slope_y)
    
    #overall_slope = np.sqrt(F_slope_x**2 + F_slope_y**2)
    #print('Overall slope:', overall_slope)

    return np.array([F_slope_x, F_slope_y])


### Test
#mass=1000
#F_slope_x, F_slope_y = sea_surface_slope_force(mass, ssh, lat)