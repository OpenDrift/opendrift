import numpy as np
import pyproj
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Sketch of an analysis package for trajectory datasets, which is independent of OpenDrift

def set_up_map(td=None, buffer=.1, corners=None, landmask='auto', lscale='auto', fast=True,
               ocean_color='white', land_color=cfeature.COLORS['land'], figsize=11):

    if corners is None:
        lonmin = td.lon.min() - buffer
        lonmax = td.lon.max() + buffer
        latmin = td.lat.min() - buffer
        latmax = td.lat.max() + buffer
    else:
        lonmin = corners[0]
        lonmax = corners[1]
        latmin = corners[2]
        latmax = corners[3]

    crs = ccrs.Mercator()
    globe = crs.globe
    gcrs = ccrs.PlateCarree(globe=crs.globe)  # For coordinate transform
    meanlat = (latmin + latmax) / 2
    aspect_ratio = float(latmax - latmin) / (float(lonmax - lonmin))
    aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))

    plt.close()  # Close any existing windows, so that they dont show up with new
    if aspect_ratio > 1:
        fig = plt.figure(figsize=(figsize / aspect_ratio, figsize))
    else:
        fig = plt.figure(figsize=(figsize, figsize * aspect_ratio))

    ax = fig.add_subplot(111, projection=crs)
    ax.set_extent([lonmin, lonmax, latmin, latmax], crs=gcrs)

    gl = ax.gridlines(gcrs, draw_labels=True)
    gl.top_labels = None

    fig.canvas.draw()
    fig.set_tight_layout(True)

    ####################
    # Landmask
    ####################
    if landmask == 'auto':
        from opendrift.readers import reader_global_landmask
        reader_global_landmask.plot_land(ax, lonmin, latmin, lonmax, latmax,
                                         fast, ocean_color, land_color, lscale, globe)

    return ax, fig, gcrs

def plot(td, background=None, show=True, trajectory_kwargs={}, map_kwargs={}):

    # td: trajectory Xarray Dataset
    # background: DataArray with background field (not yet implemented)

    ax, fig, gcrs = set_up_map(td, **map_kwargs)

    if 'trajectory' in td.dims:
        num_trajectories = len(td.trajectory)
    else:
        num_trajectories = 1  # E.g. removed by averaging

    # Default trajectory options
    if 'alpha' not in trajectory_kwargs:
        # The more trajectories, the more transparent we make the lines
        min_alpha = 0.1
        max_trajectories = 5000.0
        alpha = min_alpha**(2 * (num_trajectories - 1) / (max_trajectories - 1))
        trajectory_kwargs['alpha'] = np.max((min_alpha, alpha))
    if 'color' not in trajectory_kwargs:
        trajectory_kwargs['color'] = 'gray'

    ####################
    # Plot trajectories
    ####################
    ax.plot(td.lon.T, td.lat.T, transform=gcrs, **trajectory_kwargs)

    if show is True:
        plt.show()
    else:
        return ax, fig, gcrs


def skillscore_liu_weissberg(lon_obs, lat_obs, lon_model, lat_model, tolerance_threshold=1):
    ''' calculate skill score from normalized cumulative seperation
    distance. Liu and Weisberg 2011.

    Returns the skill score bewteen 0. and 1.
    '''

    lon_obs = np.array(lon_obs)
    lat_obs = np.array(lat_obs)
    lon_model = np.array(lon_model)
    lat_model = np.array(lat_model)
    d = distance_between_trajectories(lon_obs, lat_obs, lon_model, lat_model)
    l = distance_along_trajectory(lon_obs, lat_obs)
    s = d.sum() / l.cumsum().sum()
    if tolerance_threshold==0:
        skillscore = 0
    else:
        skillscore = max(0, 1 - s/tolerance_threshold)

    return skillscore

def distance_between_trajectories(lon1, lat1, lon2, lat2):
    '''Calculate the distances [m] between two trajectories'''

    assert len(lon1) == len(lat1) == len(lat1) == len(lat2)
    geod = pyproj.Geod(ellps='WGS84')
    azimuth_forward, a2, distance = geod.inv(lon1, lat1, lon2, lat2)

    return distance

def distance_along_trajectory(lon, lat):
    '''Calculate the distances [m] between points along a trajectory'''

    geod = pyproj.Geod(ellps='WGS84')
    azimuth_forward, a2, distance = geod.inv(lon[1:], lat[1:], lon[0:-1], lat[0:-1])

    return distance

def skillscore_along_trajectory(lon_obs, lat_obs, time_obs,
                                lon_model, lat_model, time_model,
                                method='liu-weissberg', **kwargs):
