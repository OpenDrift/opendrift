import cartopy.crs as ccrs
import pyproj
import utm


def simulation():
    """
    Simulation Spatial Reference System (SRS) or Coordinate Reference System (CRS).

    Always just 'lonlat' ellipsoid 'WGS84'
    """
    return pyproj.Proj('+proj=lonlat')


def find_utm_proj(lon, lat):
    """
    Find a suitable UTM projection (zone) for lon and lat.

    .. warning::

      UTM is only defined between 80S and 84N. Should use UPS for those regions.

    Returns:

      pyproj.Proj in `utm` projection.
    """
    _, _, zone_no, _ = utm.from_latlon(lat, lon)
    band = 'south' if lat < 0 else 'north'

    return pyproj.Proj(
        '+proj=utm +zone={zone:d} +{band} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
        .format(zone=zone_no, band=band))

def find_utm_ccrs(lon, lat):
  """
  Return Cartopy UTM CRS for lon and lat.
  """
  _, _, zone_no, _ = utm.from_latlon(lat, lon)
  return ccrs.UTM(zone_no, True if lat < 0 else False)

