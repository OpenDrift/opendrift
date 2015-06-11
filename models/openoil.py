import os
import numpy as np

from openDriftSimulation import OpenDriftSimulation
from elements import LagrangianArray

# Defining the oil element properties
class Oil(LagrangianArray):
    """Extending LagrangianArray with variables relevant for oil particles."""
    
    variables = LagrangianArray.add_variables([
        ('mass_oil', {'dtype': np.float32,
                     'unit': 'kg',
                     'default': 1}),
        ('viscosity', {'dtype': np.float32,
                       'unit': 'mm2/s (centiStokes)',
                       'default': 100}),
        ('density', {'dtype': np.float32,
                     'unit': 'kg/cm^3',
                     'default': .8}),
        ('age_seconds', {'dtype': np.float32,
                         'unit': 's',
                         'default': 0}),
        ('age_exposure_seconds', {'dtype': np.float32,
                                  'unit': 's',
                                  'default': 0}),
        ('age_emulsion_seconds', {'dtype': np.float32,
                                  'unit': 's',
                                  'default': 0}),
        ('mass_emulsion', {'dtype': np.float32,
                          'unit': 'kg',
                          'default': 0}),
        ('fraction_evaporated', {'dtype': np.float32,
                          'unit': 'fraction',
                          'default': 0}),
        ('water_content', {'dtype': np.float32,
                          'unit': 'fraction',
                          'default': 0})])



class OpenOil(OpenDriftSimulation):
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_to_direction',
                          'x_wind', 'y_wind', 'land_binary_mask']

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_to_direction': np.nan,
                       'x_wind': 0, 'y_wind': 0}

    status_categories = {'initial': 'green',
                         'active': 'blue',
                         'dispersed': 'darkviolet',
                         'evaporated': 'yellow',
                         'stranded': 'red'}

    wind_factor = 0.02

    def __init__(self, *args, **kwargs):

        # Read oil properties from file
        d = os.path.dirname(os.path.realpath(__file__))
        oilprop = np.loadtxt(d + '/prop.dat')
        self.model.fref = oilprop[1:, 1]*.01  # Evaporation as fraction
        self.model.tref = oilprop[1:, 0]*3600  # Time in seconds
        self.model.wmax = oilprop[1:, 3]  # 
        self.model.reference_wind = oilprop[0,1]
        self.model.reference_thickness = oilprop[0,0]

        # Calling general constructor of parent class
        super(OpenOil, self).__init__(*args, **kwargs)


    def update(self):
        """Update positions and properties of oil particles."""

        self.elements.age_seconds += self.elements.age_seconds \
                                        + self.time_step.total_seconds()

        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)

        ################
        ## Evaporation
        ################
        # NB: temorarily assuming all particles are at ocean surface!

        Urel = windspeed/self.model.reference_wind  # Relative wind
        h = 2  # Film thickness in mm, harcoded for now
        self.elements.age_exposure_seconds += \
            (self.model.reference_thickness/h)*Urel* \
             self.time_step.total_seconds()

        self.elements.fraction_evaporated = np.interp(
            self.elements.age_exposure_seconds,
            self.model.tref, self.model.fref)


        ##################
        # Emulsification
        ##################

        # Apparent emulsion age of particles
        self.elements.age_emulsion_seconds += \
            Urel*self.time_step.total_seconds()

        self.elements.water_content = np.interp(
            self.elements.age_emulsion_seconds,
            self.model.tref, self.model.wmax)

        ###############
        # Dispersion
        ###############
        self.elements.depth = \
            self.environment.sea_surface_wave_significant_height
        whitecap_coverage = (2.54E-4)*np.power(windspeed, 3.58)  # In percent
        wave_period = 3.85*np.sqrt(
            self.environment.sea_surface_wave_significant_height)  # In seconds
        time_between_breaking_events = 3.85/whitecap_coverage  # TBC
        rho_w = 1025  # kg/m3
        wave_period[wave_period==0] = 5  # NB - temporal workaround when no waves
        dissipation_energy = 0.0034*rho_w*9.81*(wave_period**2)
        dsize_coeff = 2100  # Wettre et al., Appendix A
        c_oil = dsize_coeff*np.power(self.elements.viscosity, -0.4)
        p = np.random.rand(len(self.elements.lon), )  # Random numbers between 0 and 1
        oil_per_unit_surface = 1
        droplet_size = np.power((p*oil_per_unit_surface)/(c_oil*
                                np.power(dissipation_energy, 0.57)),
                                1/1.17)
        #self.deactivate_elements(droplet_size<1E-7, reason='dispersed')

        # Deactivate elements on land
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

        # Deactivate evaporated elements
        self.deactivate_elements(self.elements.fraction_evaporated > 0.7,
                                 reason='evaporated')

        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Wind drag
        self.update_positions(self.environment.x_wind*self.wind_factor,
                              self.environment.y_wind*self.wind_factor)


    def seed_from_gml(self, gmlfile, num_elements=1000):
        """Read oil slick contours from GML file, and seed particles within."""

        # Specific imports
        import datetime
        import matplotlib.nxutils as nx
        from xml.etree import ElementTree
        from matplotlib.patches import Polygon
        from mpl_toolkits.basemap import pyproj

        namespaces = {'od': 'http://cweb.ksat.no/cweb/schema/geoweb/oil',
                      'gml': 'http://www.opengis.net/gml'}
        slicks = []

        with open(gmlfile, 'rt') as e:
            tree = ElementTree.parse(e)
  
        pos1 = 'od:oilDetectionMember/od:oilDetection/od:oilSpill/gml:Polygon'
        pos2 = 'gml:exterior/gml:LinearRing/gml:posList'

        # This retrieves some other types of patches, found in some files only
        # Should be combines with the above, to get all patches
        #pos1 = 'od:oilDetectionMember/od:oilDetection/od:oilSpill/gml:Surface/gml:polygonPatches'
        #pos2 = 'gml:PolygonPatch/gml:exterior/gml:LinearRing/gml:posList'


        # Find detection time
        time_pos = 'od:oilDetectionMember/od:oilDetection/od:detectionTime'
        self.start_time = datetime.datetime.strptime(
                        tree.find(time_pos, namespaces).text,
                                  '%Y-%m-%dT%H:%M:%S.%fZ')

        for patch in tree.findall(pos1, namespaces):
            pos = patch.find(pos2, namespaces).text
            c = np.array(pos.split()).astype(np.float)
            lon = c[0::2]
            lat = c[1::2]
            slicks.append(Polygon(zip(lon, lat)))

        # Find boundary and area of all patches
        lons = np.array([])
        lats = lons.copy()
        for slick in slicks:
            ext = slick.get_extents()
            lons = np.append(lons, [ext.xmin, ext.xmax])
            lats = np.append(lats, [ext.ymin, ext.ymax])
            # Make a stereographic projection centred on the polygon
        lonmin = lons.min()
        lonmax = lons.max()
        latmin = lats.min()
        latmax = lats.max()
 
        # Place n points within the polygons
        proj = pyproj.Proj('+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f' %
                           (latmin, latmax, (latmin+latmax)/2, (lonmin+lonmax)/2))
        slickarea = np.array([])
        for slick in slicks:
            lonlat = slick.get_xy()
            lon = lonlat[:,0]
            lat = lonlat[:,1]
            x,y = proj(lon, lat)

            area_of_polygon = 0.0
            for i in xrange(-1, len(x)-1):
                area_of_polygon += x[i] * (y[i+1] - y[i-1])
            area_of_polygon = abs(area_of_polygon) / 2.0
            slickarea = np.append(slickarea, area_of_polygon) # in m2

        # Make points
        deltax = np.sqrt(np.sum(slickarea)/num_elements)

        lonpoints = np.array([])
        latpoints = np.array([])
        for i, slick in enumerate(slicks):
            lonlat = slick.get_xy()
            lon = lonlat[:,0]
            lat = lonlat[:,1]
            x,y = proj(lon, lat)
            xvec = np.arange(x.min(), x.max(), deltax)
            yvec = np.arange(y.min(), y.max(), deltax)
            x, y = np.meshgrid(xvec, yvec)
            lon, lat = proj(x, y, inverse=True)
            lon = lon.ravel()
            lat = lat.ravel()
            points = np.c_[lon, lat]
            ind = nx.points_inside_poly(points, slick.xy)
            lonpoints = np.append(lonpoints, lon[ind])
            latpoints = np.append(latpoints, lat[ind])

        # Finally seed at found positions
        kwargs = {}
        kwargs['lon'] = lonpoints
        kwargs['lat'] = latpoints
        kwargs['ID'] = np.arange(len(lonpoints)) + 1
        kwargs['mass_oil'] = 1
        self.elements = self.ElementType(**kwargs)
