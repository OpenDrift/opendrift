import numpy as np
import pyproj
from matplotlib.patches import Polygon
from matplotlib.path import Path

class Seeder:
    def seed_elements(self, lon, lat, time, radius=0, number=None,
                      radius_type='gaussian', **kwargs):
        """Seed elements with given position(s), time and properties.

        Arguments:
            lon: scalar or array
                central longitude(s).
            lat: scalar or array
                central latitude(s).
            radius: scalar or array
                radius in meters around each lon-lat pair,
                within which particles will be randomly seeded.
            number: integer, total number of particles to be seeded
                If number is None, the number of elements is the
                length of lon/lat or time if these are arrays. Otherwise
                the number of elements are obtained from the config-default.
            time: datenum or list
                The time at which particles are seeded/released.
                If time is a list with two elements, elements are seeded
                continously from start/first to end/last time.
                If time is a list with more than two elements, the number of elements
                is equal to len(time) and are seeded as a time series.
            radius_type: string
                If 'gaussian' (default), the radius is the standard deviation in
                x-y-directions. If 'uniform', elements are spread evenly and
                always inside a circle with the given radius.
            kwargs:
                keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """

        if 'cone' in kwargs and kwargs['cone'] is True:
            self.logger.warning('Keyword *cone* for seed_elements is deprecated, use seed_cone() instead.')
            del kwargs['cone']
            self.seed_cone(lon, lat, time, radius, number, **kwargs)
            return

        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        radius = np.atleast_1d(radius).ravel()
        time = np.atleast_1d(time)

        if len(lon) != len(lat):
            raise ValueError('Lon and lat must have same lengths')

        if len(lon) > 1:
            if number is not None and number != len(lon):
                raise ValueError('Lon and lat have length %s, but number is %s' % (len(lon), number))
            number = len(lon)
        else:
            if number is None:
                if len(time) > 2:
                    number = len(time)  # Interpreting as time series
                else:
                    number = self.get_config('seed:number')
            lon = lon*np.ones(number)
            lat = lat*np.ones(number)

        if len(time) != number and len(time) > 1:
            if len(time) == 2:  # start -> end
                td = (time[1]-time[0])/(number-1)  # timestep between points
                time = [time[0] + i*td for i in range(number)]
            else:
                raise ValueError('Time array has length %s, must be 1, 2 or %s' % (len(time), number))

        # Add radius / perturbation
        if radius.max() > 0:
            geod = pyproj.Geod(ellps='WGS84')
            ones = np.ones(np.sum(number))
            if radius_type == 'gaussian':
                x = np.random.randn(np.sum(number))*radius
                y = np.random.randn(np.sum(number))*radius
                az = np.degrees(np.arctan2(x, y))
                dist = np.sqrt(x*x+y*y)
            elif radius_type == 'uniform':
                az = np.random.randn(np.sum(number))*360
                dist = np.sqrt(np.random.uniform(0, 1, np.sum(number)))*radius
            lon, lat, az = geod.fwd(lon, lat, az, dist, radians=False)

        # If z is 'seafloor'
        if not 'z' in kwargs or kwargs['z'] is None:
            if 'seed:seafloor' in self._config:
                if self.get_config('seed:seafloor') is True:
                    kwargs['z'] = 'seafloor'
                    self.logger.debug('Seafloor is selected, neglecting z')
        if 'z' in kwargs and isinstance(kwargs['z'], str) \
                and kwargs['z'][0:8] == 'seafloor':
            # We need to fetch seafloor depth from reader
            seafloor_constant = self.get_config('environment:constant:sea_floor_depth_below_sea_level')
            seafloor_fallback = self.get_config('environment:fallback:sea_floor_depth_below_sea_level')
            if seafloor_constant is not None:
                env = {'sea_floor_depth_below_sea_level': np.array(seafloor_constant)}
            elif ('sea_floor_depth_below_sea_level' in self.priority_list
                        ) or len(self._lazy_readers()):
                if not hasattr(self, 'time'):
                    self.time = time[0]
                env, env_profiles, missing = \
                    self.get_environment(['sea_floor_depth_below_sea_level'],
                                         time=time[0], lon=lon, lat=lat,
                                         z=0*lon, profiles=None)
            elif seafloor_fallback is not None:
                env = {'sea_floor_depth_below_sea_level': np.array(seafloor_fallback)}
            else:
                raise ValueError('A reader providing the variable '
                                 'sea_floor_depth_below_sea_level must be '
                                 'added before seeding elements at seafloor.')
            # Add M meters if given as 'seafloor+M'
            if len(kwargs['z']) > 8 and kwargs['z'][8] == '+':
                meters_above_seafloor = np.float(kwargs['z'][9::])
                self.logger.info('Seeding elements %f meters above seafloor'
                             % meters_above_seafloor)
            else:
                meters_above_seafloor = 0
            kwargs['z'] = \
                -env['sea_floor_depth_below_sea_level'].astype('float32') + meters_above_seafloor

        # Creating and scheduling elements
        elements = self.ElementType(lon=lon, lat=lat, **kwargs)
        time_array = np.array(time)
        self.schedule_elements(elements, time)

    def seed_cone(self, lon, lat, time, radius=0, number=None, **kwargs):
        """Seed elements along a transect/cone between two points/times

        Arguments:
            lon: scalar or list with 2 elements [lon0, lon1]
            lat: scalar or list with 2 elements [lat0, lat]
            time: datetime or list with 2 elements [t0, t1]
            radius: scalar or list with 2 elements [r0, r1] Unit: meters
            number: int
                The number of elements. If this is None, the number of
                elements is taken from configuration.

        Elements are seeded along a transect from
            (lon0, lat0) with uncertainty radius r0 at time t0, towards
            (lon1, lat1) with uncertainty radius r1 at time t1.
            If r0 != r1, the unceetainty radius is linearly changed along
            the transect, thus outlining a "cone".
        """

        if number is None:
            number = self.get_config('seed:number')
        if number == 1:
            raise ValueError('For a cone, the number of elements must be at least 2 or more, given is 1')

        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        radius = np.atleast_1d(radius).ravel()
        if len(lon) != len(lat):
            raise ValueError('Lon and lat must have same length (1 or 2)')
        elif len(lon) > 2:
            raise ValueError('Lon and lat must have length 1 or 2, given length is %s' % (len(lon)))
        elif len(lon) == 1:
            lon = lon*np.ones(number)
            lat = lat*np.ones(number)
        elif len(lon) == 2:  # Segment from lon0,lat1 to lon1,lat2
            geod = pyproj.Geod(ellps='WGS84')
            # Note that npts places points in-between start and end, and does not include these
            conelonlats = geod.npts(lon[0], lat[0], lon[1], lat[1],
                                    number, radians=False)
            lon, lat = zip(*conelonlats)

        if len(radius) > 2:
            raise ValueError('Seed radius must have length 1 or 2')
        elif len(radius) == 2:  # Linear increase from r0 to r1
            radius = np.linspace(radius[0], radius[1], number)

        # Forwarding calculated cone points/radii to seed_elements
        self.seed_elements(lon=lon, lat=lat, time=time, radius=radius, number=number, **kwargs)

    def seed_from_geojson(self, gjson):
        """Under development"""
        import geojson
        try:
            gj = geojson.loads(gjson)
        except:
            raise ValueError('Could not load GeoJSON string: %s' % gjson)
        if not gj.is_valid:
            raise ValueError('GeoJSON string is not valid: %s' % gj.errors())
        # Assuming temporally that g is a Feature, and not a FeatureCollection
        properties = gj['properties']
        if 'time' not in properties:
            raise ValueError('Property "time" is not available')
        kwargs = {}
        for prop in properties:
            if prop == 'time':
                t = properties['time']
                if isinstance(t, list):
                    time = [datetime.fromisoformat(t[0].replace("Z", "+00:00")),
                            datetime.fromisoformat(t[1].replace("Z", "+00:00"))]
                else:
                    time = datetime.fromisoformat(t.replace("Z", "+00:00"))
            else:
                kwargs[prop] = properties[prop]

        geometry = gj['geometry']

        if geometry['type'] == 'Polygon':
            coords = list(geojson.utils.coords(gj))
            lon, lat = zip(*[(c[0], c[1]) for c in coords])
            self.seed_within_polygon(lons=lon, lats=lat, time=time, **kwargs)
        elif geometry['type'] == 'LineString':
            coords = list(geojson.utils.coords(gj))
            lon, lat = zip(*[(c[0], c[1]) for c in coords])
            self.seed_cone(lon=lon, lat=lat, time=time, **kwargs)
        elif geometry['type'] == 'Point':
            coords = list(geojson.utils.coords(gj))
            lon, lat = zip(*[(c[0], c[1]) for c in coords])
            self.seed_elements(lon=lon, lat=lat, time=time, **kwargs)
        else:
            raise ValueError('Not yet implemented')

    def seed_repeated_segment(self, lons, lats,
                              start_time, end_time, time_interval=None,
                              number_per_segment=None,
                              total_number=None, **kwargs):
        """Seed elements repeatedly in time along a segment.

        The segment goes from lon[0],lat[0] to lon[1],lat[1].
        The number of elements should be proved as either:
        1) number_per_segment, in which case total number of elements
           is number_per_segment * len(times), or
        2) total_number, in which case the number of elements
           per segment is: total_number / len(times).
           Any extra elements are duplicated along at the first segment.

        """

        numtimes = int((end_time-start_time).total_seconds()/
                        time_interval.total_seconds() + 1)
        times = [start_time+i*time_interval for i in range(numtimes)]

        geod = pyproj.Geod(ellps='WGS84')
        if number_per_segment is None:
            number_per_segment = np.int(np.floor(total_number/numtimes))

        s_lonlats= geod.npts(lons[0], lats[0], lons[1], lats[1],
                             number_per_segment, radians=False)
        slon, slat = list(zip(*s_lonlats))
        slon = np.atleast_1d(slon)
        slat = np.atleast_1d(slat)

        lon, time = np.meshgrid(slon, times)
        lat, time = np.meshgrid(slat, times)
        lon = lon.ravel()
        lat = lat.ravel()
        time = time.ravel()

        if total_number is not None:
            additional_elements = total_number - len(lon.ravel())
            print('Repeating the %d last points, to obtain %d elements' %
                  (additional_elements, total_number))
            lon = np.concatenate((lon, lon[-additional_elements::]))
            lat = np.concatenate((lat, lat[-additional_elements::]))
            time = np.concatenate((time, time[-additional_elements::]))

        self.seed_elements(lon=lon, lat=lat, time=time, **kwargs)

    def seed_within_polygon(self, lons, lats, number=None, **kwargs):
        """Seed a number of elements within given polygon.

        Arguments:
            lon: array of longitudes

            lat: array of latitudes

            number: int, number of elements to be seeded

            kwargs: keyword arguments containing properties/attributes and
            values corresponding to the actual particle type (ElementType).
            These are forwarded to method seed_elements(). All properties
            for which there are no default value must be specified.

        """
        if number == 0:
            return

        if number is None:
            number = self.get_config('seed:number')

        lons = np.asarray(lons)
        lats = np.asarray(lats)
        if len(lons) < 3:
            self.logger.info('At least three points needed to make a polygon')
            return
        if len(lons) != len(lats):
            raise ValueError('lon and lat arrays must have same length.')
        poly = Polygon(list(zip(lons, lats)), closed=True)
        # Place N points within the polygons
        proj = pyproj.Proj('+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f '
                           '+lon_0=%f +R=6370997.0 +units=m +ellps=WGS84'
                           % (lats.min(), lats.max(),
                              (lats.min()+lats.max())/2,
                              (lons.min()+lons.max())/2))
        lonlat = poly.get_xy()
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = proj(lon, lat)
        area = 0.0
        for i in range(-1, len(x)-1):
            area += x[i] * (y[i+1] - y[i-1])
        area = abs(area) / 2

        # Make points, evenly distributed
        deltax = np.sqrt(area/number)
        lonpoints = np.array([])
        latpoints = np.array([])
        lonlat = poly.get_xy()
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = proj(lon, lat)
        xvec = np.linspace(x.min() + deltax/2, x.max() - deltax/2,
                           int((x.max()-x.min())/deltax))
        yvec = np.linspace(y.min() + deltax/2, y.max() - deltax/2,
                           int((y.max()-y.min())/deltax))
        x, y = np.meshgrid(xvec, yvec)
        lon, lat = proj(x, y, inverse=True)
        lon = lon.ravel()
        lat = lat.ravel()
        points = np.c_[lon, lat]
        ind = Path(poly.xy).contains_points(points)
        if not any(ind):  # No elements are inside, we seed on border
            lonpoints = np.append(lonpoints, lons[0:number])
            latpoints = np.append(latpoints, lats[0:number])
        else:
            lonpoints = np.append(lonpoints, lon[ind])
            latpoints = np.append(latpoints, lat[ind])
        if len(ind) == 0:
            self.logger.info('Small or irregular polygon, using center point.')
            lonpoints = np.atleast_1d(np.mean(lons))
            latpoints = np.atleast_1d(np.mean(lats))
        # Truncate if too many
        # NB: should also repeat some points, if too few
        lonpoints = lonpoints[0:number]
        latpoints = latpoints[0:number]
        if len(lonpoints) < number:
            # If number of positions is smaller than requested,
            # we duplicate the first ones
            missing = number - len(lonpoints)
            lonpoints = np.append(lonpoints, lonpoints[0:missing])
            latpoints = np.append(latpoints, latpoints[0:missing])

        # Finally seed at calculated positions
        self.seed_elements(lonpoints, latpoints, number=number,
                           **kwargs)

    def seed_from_wkt(self, wkt, number=None, **kwargs):
        """Seeds elements within (multi)polygons from WKT"""

        try:
            from osgeo import ogr, osr
        except Exception as e:
            self.logger.warning(e)
            raise ValueError('OGR library is needed to parse WKT')

        if number is None:
            number = self.get_config('seed:number')

        geom = ogr.CreateGeometryFromWkt(wkt)
        total_area = 0
        for i in range(0, geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            total_area += g.GetArea()

        self.logger.info('Total area of all polygons: %s m2' % total_area)
        num_seeded = 0
        for i in range(0, geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            num_elements = np.int(number*g.GetArea()/total_area)
            if i == geom.GetGeometryCount()-1:
                # For the last feature we seed the remaining number,
                # avoiding difference due to rounding:
                num_elements = number - num_seeded
            self.logger.info('\tSeeding %s elements within polygon number %s' %
                         (num_elements, str(i)))
            try:
                g.Transform(coordTrans)
            except:
                pass
            b = g.GetBoundary()
            if b is not None:
                points = b.GetPoints()
                lons = [p[0] for p in points]
                lats = [p[1] for p in points]
            else:
                # Alternative if OGR is not built with GEOS support
                r = g.GetGeometryRef(0)
                lons = [r.GetX(j) for j in range(r.GetPointCount())]
                lats = [r.GetY(j) for j in range(r.GetPointCount())]

            self.seed_within_polygon(lons, lats, num_elements, **kwargs)
            num_seeded += num_elements

    def seed_from_shapefile(self, shapefile, number,
                            layername=None, featurenum=None, **kwargs):
        """Seeds elements within contours read from a shapefile"""

        try:
            from osgeo import ogr, osr
        except Exception as e:
            self.logger.warning(e)
            raise ValueError('OGR library is needed to read shapefiles.')

        if 'timeformat' in kwargs:
            # Recondstructing time from filename, where 'timeformat'
            # is forwarded to datetime.strptime()
            kwargs['time'] = datetime.strptime(os.path.basename(shapefile),
                                               kwargs['timeformat'])
            del kwargs['timeformat']

        num_seeded_before = self.num_elements_scheduled()

        targetSRS = osr.SpatialReference()
        targetSRS.ImportFromEPSG(4326)
        try:
            s = ogr.Open(shapefile)
        except:
            s = shapefile

        for layer in s:
            if layername is not None and layer.GetName() != layername:
                self.logger.info('Skipping layer: ' + layer.GetName())
                continue
            else:
                self.logger.info('Seeding for layer: %s (%s features)' %
                             (layer.GetDescription(), layer.GetFeatureCount()))

            coordTrans = osr.CoordinateTransformation(layer.GetSpatialRef(),
                                                      targetSRS)

            if featurenum is None:
                featurenum = range(1, layer.GetFeatureCount() + 1)
            else:
                featurenum = np.atleast_1d(featurenum)
            if max(featurenum) > layer.GetFeatureCount():
                raise ValueError('Only %s features in layer.' %
                                 layer.GetFeatureCount())

            # Loop first through all features to determine total area
            layer.ResetReading()
            area_srs = osr.SpatialReference()
            area_srs.ImportFromEPSG(3857)
            areaTransform = osr.CoordinateTransformation(layer.GetSpatialRef(), area_srs)

            areas = np.zeros(len(featurenum))
            for i, f in enumerate(featurenum):
                feature = layer.GetFeature(f - 1)  # Note 1-indexing, not 0
                if feature is not None:
                    gom = feature.GetGeometryRef().Clone()
                    gom.Transform(areaTransform)
                    areas[i] = gom.GetArea()

            total_area = np.sum(areas)
            layer.ResetReading()  # Rewind to first layer
            self.logger.info('Total area of all polygons: %s m2' % total_area)
            # Find number of points per polygon
            numbers = np.round(number*areas/total_area).astype(int)
            numbers[numbers.argmax()] += np.int(number-sum(numbers))

            for i, f in enumerate(featurenum):
                feature = layer.GetFeature(f - 1)
                if feature is None:
                    continue
                num_elements = numbers[i]
                geom = feature.GetGeometryRef()
                self.logger.info('\tSeeding %s elements within polygon number %s' %
                             (num_elements, featurenum[i]))
                try:
                    geom.Transform(coordTrans)
                except Exception as e:
                    logging.warning('Could not transform coordinates:')
                    logging.warning(e)
                    pass
                #b = geom.GetBoundary()
                #if b is not None:
                #    points = b.GetPoints()
                #    lons = [p[0] for p in points]
                #    lats = [p[1] for p in points]
                #else:
                # Alternative if OGR is not built with GEOS support
                r = geom.GetGeometryRef(0)
                lons = [r.GetY(j) for j in range(r.GetPointCount())]
                lats = [r.GetX(j) for j in range(r.GetPointCount())]

                self.seed_within_polygon(lons, lats, num_elements, **kwargs)

    def seed_from_ladim(self, ladimfile, roms):
        """Seed elements from ladim \*.rls text file: [time, x, y, z, name]"""

        data = np.loadtxt(ladimfile,
            dtype={'names': ('time', 'x', 'y', 'z'),
                   'formats': ('S20', 'f4', 'f4', 'f4')},
            usecols=(0,1,2,3))

        time = [datetime.strptime(t, "%Y-%m-%dT%H")
                for t in data['time']]
        time = np.array(time)

        lon, lat = roms.xy2lonlat(data['x'], data['y'])
        z = -data['z']

        self.logger.info('Seeding %i elements from %s:' % (len(lon), ladimfile))
        self.logger.info('    Lons: %f to %f' % (lon.min(), lon.max()))
        self.logger.info('    Lats: %f to %f' % (lat.min(), lat.max()))
        self.logger.info('    Depths: %f to %f' % (z.min(), z.max()))
        self.logger.info('    Time: %s to %s' % (time.min(), time.max()))
        elements = self.ElementType(lon=lon, lat=lat, z=-z)

        self.schedule_elements(elements, time)

