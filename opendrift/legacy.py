# This is a collection of methods that have been removed from various places in the OpenDrift codebase.


# See earlier version of reader_ROMS_native how the GLS parameters were obtained
def gls_tke(windstress, depth, sea_water_density,
            tke, generic_length_scale, gls_parameters=None):
    '''From LADIM model, based on ROMS files.'''

    g = 9.81
    f0 = 0.1  # mean wave frequency
    c_w = 4.0  # wave mixing parameter
    c_i = 0.2  # coefficient for the interior
    if gls_parameters is None:
        # GLS parameters from ROMS, k-omega closure (see ocean.in)
        p = 0.0
        m = 1.0
        n = 1.0
        cmu0 = 0.5477  # for KANTHA_CLAYSON stability function
    else:
        p = gls_parameters['gls_p']
        m = gls_parameters['gls_m']
        n = gls_parameters['gls_n']
        cmu0 = gls_parameters['gls_cmu0']

    phi = 100. * (windstress/sea_water_density)**(3./2.)

    # dissipation and turbulent length scale for interiour of mixed layer
    eps = cmu0**(3.+p/n)*tke**(3./2.+m/n)*generic_length_scale**(-1./n)
    l_i = c_i * tke**(3./2.) * eps**(-1.)

    # diffusivity for interior of mixed layer
    # c_i = sqrt(2.) * cmu0**3
    ki = c_i * (2.*tke)**0.5 * l_i

    # length scale and diffusivity of wave-enhanced layer
    l_w = np.sqrt(phi / (g*f0))
    kwave = c_w * (2*tke)**0.5 * l_w
    kmix = ki + kwave

    K, N = np.meshgrid(kmix, depths)

    return K

#@require_mode(mode=Mode.Ready)
def seed_from_shapefile_old(self,
                        shapefile,
                        number,
                        layername=None,
                        featurenum=None,
                        **kwargs):
    """Seeds elements within contours read from a shapefile

    Obsolete, as new method based on geopandas is simpler"""

    try:
        from osgeo import ogr, osr
    except Exception as e:
        logger.warning(e)
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
            logger.info('Skipping layer: ' + layer.GetName())
            continue
        else:
            logger.info('Seeding for layer: %s (%s features)' %
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
        areaTransform = osr.CoordinateTransformation(
            layer.GetSpatialRef(), area_srs)

        areas = np.zeros(len(featurenum))
        for i, f in enumerate(featurenum):
            feature = layer.GetFeature(f - 1)  # Note 1-indexing, not 0
            if feature is not None:
                gom = feature.GetGeometryRef().Clone()
                gom.Transform(areaTransform)
                areas[i] = gom.GetArea()

        total_area = np.sum(areas)
        layer.ResetReading()  # Rewind to first layer
        logger.info('Total area of all polygons: %s m2' % total_area)
        # Find number of points per polygon
        numbers = np.round(number * areas / total_area).astype(int)
        numbers[numbers.argmax()] += int(number - sum(numbers))

        for i, f in enumerate(featurenum):
            feature = layer.GetFeature(f - 1)
            if feature is None:
                continue
            num_elements = numbers[i]
            geom = feature.GetGeometryRef()
            logger.info(f'\tSeeding {num_elements} elements within polygon number {featurenum[i]} of area {areas[i]} m3')
            try:
                geom.Transform(coordTrans)
            except Exception as e:
                logger.warning('Could not transform coordinates:')
                logger.warning(e)
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

            self.seed_within_polygon(lons=lons,
                                        lats=lats,
                                        number=num_elements,
                                        **kwargs)