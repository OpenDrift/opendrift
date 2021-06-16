# Som valid (but extreme) ranges for checking that values are reasonable
standard_names = {
    'x_wind': {'valid_min': -50, 'valid_max': 50, 'units': 'm/s',
               'long_name': 'Component of wind along x-direction '
                            '(eastwards if projection is lonlat/Mercator)'},
    'y_wind': {'valid_min': -50, 'valid_max': 50, 'units': 'm/s',
        'long_name': 'Component of wind along y-direction (northwards if '
                      'projection is lonlat/Mercator)'},

    'x_sea_water_velocity': {'valid_min': -15, 'valid_max': 15,
        'units': 'm/s',
        'long_name': 'Component of ocean current along x-direction '
                     ' (eastwards if projection is lonlat/Mercator)'},
    'y_sea_water_velocity': {'valid_min': -15, 'valid_max': 15, 'units': 'm/s',

        'long_name': 'Component of ocean current along y-direction '
                      '(northwards if projection is lonlat/Mercator)'},
    'land_binary_mask': {'valid_min': 0, 'valid_max': 1,
                         'long_name': '1 is land, 0 is sea'},
    'sea_floor_depth_below_sea_level': {'valid_min': 0, 'valid_max': 12000,
                         'long_name': 'Depth of seafloor'},
    'ocean_vertical_diffusivity': {'valid_min': 0, 'valid_max': 1}}

# Identify x-y vector components/pairs for rotation (NB: not east-west pairs!)
vector_pairs_xy = [
    ['x_wind', 'y_wind'],
    ['sea_ice_x_velocity', 'sea_ice_y_velocity'],
    ['x_sea_water_velocity', 'y_sea_water_velocity'],
    ['sea_surface_wave_stokes_drift_x_velocity',
     'sea_surface_wave_stokes_drift_y_velocity']
    ]
