# Som valid (but extreme) ranges for checking that values are reasonable
standard_names = {
    'x_wind': {'valid_min': -50, 'valid_max': 50},
    'y_wind': {'valid_min': -50, 'valid_max': 50},
    'x_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'y_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'ocean_vertical_diffusivity': {'valid_min': 0, 'valid_max': 1}}

# Identify x-y vector components/pairs for rotation (NB: not east-west pairs!)
vector_pairs_xy = [
    ['x_wind', 'y_wind'],
    ['sea_ice_x_velocity', 'sea_ice_y_velocity'],
    ['x_sea_water_velocity', 'y_sea_water_velocity'],
    ['sea_surface_wave_stokes_drift_x_velocity',
     'sea_surface_wave_stokes_drift_y_velocity']
    ]


