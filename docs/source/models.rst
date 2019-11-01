Models
=======

.. plot::

  import matplotlib.pyplot as plt
  import numpy as np

  x = np.random.randn(1000)
  plt.hist( x, 20)
  plt.grid()
  plt.title(r'Normal: $\mu=%.2f, \sigma=%.2f$'%(x.mean(), x.std()))
  plt.show()



.. plot::

  from datetime import datetime, timedelta
  from opendrift.models.openoil import OpenOil
  from opendrift.readers import reader_netCDF_CF_generic

  o = OpenOil(loglevel=0)

  #reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')
  #reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

  reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
  reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

  o.add_reader([reader_norkyst, reader_arome])

  lon = 4.6; lat = 60.0; # Outside Bergen

  time = [reader_arome.start_time, reader_arome.start_time + timedelta(hours=30)]

  # Seed oil elements at defined position and time
  o.seed_elements(lon, lat, radius=50, number=3000, time=time,
                  wind_drift_factor=.02)

  # Adjusting some configuration
  o.set_config('processes:dispersion', True)
  o.set_config('processes:evaporation', True)
  o.set_config('processes:emulsification', True)
  o.set_config('drift:current_uncertainty', .1)
  o.set_config('drift:wind_uncertainty', 1)

  # Running model
  o.run(steps=5, time_step=1800)

  # Print and plot results
  print(o)
  o.plot(fast = True, background=['x_sea_water_velocity', 'y_sea_water_velocity'], buffer=.5)
  #o.animation(fast = True, buffer=.5)

