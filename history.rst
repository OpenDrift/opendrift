History
=======

Next release
---------------------------
* A new model :mod:`sealice <opendrift.models.sealice>` has been added, written by `Julien Moreau <https://github.com/Boorhin>`_.

2021-11-08 / Release v1.7.3
---------------------------
* reader_from_url is now using requests instead of urllib, fixing problem with add_readers_from_list and .netrc authentication.
* Hidden feature for ``reader_netCDF_CF_generic``: if attributes ``shift_x`` and ``shift_y`` are defined, the returned fields are shifted this many meters in the x/y (or east/north) directions
* parameter ``show_particles`` to plot() is now renamed to ``show_elements``, as for animation()
* Map bounds are now extended to cover also comparison simulations and any trajectory_dicts.
* ``skip`` and ``scale`` as input to plot() and animation() are now None, so that density and length and arrows are determined by matplotlib/quiver, unless overridden by user.
* New method (``distance_between_trajectories``) to calculate distances between two trajectories, position by position.
* Updates to ``ChemicalDrift`` model

2021-10-27 / Release v1.7.2
---------------------------
* Fix bugs in selafin reader.
* Several improvements to the SCHISM reader.
* Add method for tuning windrift factor from observed drift.
* Add method to retrieve environment variables (from given readers) along a given trajectory (e.g. a drifter).
* Improved dateline handling in readers.
* Fix dateline bug in landmask.
* ``reader_netCDF_CF_generic``: if x, and y-coordinates are integer sequences, these are not anymore interpreted as projection coordinates.
* ``reader_netCDF_CF_generic``: taking calendar name into acount when decoding time.
* Leeway model: max_speed is increased to 5 m/s, avoiding obtaining too small data-blocks readers.
* Leeway model ASCII export: if all elements are deactivated, write previous mean position, instead of NaN.
* Improved Xarray-postprocessing (based on `opendrift.open_xarray`), as demonstrated in ``example_river_runoff.py``. Aotomatic ``analysis_file`` is omitted.
* Fixed problem related to mutating dictionary of readers when discarding.
* Added ``mixed_layer_depth`` (default 50m) as environment variable of OceanDrift (and subclasses). This is used if ``Sundby`` or ``Large`` parameterizations of vertical turbulence is activated. A new config setting defines background diffusivity (default: ``1.2e-5 m2-s``)
* ``origin_marker_name`` can now be specified when seeding, and is stored as attributes ``flag_meanings`` to output variable ``origin_marker``.
* Quiver plots are now centered on pixels/grid, instead of at corner.

2021-09-01 / Release v1.7.1
---------------------------
* Using OilLibrary v4+noaa1.1.3

2021-08-30 / Release v1.7.0
---------------------------
* New method ``reader.shift_start_time(start_time)`` to shift time coverage of reader
* Density arrays calculated with method "get_density" for files opened with `open_xarray` can now be weighted with any property, or a user provided array. `origin_marker is now a dimension of the arrays stored in analysis netCDF file. Made new method `get_density_timeseries`
* ROMS native reader now accepts datasets where lon and lat are 1-dimensional
* Fixed bug related to extrapolating 3D data to seafloor
* Fixed bug with interpolation where latitude/y-coordinate is decreasing and not increasing (flipped upside down). Also fixed small inaccuracy of structured interpolation.
* Fixed horizontal diffusion for backwards simulations
* Enable the use of `roaring-landmask <https://github.com/gauteh/roaring-landmask>`_ as landmask reader, if installed.
* Add Telemac / Selafin reader (requires telemac python scripts).

2021-05-03 / Release v1.6.0
-----------------------------
* Reader environment mappings (deriving variables from others) can be activated with >>> o.activate_environment_mapping(<mapping_name>). Method to derive wind components from ``wind_speed`` and ``wind_from_direction / wind_to_direction`` is activated by default.
* New unstructured reader for SHYFEM model output
* ``animation`` and ``animation_profile`` methods may now use legend instead of colorbar for element properties
* Arguments ``color`` to ``animation()`` and ``linecolor`` to ``plot()`` can now be arrays of length equal to the number of elements.
* Improved mechanism for drifter/trajectory overlay on animations, as illustraded by :doc:`example_current_from_drifter <gallery/example_current_from_drifter>`
* Several improvements to module ChemicalDrift
* For PlastDrift model, config ``drift:vertical_mixing=False`` still gave vertical entrainment for ``mixingmodel=analytical``, but this is now changed. Sundby83 is now default model for vertical diffusivity in PlastDrift (was Large1994)
* Increased valid range of current velocity components from 10 m/s to 15 m/s
* Rotated pole projection (ob_tran) is now parsed from CF attributes by reader_netCDF_CF_generic.
* Leeway jibing probability is calculated with exponential, giving more precise results for larger time steps. Generic arguments are removed from Leeway seeding method.
* lon, lat are now positional arguments also in Leeway.seed_elements method. Leeway.seed_from_shapefile did nor work before this fix.
* Config option ``drift:lift_to_seafloor`` is replaced by ``general:seafloor_action``, analoguos to ``general:coastline_action``.
  Available options are ``none``, ``deactivate``, ``lift_to_seafloor`` as well as new option ``previous`` - moving elements back to previous position.
* New method ``get_trajectory_lengths`` to calculate length and speeds along trajectories
* Basemodel class does not anymore have a projection, internal coordinates are now always lon, lat
* Color of ocean and landmask may now be overridden in plot- and animation methods with new input variables ``land_color`` and ``ocean_color``. A new input dictionary ``text`` allows map annotations.
* opendrift-landmask-data only loads mask once for each python process, reducing memory usage and improves performance where you run opendrift multiple times in the same script and process.

2021-02-15 / Release v1.5.6
-----------------------------
* New parallelisation of lonlat2xy for unprojected readers. The flag ``<reader>.multiprocessing_fail`` is replaced with ``<reader>.__parallel_fail__``
* plot_property() can now save figure to file if filename is provided
* netCDF attribute seed_geojson is now a GeoJSON FeatureCollection.
* reader_netCDF_CF_generic does not anymore read 2D lon/lat variables if 1D x/y variables are detected, giving much faster initialisation.
* General replacement of ``np.float`` and ``np.int`` with either ``float``, ``int`` or ``np.float32/64`` and ``np.int32/64``. np.float and np.int are deprecated in numpy 1.20.
* Fixed bug occuring when interpolating environment_profiles in time, and the number of vertical layers in the ocean-model-block is larger at time1 than at time2

2021-01-26 / Release v1.5.5
---------------------------
* New module LarvalFish, for fish eggs hatching into larvae with swimming behaviour
* Sundby83 parameterisation of vertical diffusivity is now set to 0 below mixed layer depth (default 50m)
* Deprecating seed argument `oiltype` in favor of `oil_type` in OpenOil. Warning is issued, but later this will become an error
* Fixed problem with convolution of reader fields
* Fixed newly introduced bug with Leeway ascii output file
* Cleaned up some metadata output, and seeding arguments are written as list of GeoJSON strings to attribute `seed_geojson`

2021-01-18 / Release v1.5.4
---------------------------
* seed_cone also accepts time as list with single element
* Min/max values are checked/masked also for ensemble data
* reader_netCDF_CF_generic now detects lon/lat arrays also if their variable name equals lon/lat or longitude/latitude

2021-01-15 / Release v1.5.3
---------------------------
* Fixed bug related to derived_variables (e.g. calculating x_wind, y_wind from windspeed, winddirection)

2021-01-14 / Release v1.5.2
---------------------------
* Fixed problem with double or missing logging output
* ShipDrift model now gives warning and not error if input parameter are outside bounds, and parameters are clipped to boundary values
* Fixed problem with multiprocessing/parallelization of lonlat2xy for unprojected readers

2021-01-05 / Release v1.5.1
---------------------------
* OilLibrary updated to version 1.1.3. Slightly different weathering results, and * is removed from oil names starting with GENERIC

2021-01-04 / Release v1.5.0
---------------------------
* Major restructuring of Basereader class. Readers now are sublasses of Structured, Unstructured or Continuous.
* Built in GUI is improved with posibillity to adjust all config settings.
* Some Leeway parameters are renamed from camelCase to camel_case, including: ``jibeProbability`` -> ``jibe_probability`` and ``objectType`` -> ``object_type``
* Renamed config setting ``drift:scheme`` -> ``drift:advection_scheme``

2020-11-01 / Release v1.4.2
---------------------------

* Fixed bug in v1.4.1 that OpenOil and SedimentDrift had fallback_value of 0 for `land_binary_mask`, this shall be `None`.

2020-10-31 / Release v1.4.1
---------------------------

* Built in GUI is improved with docstrings and less hardcoding, based on new config mechanism, including a new bool setting ``seed:seafloor``.
* ``model.required_variables`` is now a dictionary, which also includes the earlier ``fallback_values``, ``desired_variables`` and ``required_profiles``. Instead of providing fallback values directly in a dictionary, these shall now be provided through the config mechanism: ``o.set_config('environment:fallback:<variable>', <value>)``. Correspondingly, config setting ``environment:constant:<variable>`` may be used to specify constant values for the same variables (overriding any other readers).
* `seed_elements <https://opendrift.github.io/autoapi/opendrift/models/basemodel/index.html#opendrift.models.basemodel.OpenDriftSimulation.seed_elements>`_ is simplified, by factoring out a new method `seed_cone <https://opendrift.github.io/autoapi/opendrift/models/basemodel/index.html#opendrift.models.basemodel.OpenDriftSimulation.seed_cone>`_

2020-10-27 / Release v1.4.0
---------------------------

* New internal config mechanism, and configobj package is no longer needed. The user API (``get_config()``, ``set_config()``) is unchanged, but model developers must use the `new mechanism <https://opendrift.github.io/autoapi/opendrift/models/basemodel/index.html#opendrift.models.basemodel.OpenDriftSimulation._add_config>`_ to add configuration settings.
* Added new reader for static 2D fields (``reader_constant_2d.py``)
* Xarray, Dask and Xhistogram are new requirements. New method ``opendrift.open_xarray`` to open an output netCDF file lazily, with possibility to e.g. calculate density arrays/plots from datasets to large to fit in memory.
* New model chemicaldrift

2020-10-15 / Release v1.3.3
---------------------------

* New seed method ``seed_repeated_segment()``
* New method ``animate_vertical_distribution()``
* Vertical mixing scheme is greatly simplified, and should be faster for large number of elements.
* Vertical mixing is now disabled by default in OceanDrift, but enabled in all submodules (PelagicEggDrift, SedimentDrift, RadionuclideDrift, OpenOil)
* Vertical diffusivity option `zero` is replaced with ``constant``, which means using the fallback value.
* New config setting ``drift:horizontal_diffusivity``, providing time-step independent diffusion, in contrast to ``drift:current_uncertainty`` and ``drift:wind_uncertainty``
* Readers may be initialised from a JSON string, where `reader` is name of reader module, and other parameters are forwarded to reader constructor, e.g.: `{"reader": "reader_cmems", "dataset": "global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh"}`
* CMEMS reader now obtains username/password from .netrc instead of environment variables. CMEMS-motuclient is added to environment.yml
* CMEMS reader now takes dataset name and not product name as input, and it is possible to provide variable mapping.
* NOAA ADIOS is now default (and only) option for oil weathering, as the "built in" oil weathering module ("basic") is removed.
* GUI is generalised, to be usable for any modules. This includes taking default seed options from `config:seed:` (e.g. m3_per_hour for OpenOil)

2020-08-21 / Release v1.3.2
---------------------------

* Fixed ``vmax`` value for animations with density array
* Fixed animation marker color for deactivated elements
* Introduced mechanism to store environment variables from previous time step
* New element property ``moving``, giving possibility to temporarily freeze elements, e.g. used for sedimentation and resuspension in SedimentDrift module
* Improved robustness using Xarray in netCDF-readers. Xarray is still optional dependency, but is now tested on Travis
* nc-time-axis is new dependency, providing support for cftime axis in matplotlib

2020-07-03 / Release v1.3.1
---------------------------

* NOAA oil weathering model is now default choice in OpenOil
* Bugfix in reader_netCDF_CF_generic for particles with negative longitudes combined with global datasets with longitudes from 0-360
* Added module ``SedimentDrift``
* Removed two options from OpenOil, with corresponding config parameters:

  * Tkalich(2002) entrainment rate

  * Exponential droplet size distribution

* Renamed two config settings:

  * ``processes:turbulentmixing`` -> ``drift:vertical_mixing``

  * ``processes:verticaladvection``-> ``drift:vertical_advection``

2020-06-24 / Release v1.3.0
------------------------------
* OceanDrift3D and OpenDrift3D have been merged into OceanDrift, and OpenOil3D has been merged into OpenOil. Thus OpenOil and OceanDrift are now 3D modules, but can still be configured for 2D drift.

2020-06-12 / Release v1.2.3
---------------------------

* Seed from shapefile: GDAL (ogr / osr) changed coordinate order, updates dependencies and call.

2020-06-08 / Release v1.2.2
---------------------------

* `Victor de Aguiar <https://github.com/vic1309>`_: :mod:`Oil drift in sea ice <opendrift.models.openoil>` following Nordam et al., 2019, doi:10.1016/j.marpolbul.2019.01.019 (Sponsored by the Fram Centre in Tromsø, through the MIKON/OSMICO project).
* OpenBerg module available from the GUI.
* A generic shape reader for landmasks (use polygons directly or convenience method using shp files).
* Drop rasterio dependency and include some significant thread-safety fixes for landmask-data.

2020-05-14 / Release v1.2.1
---------------------------

* Specifying a positive time step with a negative duration is now an error. Duration should
  always be specified positive.

2020-01-08 / Release v1.2.0
---------------------------

* Basemap reader and basemap plotting removed
* Minor improvements and bug fixes
* Example scripts are now available in online :doc:`gallery <gallery/index>`
* Only a single conda environment (named "opendrift"). Fresh :doc:`installation <install>` is recommended.

2019-11-27 / Release v1.1.1
---------------------------

* Cartopy is used for plotting (with fast option only using raster, see :meth:`opendrift.models.basemap.plot`)
* GSHHS full is used for a dedicated landmask reader (using full resolution always, possibly to :mod:`only use mask <opendrift.readers.reader_global_landmask>` for faster simulations)
* New documentation at https://opendrift.github.io
* Conda packages at https://anaconda.org/OpenDrift/repo
* Pypi packages
* OilLibrary also ported to Python 3
* `Python 2 support dropped <https://github.com/python/devguide/pull/344>`_ (but may still work for a while)
