History
=======

2023-05-02 / Release v1.10.7
----------------------------
* CF projection info is now parsed with pyproj.CF.from_cf()
* Fixed bug in rotate_variable_dict for rotated pole projection
* netCDF generic reader now accepts Xarray Datasets in addition to filenames or URLs
* ROMS reader now accepts also time variable named 'bulk_time' and unit of days. Added uwnd,uwind,vwnd,wvind,tair,wspd to mapping variables

2023-03-29 / Release v1.10.6
----------------------------
* Added five new oils to OpenOil/ADIOS. Mapped NJORD 1997 to NJORD 2002.
* Temporary hack to let reader_netCDF_CF_generic read Zarr datasets
* Particles in air (z>0) are not stranded/deactivated when land_binary_mask==1
* Updated Thredds URL to CMEMS wave model
* Not dropping Vtransform in reader_ROMS_native when using MFDataset (wildcards). Thanks to Tianning Wu for spotting bug
* GUI: Timezone CET can be chosen, and added button to copy netCDF outfile to selected folder

2023-01-26 / Release v1.10.5
----------------------------
* Multiple improvements to the chemicaldrift model.
* Fix issue where oil type alias for 'EKOFISK BLEND 2002' did not work.
* Leeway: number of elements now equal to length of lon,lat input array (if number not given).
* Leeway: ASCII output gives small numerical difference on different platforms, presumably because of numerical errors.
* Fixing bug in get_environment, where unmasked arrays of nan did not lead to call for more readers.
* Add trajan as dependency.

2022-11-16 / Release v1.10.4
----------------------------
* Workaround in reader_netCDF_CF_generic to prevent wrong wind field from ECMWF model to be selected

2022-11-16 / Release v1.10.3
----------------------------
* Fix paths in opendrift_gui.

2022-11-16 / Release v1.10.2
----------------------------
* Optimizations to reading results files.
* ROMS reader improvements.
* ChemD: many improvements.
* Bugfixes.

2022-09-27 / Release v1.10.1
----------------------------
* Using cartopy shapes for full resolution again because of performance issues.
* Unit of oil viscosity (which is kinematic viscosity) is now consistent.
* When importing a subset in time, the number of actual active elements is now detected and used for initialization.

2022-09-26 / Release v1.10.0
----------------------------
* OpenDrift and roaring-landmask is now available as conda packages in conda-forge.
* Roaring landmask is now the only standard landmask provider. The `extent` and corners arguments
  have been removed from the global_landmask reader. They have not been in use when roaring-landmask
  was installed.
* The land shapes included with roaring-landmask is used if full resolution is used during plotting. Otherwise the cartopy provider is used.
* `Two bugs in OpenOil fixed by Giles Fearon <https://github.com/OpenDrift/opendrift/commit/78f2bd491ddc554d018e8527f97430211aafbba4>`__: in vertical mixing procedure, Temperature has wrong unit when calculating water density, and diameter was used instead of radius to calculate terminal velocity. This lead to moderate errors in vertical distribution of oil droplets: https://github.com/OpenDrift/opendrift/commit/457ed0ff263fb2cd51125cbc3df8c972e0b16fe7
* Fixed offset error in plotting of background fields on animations, which arose due to recent updates of matplotlib.
* Added fix (suggested by user lyingTree) for problem when seeding small number of elements within polygons.
* `figsize` is new optional argument to plot and animation methods (default is 11 inches).
* Possible to specify custom title for animation method.
* Oil type is now retrieved from stored netCDF files from OpenOil simulations.
* Fixed bug for readers with property `always_valid=True`
* Added boolean option show_trajectories to `plot` method.
* `reader_netCDF_CF_generic` does now only detect 1D-variables as x- and y-coordinates.
* For animated drifters, trajectory is now shown only up to current time step.
* Variables may now also be specified for `add_readers_from_list`.
* Allowing more than one drifter-dictionary to be animated, if keyword `drifter` (previously named `trajectory_dict`)  is a list instead of dict.
* New convenience method for structured readers to calculate ocean depth, area and volume within given coordinates.
* Generic netCDF reader now raises an error of file/URL is (apparently) raw ROMS output.
* ROMS native reader is now not rotating vectors with east/north in either variable or standard-name.
* Updates to ROMS native reader: standard_name_mapping may be provided by user, and mask, coordinates and angle may all be read from eventual gridfile.
* Added option to chose ensemble member in `reader_netCDF_CF_generic` (by user `mateuszmatu`).
* An experimental drift model based on the Eulerian modeling scheme has been added.
* It is now possible to combine readers using operators, e.g. to take the mean of two readers, or tune the intensity of a variable. See the `example_reader_operators.py` for an example.


2022-03-18 / Release v1.9.0
---------------------------
* Now using Cartopy >= 0.20. Cartopy < 0.20 is longer supported.
* Updated thredds URL to Barents2.5 ocean model
* ROMS native reader now detects variables having standard_name attribute
* Using more explicit exceptions internally, e.g. OutsideSpatialCoverageError, CouldNotInitializeReaderError etc.
* Added 7 Norwegian oils
* roaring_landmask (written in Rust) is now installed as default (faster landmask checking)


2022-02-28 / Release v1.8.4
---------------------------
* Fixed discarding of irrelevant readers, which was not working properly. Readers are now discarded if they do not cover simuation temporal or spatial coverage, or do not contain relevant variables
* Updating/renaming global CMEMS MERCATOR thredds URL. Removing obsoleted CMEMS reader
* Config setting drift:horizontal_diffusivity is changed from ADVANCED to BASIC, so that it is configurable from e.g. Drifty
* Fixed bug preventing export of final time step if the final time_step output is not completed
* Fixed bug in ShipDrift model: beta2 was not updated in loop, giving minor directional error
* Fixed bug in ShipDrift model: left and right directions were swapped

2022-01-31 / Release v1.8.3
---------------------------
* Removing duplicate oils in OpenOil

2022-01-31 / Release v1.8.2
---------------------------
* Re-inserted missing oil UTGARD CONDENSATE 2021, and added mapping from EKOFISK BLEND 2002 to 2000

2022-01-27 / Release v1.8.1
---------------------------
* Fixed bug in ShipDrift: erroneous direction used for wave forcing when Stokes drift was provided as forcing.
* New methods to calculate Liu-Weissberg and DARPA skillscores
* Blit is now an input parameter to animation, defaulting to False, as blitting destroys zorder (background field is always overlaid landmask)

2022-01-06 / Release v1.8.0
---------------------------
* The oil-library has been replaced with the new ADIOS database. Oils are
  retrieved from `adios.orr.noaa.gov <https://adios.orr.noaa.gov/>`_, but
  shipped with OpenDrift. They will be updated occasionally. Additional oils
  not yet included in ADIOS are also supplied with OpenDrift.
* A custom oil can be specified to OpenOil as a JSON string in the format of
  ADIOS. This means that if you want to use a new or updated oil from the ADIOS
  database, you can download it as JSON and specify it manually.
* The dependency on the oillibrary is now removed, and we should no longer have
  any conda-specific package dependencies.
* Faster writing of animations to file (mp4 and gif) using grab_frame and saving methods in matplotlib.animation writers
* New element property `current_drift_factor` (default 1) to OceanDrift and submodels - allowing to move particles with a fraction of ocean current.
* OpenOil and PlastDrift now inherits ElementType class from OceanDrift, instead of from Elements.PassiveTracer
* Fixed `bug <https://github.com/OpenDrift/opendrift/commit/7c49edaea55a65f3781363457b504c5dd86f55b2>`__ for vertical mixing with depths below 255m
* A new model :mod:`sealice <opendrift.models.sealice>` has been added, written by `Julien Moreau <https://github.com/Boorhin>`_.
* `Machine learning correction <https://opendrift.github.io/_modules/opendrift/models/oceandrift.html#OceanDrift.machine_learning_correction>`__ in OceanDrift model. Used for DARPA FFT Challenge, with machine learning data generated by Jean Rabault. Will be made avaiable for general use in future release.

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
