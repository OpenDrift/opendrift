History
=======

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
* Conda packages at https://anaconda.org/OpenDrift/repo (see :ref:`miniconda_install`)
* Pypi packages (see :ref:`pip_install`)
* OilLibrary also ported to Python 3
* `Python 2 support dropped <https://github.com/python/devguide/pull/344>`_ (but may still work for a while)
