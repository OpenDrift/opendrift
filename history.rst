History
=======

Next release
--------------------------

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
* New documentation at http://opendrift.github.io
* Conda packages at https://anaconda.org/OpenDrift/repo (see :ref:`miniconda_install`)
* Pypi packages (see :ref:`pip_install`)
* OilLibrary also ported to Python 3
* `Python 2 support dropped <https://github.com/python/devguide/pull/344>`_ (but may still work for a while)
