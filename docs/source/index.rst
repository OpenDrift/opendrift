.. toctree::
   :hidden:

   self

Introduction to OpenDrift
=====================================

OpenDrift is a software package for modeling the trajectories and fate of
objects or substances drifting in the ocean, or even in the atmosphere.

OpenDrift is open source (available on `GitHub <https://github.com/OpenDrift/opendrift/>`_), and is programmed in Python. As the software is very
generic, it is rather a "framework" than a "trajectory model" in the
traditional sense. Trajectory models for specific purposes (e.g. :py:mod:`oil
drift <opendrift.models.openoil>`, :py:mod:`search and rescue
<opendrift.models.leeway>`, :py:mod:`larvae drift
<opendrift.models.pelagicegg>` etc) may reuse all common functionality from the
:py:mod:`core model <opendrift.models.basemodel>`, and need only implement a Python Class describing the
purpose-specific processes (physics/biology etc). See
:doc:`theory/specification` and :doc:`theory/data_model` for more detailed
information.

A `journal paper about OpenDrift <https://www.geosci-model-dev.net/11/1405/2018/>`_ is published in Geoscientific Model Development. A 'paper describing the technical details of the oil spill module `OpenOil <https://www.ocean-sci.net/14/1581/2018/>`_ is published in Ocean Science.

.. image:: https://dl.dropboxusercontent.com/s/u9apyh7ci1mdowg/opendrift.gif?dl=0

Some key features of OpenDrift are:

* Open source (GPLv2): providing full transparency.
* Fast: typical simulation time is about 30 seconds for a 66 hour hour simulation of 1000 particles.
* Modular: may simulate transport and fate of any kind of particles (oil, ships, persons, icebergs)
* May use input forcing data (current, wind and waves) from any model, in many file formats and any map projection.
* May use backup data sources if first choice is temporarily unavailable.
* Can simulate backwards in time (specify a negative time step).
* Output is saved to CF-compliant netCDF files.
* Basic graphical user interface.
* Can use input from ensemble models.

.. plot::

  from datetime import timedelta
  from opendrift.models.openoil import OpenOil
  from opendrift.readers import reader_netCDF_CF_generic
  o = OpenOil()

  reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
  reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

  o.add_reader([reader_norkyst, reader_arome])

  lon = 4.6; lat = 60.0; # Outside Bergen

  time = [reader_arome.start_time, reader_arome.start_time + timedelta(hours=30)]

  # Seed oil elements at defined position and time
  o.seed_elements(lon, lat, radius=50, number=3000, time=time,
                  wind_drift_factor=.02)

  # Adjusting some configuration
  o.set_config('processes:dispersion', False)
  o.set_config('processes:evaporation', False)
  o.set_config('processes:emulsification', True)
  o.set_config('drift:current_uncertainty', .1)
  o.set_config('drift:wind_uncertainty', 1)

  # Running model
  o.run(steps=60, time_step=1800)

  # Print and plot results
  o.plot()


Once you have OpenDrift :doc:`installed <install>`, take a look at the
:doc:`tutorial` on how to get started, or check out the :doc:`gallery
<gallery/index>` for some examples. The details and physics of each model
is documented in the reference, along with specific examples for that model.
Models can be configured in nuanced ways which are important for more realistic
simulation (e.g. diffusion for :py:mod:`oil drift <opendrift.models.openoil>`
simulations). These should also be documented under the reference for the
particular model. See :py:mod:`opendrift.models` for an overview.

.. note::
   If you found OpenDrift useful for your study, please cite it as:

   * Dagestad, K.-F., Röhrs, J., Breivik, Ø., and Ådlandsvik, B.: `OpenDrift v1.0: a generic framework for trajectory modelling <https://www.geosci-model-dev.net/11/1405/2018/gmd-11-1405-2018.pdf>`_, Geosci. Model Dev., 11, 1405-1420, `<https://doi.org/10.5194/gmd-11-1405-2018>`_, 2018.

   For the oil spill module, please cite in addition to the above:

   * Röhrs, J., Dagestad, K.-F., Asbjørnsen, H., Nordam, T., Skancke, J., Jones, C. E., and Brekke, C.: `The effect of vertical mixing on the horizontal drift of oil spills <https://www.ocean-sci.net/14/1581/2018/>`_, Ocean Sci., 14, 1581-1601, `<https://doi.org/10.5194/os-14-1581-2018>`_, 2018.

Contents
===============

.. toctree::
   :maxdepth: 3
   :glob:

   history_link
   install
   tutorial
   theory/index
   choosing_a_model
   writing_a_new_model
   gallery/index
   oil_types
   interaction_with_coastline
   gui
   references
   autoapi/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
* :doc:`autoapi/index`

.. *:ref:`modindex`

