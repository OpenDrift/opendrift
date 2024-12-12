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

.. |ex1| image:: /gallery/images/thumb/sphx_glr_example_depth_thumb.png
   :width: 250px
   :target: /gallery/example_depth.html

.. |ex2| image:: /gallery/animations/example_oilspill_seafloor_0.gif
   :width: 250px
   :target: /gallery/example_oilspill_seafloor.html

.. |ex3| image:: /gallery/images/thumb/sphx_glr_example_oil_budget_noaa_thumb.png
   :width: 250px
   :target: /gallery/example_oil_budget_noaa.html

.. |ex4| image:: /gallery/images/thumb/sphx_glr_example_leeway_thumb.png
   :width: 250px
   :target: /gallery/example_leeway.html

+--------------------------+-------------------------+
| |ex1|                    | |ex2|                   |
|                          |                         |
| Drift at different depths| Buoyancy and turbulence |
+--------------------------+-------------------------+
| |ex3|                    | |ex4|                   |
|                          |                         |
| Oil budget               | Search and Rescue       |
+--------------------------+-------------------------+

Some key features of OpenDrift are:

* Open source (GPLv2): providing full transparency.
* Fast: optimised and vectorised Python code. Reading forcing data from disk/internet is normally the bottleneck.
* Modular: may simulate transport and fate of any kind of particles (oil, ships, persons, icebergs, micro plastics...)
* May use input forcing data (current, wind and waves) from any model, in many file formats and any map projection.
* May use backup data sources if first choice is temporarily unavailable.
* Can simulate backwards in time (specify a negative time step).
* Output is saved to CF-compliant netCDF files.
* Basic graphical user interface.
* Can use input from ensemble models.

Running a simulation can be as simple as:

.. code::

   >>> from datetime import datetime, timedelta
   >>> from opendrift.models.oceandrift import OceanDrift
   >>>
   >>> o = OceanDrift()
   >>> o.add_readers_from_list(
   >>>     ['https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])
   >>> o.seed_elements(lon=4.85, lat=60, time=datetime.now(), number=10000, radius=1000)
   >>>
   >>> o.run(duration=timedelta(hours=24))
   >>> o.animation(filename='animation.mp4')

.. image:: /gallery/animations/example_plot_0.gif

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
   performance
   tutorial
   theory/index
   choosing_a_model
   writing_a_new_model
   gallery/index
   oil_types
   interaction_with_coastline
   docker
   gui
   references
   services
   autoapi/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
* :doc:`autoapi/index`

.. *:ref:`modindex`

